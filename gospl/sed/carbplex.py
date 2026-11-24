import numpy as np
import petsc4py
from mpi4py import MPI
from time import process_time

from gospl.tools.constants import CARB_SUBSTEP_GROWTH

MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()


class CarbMesh(object):
    """
    Benthic carbonate production (opt-in ``carbonate:`` block,
    docs/DESIGN_CARBONATES.md, phase C1).

    Each factory produces below sea level at the Bosscher & Schlager (1992)
    rate

    .. math::
        P_f(w) = g_f \\tanh\\left(\\frac{I_0 e^{-k_f w}}{I_{k,f}}\\right)
                 M_f \\prod_j R_{fj}

    with ``w`` the water depth, ``M_f`` an optional habitat map and ``R_fj``
    response tables on the controls ``depth``, ``slope`` and ``sedrate``
    (the clastic marine deposition rate of the previous step). Production is
    integrated per node over the step with the water depth it creates
    (``dh/dt = sum_f P_f(sealevel - h)``), so a node never builds above
    ``sealevel - min_depth`` and keep-up, catch-up and give-up follow from
    the law. The ``framework`` part of each factory is built in place
    (lithology ``reef``); the loose part (``carb_sand``, ``carb_mud``) joins
    the marine routing of the same step (``seaplex.seaChange``).

    Everything here is pointwise (no collective apart from the budget
    reduction), so production is partition-exact by construction.
    """

    def __init__(self):
        """
        Allocate the carbonate state when ``carbOn``: production rate output,
        the loose sediment handed to the marine step, the previous clastic
        deposition rate and the habitat maps (local slices).
        """
        self.carbProd = None
        self.carbBudget = None
        self._carbDet = None
        # cached screened-Poisson operator and solver of the open-water index
        # (registered in destroy_DMPlex)
        self._carbMat = None
        self._ksp_carb = None
        self.pelDep = None
        self.pelBudget = None
        if getattr(self, "pelOn", False):
            self._initPelagic()
        if not getattr(self, "carbOn", False):
            return
        if not getattr(self, "seaDepo", False) or getattr(self, "nodep", False):
            raise ValueError(
                "carbonate: needs marine deposition (domain: seadepo: True, "
                "nodep: False): the loose carbonate sediment is spread by the "
                "marine routing."
            )
        self.carbProd = np.zeros(self.lpoints, dtype=np.float64)
        self._carbSedRate = np.zeros(self.lpoints, dtype=np.float64)
        self._carbHab = []
        for f in self.carbFact:
            if f["habitat"] is None:
                self._carbHab.append(None)
                continue
            d = np.load(str(f["habitat"][0]) + ".npz")
            m = np.asarray(d[f["habitat"][1]], dtype=np.float64)
            if m.shape[0] != self.mpoints:
                raise ValueError(
                    "carbonate: habitat map %s has %d values, the mesh %d."
                    % (f["habitat"], m.shape[0], self.mpoints))
            self._carbHab.append(np.clip(m[self.locIDs], 0.0, None))
        self._carbProduced = 0.0
        self._carbFramework = 0.0
        self._carbLoose = 0.0
        return

    def _initPelagic(self):
        """Allocate the pelagic rain state: deposit rate output, the per-node
        carry-over of sub-floor deposits, budgets and the event maps (local
        slices, m/yr)."""
        self.pelDep = np.zeros(self.lpoints, dtype=np.float64)
        self._pelCarry = np.zeros((self.lpoints, 2), dtype=np.float64)
        self._pelRain = 0.0
        self._pelPreserved = 0.0
        self._pelDissolved = 0.0
        self._pelClay = 0.0
        self._pelMaps = []
        for e in self.pelEvents:
            m = {}
            for key in ("carb", "clay"):
                v = e[key]
                if isinstance(v, tuple):
                    d = np.load(v[0] + ".npz")
                    arr = np.asarray(d[v[1]], dtype=np.float64)
                    if arr.shape[0] != self.mpoints:
                        raise ValueError(
                            "pelagic: map %s has %d values, the mesh %d."
                            % (list(v), arr.shape[0], self.mpoints))
                    m[key] = np.clip(arr[self.locIDs], 0.0, None) * 1.0e-6
                else:
                    m[key] = float(v)
            self._pelMaps.append(m)

    def _pelEvent(self):
        """Index of the pelagic event in force at ``tNow`` (the last one that
        has started; the first one before it starts)."""
        k = 0
        for i, e in enumerate(self.pelEvents):
            if e["start"] <= self.tNow + 1.0e-6:
                k = i
        return k

    def _pelPreservedFraction(self, depth, lys, ccd):
        """Fraction of the pelagic carbonate rain preserved at water depth
        ``depth``: 1 above the lysocline, 0 below the CCD, linear between."""
        if ccd <= lys:
            return np.where(depth <= lys, 1.0, 0.0)
        return np.clip((ccd - depth) / (ccd - lys), 0.0, 1.0)

    def pelagicRain(self):
        """
        Settle the pelagic carbonate and the hemipelagic clay of the step on
        the marine nodes (``lFill`` at or below sea level).

        The carbonate rain (event rate, uniform or map) is preserved by the
        fraction set by the water depth (``_pelPreservedFraction``), the rest
        dissolves (budget only); the clay settles entirely. The deposit never
        builds above sea level. Both are laid into the stratigraphy (``carb_mud``
        and ``fine``, so ``depoFineFrac = 1``) and raise the elevation. A
        node's deposit thinner than the stratigraphic floor of ``deposeStrat``
        (0.1 mm) is carried over to the next steps (``_pelCarry``) instead of
        being lost.

        Scratch Vecs: ``self.tmp`` / ``self.tmpL`` (the deposit on exit).
        Collective (halo sync, budget).
        """
        ev = self._pelMaps[self._pelEvent()]
        e = self.pelEvents[self._pelEvent()]
        hl = self.hLocal.getArray().copy()
        depth = self.sealevel - hl
        marine = (self.lFill <= self.sealevel) & (depth > 0.0)
        rain = np.where(marine, ev["carb"] * self.dt, 0.0) * np.ones(self.lpoints)
        keep = self._pelPreservedFraction(depth, e["lys"], e["ccd"])
        carb0 = rain * keep
        dissolved = rain - carb0
        clay0 = np.where(marine, ev["clay"] * self.dt, 0.0) * np.ones(self.lpoints)
        # never above sea level: what does not fit is bypassed (budget)
        tot = carb0 + clay0
        scale = np.divide(np.minimum(tot, np.maximum(depth, 0.0)), tot,
                          out=np.zeros_like(tot), where=tot > 0.0)
        carb = carb0 * scale
        clay = clay0 * scale
        bypass = tot - carb - clay

        # carry-over of sub-floor deposits
        self._pelCarry[:, 0] += carb
        self._pelCarry[:, 1] += clay
        tot = self._pelCarry.sum(axis=1)
        lay = tot >= 1.0e-4
        c = np.where(lay, self._pelCarry[:, 0], 0.0)
        k = np.where(lay, self._pelCarry[:, 1], 0.0)
        self._pelCarry[lay] = 0.0
        depo = c + k
        self.pelDep = depo / self.dt

        self.tmpL.setArray(depo)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self.cumED.axpy(1.0, self.tmp)
        self.hGlobal.axpy(1.0, self.tmp)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal)
        self.dm.globalToLocal(self.hGlobal, self.hLocal)
        self.tmpL.setArray(depo / self.dt)
        self.EbLocal.axpy(1.0, self.tmpL)
        if self.stratNb > 0:
            ff = self.depoFineFrac.copy()
            sh = None if self.litShare is None else self.litShare.copy()
            self.depoFineFrac[:] = 1.0
            if self.litShare is not None:
                self.litShare[:] = 0.0
                cm = self.carbLit.get("carb_mud") if hasattr(self, "carbLit") else None
                if cm is not None:
                    self.litShare[cm] = np.divide(c, depo, out=np.zeros_like(depo),
                                                  where=depo > 0.0)
                self.litShare[0] = np.divide(k, depo, out=np.zeros_like(depo),
                                             where=depo > 0.0)
            self.tmpL.setArray(depo)
            self.dm.localToGlobal(self.tmpL, self.tmp)
            self.deposeStrat()
            self.depoFineFrac = ff
            if sh is not None:
                self.litShare = sh

        own = self.inIDs == 1
        a = self.larea
        b = MPI.COMM_WORLD.allreduce(
            np.array([(rain * a)[own].sum(), (carb * a)[own].sum(),
                      (dissolved * a)[own].sum(), (clay * a)[own].sum(),
                      (bypass * a)[own].sum()]), op=MPI.SUM)
        # (carbonate rain, carbonate settled, dissolved, clay settled,
        #  bypassed above sea level), m3 over the step
        self.pelBudget = b
        self._pelRain += b[0]
        self._pelPreserved += b[1]
        self._pelDissolved += b[2]
        self._pelClay += b[3]
        self._pelBypass = getattr(self, "_pelBypass", 0.0) + b[4]
        if MPIrank == 0 and self.verbose:
            print(
                "[pelagic] carbonate rain %0.4g km3 (preserved %0.4g, dissolved "
                "%0.4g), clay %0.4g km3" % (b[0] * 1e-9, b[1] * 1e-9, b[2] * 1e-9,
                                            b[3] * 1e-9),
                flush=True,
            )
        return

    def _carbSlope(self, h):
        """Seafloor gradient magnitude per node: the steepest drop or rise to
        a neighbour (local; owned rows see their full neighbourhood)."""
        ngb = self.FVmesh_ngbID
        nb = np.where(ngb >= 0, ngb, np.arange(self.lpoints)[:, None])
        d = np.linalg.norm(self.lcoords[nb] - self.lcoords[:, None, :], axis=2)
        g = np.divide(np.abs(h[nb] - h[:, None]), d, out=np.zeros_like(d), where=d > 0.0)
        return g.max(axis=1)

    def _carbBuildOp(self):
        """Build and cache the screened-Poisson operator ``I - L^2 lap`` and
        its solver on first use (collective)."""
        if self._carbMat is not None:
            return
        kd = np.full(self.lpoints, self.carbExpLength ** 2, dtype=np.float64)
        self._carbMat = self._buildDiffMat(kd)
        ksp = petsc4py.PETSc.KSP().create(petsc4py.PETSc.COMM_WORLD)
        ksp.setType("fgmres")
        ksp.getPC().setType("bjacobi")
        ksp.setTolerances(rtol=1.0e-10, atol=1.0e-12, max_it=2000)
        ksp.setOptionsPrefix("carbsmooth_")
        petsc4py.PETSc.Options()["carbsmooth_sub_pc_factor_shift_type"] = "nonzero"
        ksp.setFromOptions()
        ksp.setOperators(self._carbMat, self._carbMat)
        self._ksp_carb = ksp

    def _carbScreened(self, rhs):
        """
        Solve ``(I - L^2 lap) x = rhs`` with ``L = carbExpLength`` (the
        finite-volume Laplacian of ``sethillslopecoeff`` with ``Kd = L^2``;
        open edges of a flat mesh are identity rows). Operator and solver are
        built once and cached (``_carbMat``, ``_ksp_carb``, prefix
        ``carbsmooth_``: fgmres + block Jacobi, which converges when L spans
        many cells). Collective.

        :arg rhs: local right-hand side
        :return: local solution
        """
        self._carbBuildOp()
        self.tmpL.setArray(rhs)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self.tmp1.set(0.0)
        self._ksp_carb.solve(self.tmp, self.tmp1)
        if self._ksp_carb.getConvergedReason() < 0 and MPIrank == 0:
            print("[carbonate] open-water index solve did not converge (%d)"
                  % self._ksp_carb.getConvergedReason(), flush=True)
        self.dm.globalToLocal(self.tmp1, self.tmpL)
        return self.tmpL.getArray().copy()

    def _carbExposure(self, h):
        """
        Open-water index ``E = h - s`` (m), ``s`` the elevation smoothed over
        ``carbExpLength`` (``(I - L^2 lap) s = h``): positive on platform
        margins, rims and highs, negative in lagoons and behind barriers. It
        is solved for ``E`` itself, ``(I - L^2 lap) E = A h - h``, so the
        solver tolerance applies to the relief, not to the elevation. For a
        parabolic ridge ``h = -b x^2`` it is exactly ``2 b L^2``. Collective.
        Scratch: ``tmp``, ``tmpL``, ``tmp1``.
        """
        self._carbBuildOp()
        self.tmpL.setArray(h)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self._carbMat.mult(self.tmp, self.tmp1)
        self.tmp1.axpy(-1.0, self.tmp)                        # A h - h
        self.dm.globalToLocal(self.tmp1, self.tmpL)
        return self._carbScreened(self.tmpL.getArray().copy())

    def _carbStatic(self, h):
        """Per-factory factor that does not change within the step: habitat
        map times the ``slope``, ``sedrate`` and ``exposure`` responses
        (the open-water index is solved once per step, only when used)."""
        out = []
        slope = None
        expo = None
        if any(ctl == "exposure" for f in self.carbFact for ctl, _, _ in f["resp"]):
            expo = self._carbExposure(h)
        self.carbExpo = expo
        for k, f in enumerate(self.carbFact):
            m = np.ones(self.lpoints, dtype=np.float64)
            if self._carbHab[k] is not None:
                m = m * self._carbHab[k]
            for ctl, x, fv in f["resp"]:
                if ctl == "slope":
                    if slope is None:
                        slope = self._carbSlope(h)
                    m = m * np.interp(slope, x, fv)
                elif ctl == "sedrate":
                    m = m * np.interp(self._carbSedRate, x, fv)
                elif ctl == "exposure":
                    m = m * np.interp(expo, x, fv)
            out.append(m)
        return out

    def _carbRates(self, w, static):
        """Production rate of every factory (m/yr) at water depth ``w``."""
        rates = []
        on = w > self.carbMinDepth
        for f, m in zip(self.carbFact, static):
            p = f["g"] * np.tanh(self.carbI0 * np.exp(-f["k"] * np.maximum(w, 0.0))
                                 / f["Ik"]) * m
            for ctl, x, fv in f["resp"]:
                if ctl == "depth":
                    p = p * np.interp(w, x, fv)
            rates.append(np.where(on, p, 0.0))
        return rates

    def _carbSubsteps(self):
        """Sub-steps of the (second-order) production integration: about
        ``CARB_SUBSTEP_GROWTH`` of growth each at the largest possible rate. Depends only on the inputs, so it
        is the same on every rank."""
        gmax = sum(f["g"] for f in self.carbFact)
        return int(min(max(np.ceil(gmax * self.dt / CARB_SUBSTEP_GROWTH), 1), 1000))

    def carbProduction(self):
        """
        Produce carbonate over the step on the marine nodes (``lFill`` at or
        below sea level), build the framework in place and keep the loose
        sediment for the marine step.

        The framework raises the elevation (``hGlobal``, ``cumED``,
        ``EbLocal``) and is laid into the stratigraphy as lithology ``reef``;
        the loose volumes are stored in ``_carbDet`` (total), ``_carbMud``
        (fine-like) and ``_carbL`` (per lithology), read by ``seaChange``.

        Scratch Vecs: ``self.tmp`` / ``self.tmpL`` (the framework thickness on
        exit, as ``deposeStrat`` needs it). Collective (halo sync, budget).
        """
        t0 = process_time()
        hl = self.hLocal.getArray().copy()
        marine = self.lFill <= self.sealevel
        static = self._carbStatic(hl)
        nsub = self._carbSubsteps()
        dts = self.dt / nsub
        h = hl.copy()
        prod = [np.zeros(self.lpoints) for _ in self.carbFact]
        for _ in range(nsub):
            # Heun (second order): rates at the start and at the Euler end
            # point of the sub-step, averaged
            w = self.sealevel - h
            r1 = self._carbRates(w, static)
            t1 = r1[0].copy()
            for r in r1[1:]:
                t1 += r
            r2 = self._carbRates(w - t1 * dts, static)
            rates = [0.5 * (a + b) for a, b in zip(r1, r2)]
            tot = rates[0].copy()
            for r in rates[1:]:
                tot += r
            dz = tot * dts
            room = np.maximum(w - self.carbMinDepth, 0.0)
            scale = np.divide(np.minimum(dz, room), dz, out=np.zeros_like(dz),
                              where=dz > 0.0)
            scale[~marine] = 0.0
            for k, r in enumerate(rates):
                prod[k] += r * dts * scale
            h += dz * scale

        total = np.zeros(self.lpoints)
        frame = np.zeros(self.lpoints)
        sand = np.zeros(self.lpoints)
        for f, p in zip(self.carbFact, prod):
            total += p
            frame += f["phi"] * p
            sand += (1.0 - f["phi"]) * f["sand"] * p
        # Framework thinner than the stratigraphic floor of deposeStrat
        # (0.1 mm) would raise the elevation without being recorded: it joins
        # the loose sand instead (conservative).
        thin = frame < 1.0e-4
        sand[thin] += frame[thin]
        frame[thin] = 0.0
        loose = total - frame
        mud = np.maximum(loose - sand, 0.0)
        self.carbProd = total / self.dt

        # The loose sediment for the marine routing (deposit volumes).
        self._carbDet = loose * self.larea
        self._carbMud = mud * self.larea
        self._carbL = np.zeros((self.lpoints, self.litNb))
        self._carbL[:, self.carbLit["carb_sand"]] = sand * self.larea
        self._carbL[:, self.carbLit["carb_mud"]] = mud * self.larea

        # Framework: built in place, recorded as `reef`.
        if "reef" in self.carbLit:
            self.tmpL.setArray(frame)
            self.dm.localToGlobal(self.tmpL, self.tmp)
            self.cumED.axpy(1.0, self.tmp)
            self.hGlobal.axpy(1.0, self.tmp)
            self.dm.globalToLocal(self.cumED, self.cumEDLocal)
            self.dm.globalToLocal(self.hGlobal, self.hLocal)
            self.tmpL.setArray(frame / self.dt)
            self.EbLocal.axpy(1.0, self.tmpL)
            ff, sh = self.depoFineFrac.copy(), self.litShare.copy()
            self.depoFineFrac[:] = 0.0
            self.litShare[:] = 0.0
            self.litShare[self.carbLit["reef"]][:] = 1.0
            self.tmpL.setArray(frame)
            self.dm.localToGlobal(self.tmpL, self.tmp)
            self.deposeStrat()
            self.depoFineFrac, self.litShare = ff, sh

        own = self.inIDs == 1
        b = MPI.COMM_WORLD.allreduce(
            np.array([(total * self.larea)[own].sum(), (frame * self.larea)[own].sum(),
                      (loose * self.larea)[own].sum()]), op=MPI.SUM)
        self.carbBudget = b
        self._carbProduced += b[0]
        self._carbFramework += b[1]
        self._carbLoose += b[2]
        if MPIrank == 0 and self.verbose:
            print(
                "[carbonate] produced %0.4g km3 (framework %0.4g, loose %0.4g) "
                "(%0.02f seconds)" % (b[0] * 1e-9, b[1] * 1e-9, b[2] * 1e-9,
                                      process_time() - t0),
                flush=True,
            )
        return

    def carbClasticRate(self):
        """
        After ``seaChange``: the clastic (non-carbonate) part of the marine
        deposit of this step as a rate (m/yr), the ``sedrate`` control of the
        next step's production (turbidity and burial). Reads the deposit
        laid into the stratigraphy from ``self.tmp`` and its composition
        (``depoFineFrac``, ``litShare``). Rank-local apart from the halo read.
        """
        self.dm.globalToLocal(self.tmp, self.tmpL)
        dep = np.maximum(self.tmpL.getArray(), 0.0)
        sh = self._groupShares()
        carb = np.zeros(self.lpoints)
        for name, l in self.carbLit.items():
            if name != "reef":
                carb += self._litDepoFrac(l, self.depoFineFrac, sh)
        self._carbSedRate = dep * np.clip(1.0 - carb, 0.0, 1.0) / self.dt
        return
