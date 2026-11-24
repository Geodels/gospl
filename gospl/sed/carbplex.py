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

    def _carbSlope(self, h):
        """Seafloor gradient magnitude per node: the steepest drop or rise to
        a neighbour (local; owned rows see their full neighbourhood)."""
        ngb = self.FVmesh_ngbID
        nb = np.where(ngb >= 0, ngb, np.arange(self.lpoints)[:, None])
        d = np.linalg.norm(self.lcoords[nb] - self.lcoords[:, None, :], axis=2)
        g = np.divide(np.abs(h[nb] - h[:, None]), d, out=np.zeros_like(d), where=d > 0.0)
        return g.max(axis=1)

    def _carbStatic(self, h):
        """Per-factory factor that does not change within the step: habitat
        map times the ``slope`` and ``sedrate`` responses."""
        out = []
        slope = None
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
