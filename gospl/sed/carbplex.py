import numpy as np
import petsc4py
from mpi4py import MPI
from time import process_time

from gospl.tools.constants import CARB_SUBSTEP_GROWTH, DEPOSIT_FLOOR, ICE_COVER_MIN

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
        # clay decay from the coast: operator rebuilt each step (the coast
        # moves), solver cached; both in destroy_DMPlex
        self._clayMat = None
        self._ksp_clay = None
        self._sstMaps = None
        self._carbSupplyScale = 1.0
        self._carbProdCache = None
        self.carbSupplyBudget = None
        self._plumeMat = None
        self._ksp_plume = None
        if getattr(self, "sstOn", False):
            self._initSST()
        if getattr(self, "carbSupply", None) is not None or any(
                f.get("solute") for f in getattr(self, "carbFact", [])):
            self._initGeoLink()
        if getattr(self, "pelOn", False):
            self._initPelagic()
        self.karstDiss = None
        self.karstBudget = None
        self._karstSolute = None
        if getattr(self, "karstOn", False):
            self._initKarst()
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

    def _initGeoLink(self):
        """Check the geochemistry coupling (phase C6): the species named by
        ``carbonate_supply`` and by the ``solute`` responses must be tracers of
        the groundwater geochemistry, with river routing and (for the supply)
        the marine reservoir."""
        names = list(getattr(self, "gwGeoName", []) or [])
        if not getattr(self, "gwRiverLoad", False) or not names:
            raise ValueError(
                "carbonate: the solute response and carbonate_supply need the "
                "groundwater geochemistry with `river_load: True`.")
        want = [f["solute"][0] for f in getattr(self, "carbFact", []) if f.get("solute")]
        if getattr(self, "carbSupply", None) is not None:
            if not getattr(self, "gwMarineCoupling", False):
                raise ValueError(
                    "carbonate_supply: needs `marine_coupling: True` in the "
                    "groundwater geochemistry (the marine reservoir).")
            want.append(self.carbSupply["species"])
        for sp in want:
            if sp not in names:
                raise ValueError(
                    "carbonate: geochemistry species '%s' not found (species: %s)."
                    % (sp, ", ".join(names)))
        if getattr(self, "carbSupply", None) is not None:
            self._carbSupplyK = names.index(self.carbSupply["species"])

    def _carbPlume(self, species, L):
        """
        Coastal plume of a geochemistry species (``solute`` response, phase
        C6): the river dissolved load of the species delivered at the coast
        this step (``riverSoluteSp`` at the coastal exits, per unit area)
        spread offshore by ``(I - L^2 lap) c = q``, land as identity rows
        (holding their own input). Units: solute per m2 per yr. The operator is
        rebuilt each call (the coast moves); cached solver ``_ksp_plume``.
        Collective. Scratch: ``tmp``, ``tmpL``, ``tmp1``.
        """
        k = list(self.gwGeoName).index(species)
        exits = np.zeros(self.lpoints, dtype=bool)
        mi = getattr(self, "marineSoluteInput", None)
        if mi is not None:
            exits = mi > 0.0
        else:
            exits = self.lFill <= self.sealevel
        q = np.where(exits, self.riverSoluteSp[:, k], 0.0) / self.larea
        return self._carbPlumeSolve(q, L)

    def _carbTurbidity(self, L, source="transport"):
        """
        Turbid plume of the river mud (``turbidity`` response), spread along
        and off the coast by ``(I - L^2 lap) c = q``, land as identity rows.
        ``c`` is in m/yr of mud per unit area and decays as ``exp(-d / L)``
        away from its source ``q``:

        - ``source="transport"`` (default): the river mud the marine routing
          carried THROUGH each cell during the previous step (``_turbRate``,
          summed over the cascade passes by ``seaplex._distOcean``, per unit
          area and year): turbid water wherever mud is in transit, including
          over delta platforms the load crosses without settling. Zero
          before the first marine step.
        - ``source="mouths"``: the fine load the rivers deliver to the sea
          this step (``vSedFLocal`` at the marine cells the river routing
          reaches), so the plume spreads from the river mouths only.

        Collective. Scratch: ``tmp``, ``tmpL``, ``tmp1``.
        """
        sea = self.lFill <= self.sealevel
        if source == "mouths":
            q = np.where(sea, np.maximum(self.vSedFLocal.getArray(), 0.0), 0.0) / self.larea
        else:
            r = getattr(self, "_turbRate", None)
            q = np.zeros(self.lpoints) if r is None else np.where(sea, r, 0.0)
        return self._carbPlumeSolve(q, L)

    def _carbPlumeSolve(self, q, L):
        """Solve ``(I - L^2 lap) c = q`` with land as identity rows (the
        operator rebuilt each call, the coast moves; cached solver
        ``_ksp_plume``) and return the local ``max(c, 0)``. Shared by the
        ``solute`` and ``turbidity`` plumes. Collective. Scratch: ``tmp``,
        ``tmpL``, ``tmp1``."""
        land = self.hLocal.getArray() >= self.sealevel
        kd = np.full(self.lpoints, L * L, dtype=np.float64)
        if self._plumeMat is not None:
            self._plumeMat.destroy()
        self._plumeMat = self._buildDiffMat(kd, fixed=land)
        if self._ksp_plume is None:
            ksp = petsc4py.PETSc.KSP().create(petsc4py.PETSc.COMM_WORLD)
            ksp.setType("fgmres")
            ksp.getPC().setType("bjacobi")
            ksp.setTolerances(rtol=1.0e-10, atol=1.0e-30, max_it=2000)
            ksp.setOptionsPrefix("plume_")
            petsc4py.PETSc.Options()["plume_sub_pc_factor_shift_type"] = "nonzero"
            ksp.setFromOptions()
            self._ksp_plume = ksp
        self._ksp_plume.setOperators(self._plumeMat, self._plumeMat)
        self.tmpL.setArray(q)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self.tmp1.set(0.0)
        self._ksp_plume.solve(self.tmp, self.tmp1)
        self.dm.globalToLocal(self.tmp1, self.tmpL)
        return np.maximum(self.tmpL.getArray().copy(), 0.0)

    def carbSupplyLimit(self):
        """
        Supply-limited burial (``carbonate_supply``, phase C6). Before the
        production and the pelagic rain of the step: their potential burial
        (benthic production, plus the settled pelagic carbonate when the rain
        shares the supply) against what the marine reservoir of the species
        can make (``marineSolute * yield``). When the potential exceeds it,
        both are scaled down by the same factor (``_carbSupplyScale``); the
        burial actually made is removed from the reservoir as it is laid
        (``_carbDrawdown``). The benthic integration is kept and reused by
        ``carbProduction``. Collective (reductions).
        """
        k = self._carbSupplyK
        own = self.inIDs == 1
        demand = 0.0
        if getattr(self, "carbOn", False):
            prod = self._carbIntegrate()
            self._carbProdCache = prod
            tot = prod[0].copy()
            for p in prod[1:]:
                tot += p
            demand += float(np.sum((tot * self.larea)[own]))
        if getattr(self, "pelOn", False) and self.carbSupply["pelagic"]:
            e = self.pelEvents[self._pelEvent()]
            ev = self._pelMaps[self._pelEvent()]
            depth = self.sealevel - self.hLocal.getArray()
            marine = (self.lFill <= self.sealevel) & (depth > 0.0)
            rain = np.where(marine, ev["carb"] * self.dt, 0.0) * np.ones(self.lpoints)
            if e.get("temp") is not None:
                rain = rain * np.interp(self._sst(), e["temp"][0], e["temp"][1])
            carb = np.minimum(rain * self._pelPreservedFraction(depth, e["lys"], e["ccd"]),
                              np.maximum(depth, 0.0))
            demand += float(np.sum((carb * self.larea)[own]))
        demand = MPI.COMM_WORLD.allreduce(demand, op=MPI.SUM)
        avail = max(float(self.marineSolute[k]), 0.0) * self.carbSupply["yield"]
        s = 1.0 if demand <= avail or demand <= 0.0 else avail / demand
        self._carbSupplyScale = s
        self.carbSupplyBudget = np.array([avail, demand, s, float(self.marineSolute[k])])
        if MPIrank == 0 and self.verbose:
            print("[carbonate] supply: available %0.4g km3, demand %0.4g km3, scale %0.3f"
                  % (avail * 1e-9, demand * 1e-9, s), flush=True)

    def _carbDrawdown(self, vol):
        """Remove the carbonate volume ``vol`` (m3, already reduced over the
        ranks) buried this step from the marine reservoir of the supply
        species (no-op without ``carbonate_supply``)."""
        if getattr(self, "carbSupply", None) is None or vol <= 0.0:
            return
        k = self._carbSupplyK
        self.marineSolute[k] = max(float(self.marineSolute[k]) - vol / self.carbSupply["yield"],
                                   0.0)

    def _initSST(self):
        """Load the sea-surface temperature event maps (local slices)."""
        self._sstMaps = []
        for e in self.sstEvents:
            v = e["val"]
            if isinstance(v, tuple):
                d = np.load(v[0] + ".npz")
                arr = np.asarray(d[v[1]], dtype=np.float64)
                if arr.shape[0] != self.mpoints:
                    raise ValueError("sst: map %s has %d values, the mesh %d."
                                     % (list(v), arr.shape[0], self.mpoints))
                self._sstMaps.append(arr[self.locIDs])
            else:
                self._sstMaps.append(float(v))

    def _sst(self):
        """
        Sea-surface temperature (deg C) at ``tNow`` on the local nodes: the
        event in force (the last that has started; the first before it
        starts), plus the curve interpolated in time as an anomaly; or the
        curve alone (uniform). Rank-local.
        """
        out = np.zeros(self.lpoints, dtype=np.float64)
        if self.sstEvents:
            k = 0
            for i, e in enumerate(self.sstEvents):
                if e["start"] <= self.tNow + 1.0e-6:
                    k = i
            out = out + self._sstMaps[k]
        if self.sstCurve is not None:
            out = out + float(np.interp(self.tNow, self.sstCurve[0], self.sstCurve[1]))
        return out

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
        self._pelHand = None
        self._pelHanded = False
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

    def _clayDecay(self, L):
        """
        Fraction ``c`` of the coastal clay rain reaching each node, from
        ``(I - L^2 lap) c = 0`` on the sea with ``c = 1`` on land (``h`` at or
        above sea level, identity rows of the operator; open edges of a flat
        mesh are identity rows with ``c = 0``, the open ocean). In 1-D it is
        ``exp(-x/L)`` at distance ``x`` from the coast (exactly ``lambda^i``
        on a lattice, ``cosh(dx/L') = 1 + dx^2 / (2 L^2)``), and it goes around
        headlands and islands. All sea, no land: ``c = 0``. The operator is
        rebuilt every call (the coast moves); the solver is cached
        (``_ksp_clay``, prefix ``claydecay_``). Collective. Scratch: ``tmp``,
        ``tmpL``, ``tmp1``.

        :arg L: decay length (m)
        :return: local ``c`` in [0, 1]
        """
        land = self.hLocal.getArray() >= self.sealevel
        kd = np.full(self.lpoints, L * L, dtype=np.float64)
        if self._clayMat is not None:
            self._clayMat.destroy()
        self._clayMat = self._buildDiffMat(kd, fixed=land)
        if self._ksp_clay is None:
            # A decay length of tens of cells makes the operator stiff:
            # block-Jacobi needed hundreds of iterations (0.4 s a step, 15 % of
            # the step on a 47k-node margin). Algebraic multigrid (hypre, as
            # the groundwater head) where PETSc has it, at 1e-8 (the benchmark
            # checks the exact discrete solution to that, 2026-11-25).
            ksp = petsc4py.PETSc.KSP().create(petsc4py.PETSc.COMM_WORLD)
            ksp.setType("fgmres")
            if petsc4py.PETSc.Sys.hasExternalPackage("hypre"):
                ksp.getPC().setType("hypre")
                ksp.getPC().setHYPREType("boomeramg")
            else:
                ksp.getPC().setType("bjacobi")
                petsc4py.PETSc.Options()["claydecay_sub_pc_factor_shift_type"] = "nonzero"
            ksp.setTolerances(rtol=1.0e-8, atol=1.0e-14, max_it=2000)
            ksp.setOptionsPrefix("claydecay_")
            ksp.setFromOptions()
            self._ksp_clay = ksp
        self._ksp_clay.setOperators(self._clayMat, self._clayMat)
        self.tmpL.setArray(np.where(land, 1.0, 0.0))
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self.tmp1.set(0.0)
        self._ksp_clay.solve(self.tmp, self.tmp1)
        if self._ksp_clay.getConvergedReason() < 0 and MPIrank == 0:
            print("[pelagic] clay decay solve did not converge (%d)"
                  % self._ksp_clay.getConvergedReason(), flush=True)
        self.dm.globalToLocal(self.tmp1, self.tmpL)
        return np.clip(self.tmpL.getArray().copy(), 0.0, 1.0)

    def _pelPreservedFraction(self, depth, lys, ccd):
        """Fraction of the pelagic carbonate rain preserved at water depth
        ``depth``: 1 above the lysocline, 0 below the CCD, linear between."""
        if ccd <= lys:
            return np.where(depth <= lys, 1.0, 0.0)
        return np.clip((ccd - depth) / (ccd - lys), 0.0, 1.0)

    def _pelResuspend(self):
        """True when the event in force hands the rain to the marine
        diffusion (``resuspend``) and the marine step runs."""
        e = self.pelEvents[self._pelEvent()]
        return bool(e["resuspend"]) and getattr(self, "seaDepo", False) \
            and not getattr(self, "nodep", False)

    def pelagicRain(self, resuspend=False):
        """
        Pelagic carbonate and hemipelagic clay of the step on the marine nodes
        (``lFill`` at or below sea level).

        The carbonate rain (event rate, uniform or map) is preserved by the
        fraction set by the water depth (``_pelPreservedFraction``), the rest
        dissolves (budget only); the clay (optionally decaying away from the
        coast, ``_clayDecay``) settles entirely. The rain never builds above
        sea level (the excess is the ``bypassed`` budget item).

        Settled (default): called after the marine step; laid into the
        stratigraphy (``carb_mud`` and ``fine``, ``depoFineFrac = 1``) and the
        elevation. A node's deposit thinner than the deposit floor
        (``DEPOSIT_FLOOR``, round-off) is carried over to the next steps
        (``_pelCarry``).

        Resuspended (``resuspend`` event, ``resuspend=True``): called before
        the marine step; the volumes are kept in ``_pelHand`` and
        ``seaChange`` adds them to its diffusion input, each part in its
        lithology column. The call after the marine step is then skipped.

        Scratch Vecs: ``self.tmp`` / ``self.tmpL``. Collective (halo sync,
        budget, the clay decay solve).
        """
        if not resuspend and getattr(self, "_pelHanded", False):
            self._pelHanded = False
            return
        ev = self._pelMaps[self._pelEvent()]
        e = self.pelEvents[self._pelEvent()]
        hl = self.hLocal.getArray().copy()
        depth = self.sealevel - hl
        marine = (self.lFill <= self.sealevel) & (depth > 0.0)
        rain = np.where(marine, ev["carb"] * self.dt, 0.0) * np.ones(self.lpoints)
        if e.get("temp") is not None:
            # pelagic production scaled by the sea-surface temperature
            rain = rain * np.interp(self._sst(), e["temp"][0], e["temp"][1])
        if getattr(self, "carbSupply", None) is not None and self.carbSupply["pelagic"]:
            rain = rain * self._carbSupplyScale
        keep = self._pelPreservedFraction(depth, e["lys"], e["ccd"])
        carb0 = rain * keep
        dissolved = rain - carb0
        clay0 = np.where(marine, ev["clay"] * self.dt, 0.0) * np.ones(self.lpoints)
        if e["decay"] is not None:
            # clay rain decreasing away from the coast
            clay0 = clay0 * self._clayDecay(e["decay"])
        # never above sea level: what does not fit is bypassed (budget)
        tot = carb0 + clay0
        scale = np.divide(np.minimum(tot, np.maximum(depth, 0.0)), tot,
                          out=np.zeros_like(tot), where=tot > 0.0)
        carb = carb0 * scale
        clay = clay0 * scale
        bypass = tot - carb - clay

        if resuspend:
            # handed to this step's marine diffusion (seaChange)
            self._pelHand = (carb * self.larea, clay * self.larea)
            self._pelHanded = True
            self.pelDep = (carb + clay) / self.dt
        else:
            self._pelSettle(carb, clay)

        own = self.inIDs == 1
        a = self.larea
        b = MPI.COMM_WORLD.allreduce(
            np.array([(rain * a)[own].sum(), (carb * a)[own].sum(),
                      (dissolved * a)[own].sum(), (clay * a)[own].sum(),
                      (bypass * a)[own].sum()]), op=MPI.SUM)
        # (carbonate rain, carbonate settled or resuspended, dissolved, clay,
        #  bypassed above sea level), m3 over the step
        self.pelBudget = b
        if getattr(self, "carbSupply", None) is not None and self.carbSupply["pelagic"]:
            self._carbDrawdown(b[1])
        self._pelRain += b[0]
        self._pelPreserved += b[1]
        self._pelDissolved += b[2]
        self._pelClay += b[3]
        self._pelBypass = getattr(self, "_pelBypass", 0.0) + b[4]
        if MPIrank == 0 and self.verbose:
            print(
                "[pelagic] carbonate rain %0.4g km3 (preserved %0.4g, dissolved "
                "%0.4g), clay %0.4g km3%s" % (b[0] * 1e-9, b[1] * 1e-9, b[2] * 1e-9,
                                              b[3] * 1e-9,
                                              " (resuspended)" if resuspend else ""),
                flush=True,
            )
        return

    def _pelSettle(self, carb, clay):
        """Lay the settled rain (thicknesses ``carb``, ``clay``) into the
        elevation and the stratigraphy, carrying over sub-floor deposits."""
        self._pelCarry[:, 0] += carb
        self._pelCarry[:, 1] += clay
        tot = self._pelCarry.sum(axis=1)
        lay = tot >= DEPOSIT_FLOOR
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
            self._deposeMarineProv()
            self.depoFineFrac = ff
            if sh is not None:
                self.litShare = sh

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
        map times the ``slope``, ``sedrate``, ``exposure``, ``temperature``,
        ``wave``, ``solute`` and ``turbidity`` responses
        (the open-water index is solved once per step, only when used)."""
        out = []
        slope = None
        expo = None
        sst = None
        wave = None
        plume = {}
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
                elif ctl == "temperature":
                    if sst is None:
                        sst = self._sst()
                    m = m * np.interp(sst, x, fv)
                elif ctl == "wave":
                    if wave is None:
                        wave = self.waveEnergy()
                    m = m * np.interp(wave, x, fv)
                elif ctl == "solute":
                    key = f["solute"]
                    if key not in plume:
                        plume[key] = self._carbPlume(*key)
                    m = m * np.interp(plume[key], x, fv)
                elif ctl == "turbidity":
                    key = ("turbidity",) + tuple(f["turbidity"])
                    if key not in plume:
                        plume[key] = self._carbTurbidity(*f["turbidity"])
                    m = m * np.interp(plume[key], x, fv)
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
        if self._carbProdCache is not None:
            prod = self._carbProdCache
            self._carbProdCache = None
        else:
            prod = self._carbIntegrate()
        s = self._carbSupplyScale
        if s < 1.0:
            prod = [p * s for p in prod]
        self._carbFinish(prod, t0)
        return

    def _carbIntegrate(self):
        """Integrate the production of every factory over the step (Heun
        sub-steps, sea-level cap); returns the thickness produced per factory.
        Collective only through the response fields it may solve."""
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
        return prod

    def _carbFinish(self, prod, t0):
        """Split the production ``prod`` into framework and loose sediment,
        build the framework, keep the loose part for the marine step, update
        the budgets (and draw the supply reservoir down)."""
        total = np.zeros(self.lpoints)
        frame = np.zeros(self.lpoints)
        sand = np.zeros(self.lpoints)
        for f, p in zip(self.carbFact, prod):
            total += p
            frame += f["phi"] * p
            sand += (1.0 - f["phi"]) * f["sand"] * p
        # Framework thinner than the deposit floor of deposeStrat (round-off)
        # would raise the elevation without being recorded: it joins the
        # loose sand instead (conservative).
        thin = frame < DEPOSIT_FLOOR
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
            self._deposeMarineProv()
            self.depoFineFrac, self.litShare = ff, sh

        own = self.inIDs == 1
        b = MPI.COMM_WORLD.allreduce(
            np.array([(total * self.larea)[own].sum(), (frame * self.larea)[own].sum(),
                      (loose * self.larea)[own].sum()]), op=MPI.SUM)
        self.carbBudget = b
        self._carbDrawdown(b[0])
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

    def _initKarst(self):
        """Allocate the karst state (docs/DESIGN_KARST.md): the cumulative
        dissolution ``karstDiss`` (m), the reference area of the focusing
        (the mean cell area by default) and the geochemistry species fed.
        Collective (one reduction)."""
        k = self.karst
        self.karstDiss = np.zeros(self.lpoints, dtype=np.float64)
        self._karstSolid = np.zeros(len(k["lit"]), dtype=np.float64)
        self._karstBulk = 0.0
        if k["area"] is None:
            own = self.inIDs == 1
            a = MPI.COMM_WORLD.allreduce(np.array([self.larea[own].sum(), own.sum()],
                                                  dtype=np.float64), op=MPI.SUM)
            k["area"] = float(a[0] / max(a[1], 1.0))
        self._karstSp = None
        if k["species"] is not None:
            names = list(getattr(self, "gwGeoName", []) or [])
            if not getattr(self, "gwGeochemOn", False) or \
                    not getattr(self, "gwRiverLoad", False) or not names:
                raise ValueError(
                    "karst: `species` needs the groundwater geochemistry with "
                    "`river_load: True` (the river dissolved load).")
            if k["species"] not in names:
                raise ValueError("karst: geochemistry species '%s' not found "
                                 "(species: %s)." % (k["species"], ", ".join(names)))
            self._karstSp = names.index(k["species"])

    def _karstRate(self):
        """Karst dissolution rate (m/yr of rock, bulk) per node,
        ``rate * R * max(Q / (R A0), 1)^focus * f_T`` on subaerial land
        (not the sea, lakes or ice), ``R`` the runoff (rain minus
        evaporation) and ``Q`` the discharge of this step's
        ``flowAccumulation``. Rank-local."""
        k = self.karst
        h = self.hLocal.getArray()
        rain = np.asarray(self.rainVal, dtype=np.float64) * np.ones(self.lpoints)
        evap = getattr(self, "evapVal", None)
        R = np.maximum(rain - (0.0 if evap is None else evap), 0.0)
        land = h > self.sealevel
        land &= ~self._underLake(h)
        if getattr(self, "iceOn", False) and getattr(self, "iceHL", None) is not None:
            land &= self.iceHL.getArray() <= ICE_COVER_MIN
        D = k["rate"] * R
        if k["focus"] > 0.0:
            Q = np.maximum(self.FAL.getArray(), 0.0)
            ratio = np.divide(Q, R * k["area"], out=np.ones_like(Q), where=R > 0.0)
            D = D * np.maximum(ratio, 1.0) ** k["focus"]
        if k["temp"] is not None:
            D = D * np.interp(self._sst(), k["temp"][0], k["temp"][1])
        return np.where(land, D, 0.0)

    def karstDissolution(self):
        """
        Karst (docs/DESIGN_KARST.md): dissolve the carbonate lithologies of
        the exposed stratigraphy at the rate of ``_karstRate`` (capped at
        ``max_dissolution`` per step), lithology by lithology from the top
        layer down (``dissolveStrat``: the insoluble residue stays, a residue
        ``cover`` stops it). Lowers the elevation, the cumulative erosion and
        the rate field; the dissolved solid volume (m3/yr) is kept in
        ``_karstSolute`` (times ``solute_per_m3``) for the river routing of
        the geochemistry species. Collective (halo sync, budget).
        Scratch: ``tmp``, ``tmpL``.
        """
        k = self.karst
        d = np.clip(self._karstRate() * self.dt, 0.0, k["max"])
        # halo: the owners' demand on the ghosts
        self.tmpL.setArray(d)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self.dm.globalToLocal(self.tmp, self.tmpL)
        d = self.tmpL.getArray().copy()
        removed, solid = self.dissolveStrat(d, k["lit"], k["cover"])

        self.tmpL.setArray(-removed)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self.cumED.axpy(1.0, self.tmp)
        self.hGlobal.axpy(1.0, self.tmp)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal)
        self.dm.globalToLocal(self.hGlobal, self.hLocal)
        self.tmpL.setArray(-removed / self.dt)
        self.EbLocal.axpy(1.0, self.tmpL)
        self.tmp.set(0.0)
        self.tmpL.set(0.0)
        self.karstDiss += removed

        svol = solid * self.larea[:, None]
        self._karstSolute = svol.sum(axis=1) * k["solute"] / self.dt
        own = self.inIDs == 1
        b = MPI.COMM_WORLD.allreduce(
            np.concatenate([[float((removed * self.larea)[own].sum())],
                            svol[own].sum(axis=0)]), op=MPI.SUM)
        self.karstBudget = b
        self._karstBulk += b[0]
        self._karstSolid += b[1:]
        if MPIrank == 0 and self.verbose:
            print("[karst] dissolved %0.4g km3 (solid %0.4g km3)"
                  % (b[0] * 1e-9, b[1:].sum() * 1e-9), flush=True)
        return
