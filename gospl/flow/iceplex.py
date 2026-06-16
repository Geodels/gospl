import os
import petsc4py
import numpy as np

from mpi4py import MPI
from time import process_time

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import ice_flux
    from gospl._fortran import ice_flux_limiter
    from gospl._fortran import ice_flux_rscaled
    from gospl._fortran import ice_velocity
    from gospl._fortran import mfdreceivers

# Ice density (kg/m^3), used for the SIA deformation coefficient and loading.
RHO_ICE = 910.0

# Reference thickness (m) for the accumulation substep limit: ice growing from
# zero in the accumulation zone is resolved in increments scaled to this depth so
# the SIA flux can drain it as it thickens (rather than one giant substep dumping
# the whole step's accumulation before any flow — which produces runaway ice).
ICE_HREF = 100.0

MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()


class IceMesh(object):
    r"""
    Shallow-Ice-Approximation (SIA) ice-sheet model.

    This class evolves the ice thickness :math:`H` as a non-linear diffusion of
    the ice surface :math:`s = z_\mathrm{bed} + H`,

    .. math::
        \frac{\partial H}{\partial t} = \dot{m} - \nabla \cdot \mathbf{q},

    where the ice flux :math:`\mathbf{q}` follows Glen's-law internal
    deformation plus basal sliding (the ``ice_flux`` Fortran kernel) and
    :math:`\dot{m}` is the equilibrium-line-altitude (ELA) surface mass balance.
    The thickness is integrated **explicitly** with a mass-conserving,
    positivity-preserving flux limiter (the ``ice_flux_limiter`` /
    ``ice_flux_rscaled`` kernels cap each cell's outflux to the ice it holds, so
    :math:`H\geq 0` for any substep and ice volume is conserved to machine
    precision). The goSPL timestep is split
    into accuracy-controlled substeps (``sia.cfl``); the measured stable substep
    at goSPL's km resolution is :math:`10^3`–:math:`10^4` yr, so only a handful
    of substeps are needed. The free boundary :math:`H\geq 0` (an obstacle
    problem at the ice margin) is handled natively by the limiter rather than by
    a post-hoc clamp, which would inject mass.

    From the converged thickness the model derives the basal sliding speed
    (driving velocity-based glacial abrasion :math:`E_g = K_g |u_b|^l`), the
    ablation meltwater re-injected into the river flow accumulation, glacial
    till (abraded rock transported and deposited as moraine in the ablation
    zone, conserving rock volume), and the ice load feeding the flexural
    isostasy.

    Useful links:
    - `Braun et al. (1999) <https://doi.org/10.3189/172756499781821797>`_
    - `Herman & Braun (2008) <https://doi.org/10.1029/2007JF000807>`_
    - `Egholm (2012) <https://doi.org/10.1016/j.geomorph.2011.12.019>`_
    - `Deal & Prasicek (2020) <https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2020GL089263>`_
    - `Hergarten (2021) <http://dx.doi.org/10.5194/esurf-2021-1>`_
    - `Liebl et al. (2023) <https://gmd.copernicus.org/articles/16/1315/2023/>`_
    """

    def __init__(self, *args, **kwargs):
        """
        Initialisation of the `IceMesh` class.

        This method initializes the ice-related fields (ice thickness, basal
        velocity, meltwater and flexural response) based on the configuration
        flags `iceOn` and `flexOn`. Memory is allocated for these fields.
        """

        if self.iceOn:
            self.iceHL = self.hLocal.duplicate()
            # Local meltwater field captured at the end of iceAccumulation
            # and re-injected into the river FA source term in flowplex so
            # glacial discharge is not lost downstream.
            self.iceMeltL = self.hLocal.duplicate()
            if self.flexOn:
                self.iceFlex = self.hLocal.duplicate()
            # Basal sliding speed (m/yr) from the SIA solve; the driver for the
            # velocity-based glacial abrasion law. Registered in destroy_DMPlex.
            self.iceUbL = self.hLocal.duplicate()
            # Glacial abrasion rate E_g = Kg|u_b|^l (m/yr); a diagnostic output
            # field, populated every step from the basal velocity (zero where
            # abrasion is off, i.e. Kg = 0).
            self.iceAbrL = self.hLocal.duplicate()
            self.iceHL.set(0.0)
            self.iceMeltL.set(0.0)
            self.iceUbL.set(0.0)
            self.iceAbrL.set(0.0)
            # Seed a pre-existing ice thickness on a fresh start; the SIA solve
            # evolves it from there. On a restart (rStep > 0) readData restores
            # the evolved iceH instead, so the seed is skipped. The flexure
            # reference is taken AFTER this seed (first iceAccumulation), so a
            # pre-existing ice load is not applied as a transient shock.
            if self._iceInitSpec is not None and self.rStep == 0:
                sc, spec = self._iceInitSpec
                if spec is not None:
                    h0 = self._loadIceMap(spec, "ice hinit")[self.locIDs]
                else:
                    h0 = np.full(self.lpoints, sc, dtype=np.float64)
                self.iceHL.setArray(np.maximum(h0, 0.0))
            # Glacial-till mass-balance diagnostics (m^3, owned-node running
            # totals reduced by the till-conservation test).
            self._tillEroded = 0.0
            self._tillDeposited = 0.0

        return

    def _iceSIAParams(self, elaH, iceH):
        """
        Shared SIA setup: Glen exponent `n`, deformation coeff `ad`, sliding
        coeff `as`, bed elevation `zbed`, and the ELA surface mass balance
        `mdot` (m ice/yr, no `meltfac`). Used by both SIA solvers.
        """
        n = self.sia_glen
        # Deformation coefficient ad = 2 A (rho_ice g)^n / (n+2) (SIA flux form).
        ad = 2.0 * self.sia_Aglen * (RHO_ICE * self.gravity) ** n / (n + 2.0)
        as_ = self.sia_slide
        zbed = self.hLocal.getArray().copy()              # bed elevation (fixed this step)
        # ELA mass-balance ramp (1 above the ice cap, 0 at the ELA, negative
        # below). elaH/iceH may be scalars (uniform) or per-vertex arrays
        # (spatial maps). Where iceH <= elaH (degenerate / no-ice-band cell) the
        # ramp is zeroed so the division never blows up; this is exact for the
        # uniform case, where the caller's guard already ensures iceH > elaH.
        denom = iceH - elaH
        safe = np.where(denom > 0.0, denom, 1.0)
        ramp = np.where(denom > 0.0, (zbed - elaH) / safe, 0.0)
        ramp = np.minimum(ramp, 1.0)
        mdot = self.rainVal * ramp                        # surface mass balance (m ice/yr)
        return n, ad, as_, zbed, mdot

    def _iceSIAFinalize(self, H, zbed, mdot, iceT, as_, n):
        """
        Common post-solve steps for the SIA solve: terminus clamp, store the
        ice thickness and the basal sliding speed (the abrasion driver), and
        capture ablation meltwater (m^3/yr) for the river coupling.

        The terminus floor is ``max(hterm, sea level)``: ice is removed below
        the prescribed terminus and, regardless, below the (possibly
        time-varying) sea surface — no land ice persists offshore (there is no
        marine-ice / calving model). When ``hterm`` is unprescribed (the
        ``TERMINUS_UNSET`` sentinel) the floor is simply the sea-level position.
        """
        H = np.maximum(H, 0.0)
        terminus = np.maximum(iceT, self.sealevel)
        H[zbed < terminus] = 0.0
        self.iceHL.setArray(H)
        # Basal sliding speed from the converged thickness (steepest-descent
        # SIA velocity); zero where there is no ice.
        ub = ice_velocity(self.lpoints, H, zbed, as_, n)
        ub[H <= 1.0e-2] = 0.0
        self.iceUbL.setArray(ub)
        # Glacial abrasion rate diagnostic E_g = Kg|u_b|^l (m/yr), subaerial
        # only; matches the incision applied by _glacialAbrasion / glacialTill.
        abr = self.ice_Kg * np.power(np.maximum(ub, 0.0), self.ice_abr_l)
        abr[self.seaID] = 0.0
        self.iceAbrL.setArray(abr)
        melt_local = np.maximum(-mdot, 0.0) * self.larea * (H > 1.0e-2)
        self.iceMeltL.setArray(melt_local)
        return

    def _iceFlowSIA(self, elaH, iceH, iceT):
        r"""
        SIA ice-thickness solve: explicit, mass-conserving, positivity-preserving
        integration of ``∂H/∂t = ṁ − ∇·q``.

        Each forward substep is ``H ← H + Δt_sub (ṁ − ∇·q)``. ``ice_flux_limiter``
        builds the per-cell factor capping its outflux to the ice it holds; it is
        halo-synced so a partition-boundary cell scales identically on both ranks,
        then ``ice_flux_rscaled`` assembles the R-scaled divergence. The limiter
        (i) conserves ice volume to machine precision — the scaled edge flux is
        applied equal-and-opposite to the two cells of a shared face — and
        (ii) keeps ``H ≥ 0`` for **any** substep, so the free
        boundary at the ice margin (an obstacle problem) needs no mass-injecting
        post-clamp. A clamp is applied only after surface ablation, which can
        physically empty a cell (this removes over-melt, it never creates ice).

        Because positivity is unconditional, the substep size ``Δt_sub`` is an
        *accuracy* choice (``sia.cfl``), capped to ``sia.max_substeps`` substeps.
        It resolves two rates: the per-cell time to lose ice to flux+ablation
        (capped — the limiter/clamp keep those safe), and the time to *gain* ice
        by accumulation relative to a reference depth ``ICE_HREF`` (NOT capped —
        ice growing from zero must thicken gradually so the SIA flux can drain it
        as it builds; otherwise one giant substep dumps a whole step's
        accumulation before any flow, producing runaway km-thick ice). For modest
        accumulation only a handful of substeps are needed; high accumulation
        (m/yr) over a large ``dt`` is genuinely stiff and needs many. Validated
        against the analytical SIA dome (volume conservation under zero balance).
        """
        n, ad, as_, zbed, mdot = self._iceSIAParams(elaH, iceH)

        H = self.iceHL.getArray().copy()
        cfl = self.sia_cfl
        maxsub = self.sia_max_substeps
        Hmin = 1.0  # m: thin-cell floor for the substep estimate (the limiter
        #             guarantees positivity for thinner cells regardless)
        tiny = 1.0e-30

        remaining = self.dt
        nsub = 0
        while remaining > 1.0e-6 * self.dt and nsub < maxsub:
            # Halo-sync the thickness so the neighbour stencil in ice_flux* sees
            # up-to-date ghost values, then read back the synced local array.
            # NB: use the scratch global vec `self.tmp`, never `self.hGlobal`
            # (that holds the bed elevation read below as `zbed`).
            self.tmpL.setArray(H)
            self.dm.localToGlobal(self.tmpL, self.tmp)
            self.dm.globalToLocal(self.tmp, self.tmpL)
            # Copy: `tmpL` is reused below for the rfac halo-sync, which would
            # otherwise alias-clobber this thickness array.
            H = self.tmpL.getArray().copy()

            # Accuracy CFL: pick a substep so a cell changes by at most `cfl` of
            # its ice per substep. Two contributions:
            raw = ice_flux(self.lpoints, H, zbed, ad, as_, n)
            thick = H > Hmin
            # (1) flux outflux + ablation loss. A cell can lose at most its own
            #     ice in a step — the flux limiter caps the outflux and the melt
            #     clamp caps ablation — so a cell shedding >=100% of its ice is
            #     already handled conservatively and must NOT drive the substep to
            #     zero (a real stall: thin cells on steep slopes / strong ablation
            #     gave rate->inf). Cap this rate against the FULL `dt` (not
            #     `remaining`, which would shrink and crawl geometrically).
            loss = np.where((raw > 0.0) & thick, raw, 0.0)
            melt = np.where((mdot < 0.0) & thick, -mdot, 0.0)
            loss_rate = np.minimum((loss + melt) / np.maximum(H, Hmin), 1.0 / self.dt)
            # (2) accumulation gain. This is NOT capped: a thin/zero cell in the
            #     accumulation zone must add ice gradually so the growing ice can
            #     flow and drain as it thickens — otherwise one giant substep dumps
            #     the whole step's accumulation before any drainage (the cause of
            #     runaway km-thick ice). The reference thickness `ICE_HREF` sets how
            #     finely growth from zero is resolved.
            acc = np.where(mdot > 0.0, mdot, 0.0)
            acc_rate = acc / np.maximum(H, ICE_HREF)
            rate = loss_rate + acc_rate
            gmax = MPI.COMM_WORLD.allreduce(
                float(rate.max()) if rate.size else 0.0, op=MPI.MAX
            )
            sub = remaining if gmax <= tiny else min(remaining, cfl / gmax)

            # Flux limiter factor (caps each cell's outflux to its ice over the
            # substep). Only owned cells have the full neighbour stencil, so
            # halo-sync rfac to ghosts before the divergence — otherwise a
            # limited (rfac<1) cell on a partition boundary would break
            # conservation (the two ranks would scale the shared face differently).
            rfac = ice_flux_limiter(self.lpoints, H, zbed, ad, as_, n, sub)
            self.tmpL.setArray(rfac)
            self.dm.localToGlobal(self.tmpL, self.tmp)
            self.dm.globalToLocal(self.tmp, self.tmpL)
            rfac = self.tmpL.getArray()

            # Conservative, positivity-preserving flux divergence for this substep.
            val = ice_flux_rscaled(self.lpoints, H, zbed, ad, as_, n, rfac)
            H = H + sub * (mdot - val)
            # Only ablation can drive a cell negative (flux cannot); clamp it —
            # this stops over-melt, it does not create ice.
            H = np.maximum(H, 0.0)

            remaining -= sub
            nsub += 1

        # Substeps taken this step (diagnostic; the limiter bounds this at ~1/cfl
        # so it must not approach max_substeps — a regression guard for the
        # thin-cell / strong-ablation substep stall).
        self._sia_nsub = nsub
        if nsub >= maxsub and remaining > 1.0e-6 * self.dt and MPIrank == 0:
            print(
                "SIA explicit ice solve hit the substep cap (%d) with %.3g yr of "
                "the %.3g yr step left; raise ice.sia.max_substeps or lower "
                "ice.sia.cfl." % (maxsub, remaining, self.dt),
                flush=True,
            )

        self._iceSIAFinalize(H, zbed, mdot, iceT, as_, n)

        return

    def iceAccumulation(self):
        """
        Main Ice Accumulation Calculation.

        This method evolves the ice thickness with the SIA solve
        (:meth:`_iceFlowSIA`) using the ice-cap altitude, the equilibrium-line
        altitude (ELA) and the glacier terminus. Each of these may be a uniform
        scalar (time-varying) or a per-vertex map — the latter for global models
        where the ELA varies strongly with latitude. It also snapshots the ice
        load for the flexural isostasy when enabled.
        """

        ti = process_time()

        # Snapshot the current ice load for flexure (applyFlexure loads the
        # increment iceHL - iceFlex), taken BEFORE the thickness is updated.
        if self.flexOn:
            self.iceHL.copy(result=self.iceFlex)

        # Resolve the glacier-geometry altitudes: per-vertex map where supplied
        # (spatial ELA for global tropical-vs-polar runs), otherwise the uniform
        # time-varying scalar. `spatial` flags that any field is a map.
        iceH = self.iceMesh[self.locIDs] if self.iceMesh is not None else self.iceH(self.tNow)
        elaH = self.elaMesh[self.locIDs] if self.elaMesh is not None else self.elaH(self.tNow)
        iceT = self.termMesh[self.locIDs] if self.termMesh is not None else self.iceT(self.tNow)
        spatial = (
            self.iceMesh is not None
            or self.elaMesh is not None
            or self.termMesh is not None
        )

        # Uniform case keeps the cheap global short-circuit (byte-identical):
        # no ice when the whole surface is below the ELA, or for a degenerate
        # config where the ice-cap altitude is not above the ELA. With spatial
        # maps these are per-node and handled robustly in _iceSIAParams.
        if not spatial:
            max_elev = self.hGlobal.max()[1]
            if max_elev < elaH or iceH <= elaH:
                self.iceHL.set(0.)
                self.iceUbL.set(0.)
                self.iceAbrL.set(0.)
                self.iceMeltL.set(0.)
                if self.flexOn:
                    self.iceFlex.set(0.)
                return

        # Implicit SIA ice-thickness solve.
        self._iceFlowSIA(elaH, iceH, iceT)

        # At the first step, seed the flexure reference with the ice just
        # solved so the initial load is not applied as a transient.
        if self.flexOn and self.tNow == self.tStart:
            self.iceHL.copy(result=self.iceFlex)

        if MPIrank == 0 and self.verbose:
            print(
                "Glaciers Accumulation - SIA (%0.02f seconds)"
                % (process_time() - ti),
                flush=True,
            )

        return

    def glacialTill(self):
        r"""
        Glacial till — production, transport and deposition.

        Glacial abrasion (:math:`E_g = K_g |u_b|^l`) erodes rock under sliding
        ice; the abraded material becomes **till** carried by the ice and
        deposited where the ice **melts out** — the ablation zone / terminal
        moraine. The bed is therefore **lowered** under fast ice and **raised**
        in the ablation zone: a transport, not a source/sink.

        Two modes, selected by whether stratigraphy is on:

        - **Bulk bed** (no stratigraphy): the abraded volume ``Vtot`` is
          redistributed to the ablation cells weighted by the meltwater rate
          (``iceMeltL``) as a bed-to-bed transport, so ``Σ deposited == Vtot``
          and the net bed-volume change is zero by construction.
        - **Stratigraphic** (:meth:`_glacialTillStrata`): the till is removed
          from the layers it was abraded from and re-deposited as a fresh
          moraine layer in the ablation zone — split into the coarse/fine
          lithology fractions under dual lithology. Conservation is on the
          *solid* phase (per fraction), so the bed bulks up by the porosity
          contrast (uncompacted till vs compacted source rock).

        Both modes are guarded by a dedicated till-conservation test (the
        dual-lithology lesson: a volume the total budget can't see needs its
        own guard).

        Active only with ``till.on`` and ``Kg > 0``; a no-op otherwise.

        **Deposition distribution.** By default the till is spread across the
        whole ablation zone weighted by the meltwater rate. With
        ``till.route: True`` it is instead **routed down the ice-surface flow
        network** and melts out toward each catchment's terminus
        (:meth:`_routeTill`) — appropriate for high-resolution regional runs
        where individual glacier catchments and termini are resolved. Both
        produce a per-cell deposition weight (summing to one) consumed
        identically by the bulk and stratigraphic paths, so conservation is
        unchanged.
        """
        if not (self.iceOn and self.ice_till_on) or self.ice_Kg <= 0.0:
            return

        ti = process_time()
        ub = self.iceUbL.getArray()
        Eg = self.ice_Kg * np.power(np.maximum(ub, 0.0), self.ice_abr_l)  # m/yr
        Eg[self.seaID] = 0.0
        Vero = Eg * self.dt * self.larea                  # abraded volume (m^3) per cell
        owned = self.inIDs == 1
        Vtot = MPI.COMM_WORLD.allreduce(float(np.sum(Vero[owned])), op=MPI.SUM)

        # No abrasion -> nothing happens (skipping the erosion too keeps it
        # conservative: no orphaned till).
        if Vtot <= 0.0:
            return

        # Per-cell deposition weight (local array, Σ over owned nodes = 1):
        # melt-weighted spreading, or catchment-routed melt-out to the termini.
        if self.ice_till_route:
            dep_w = self._routeTill(Vero, Vtot, owned)
        else:
            melt = self.iceMeltL.getArray()               # ablation volume (m^3/yr)
            Wtot = MPI.COMM_WORLD.allreduce(
                float(np.sum(melt[owned])), op=MPI.SUM
            )
            # No ablation zone to receive the melt-spread till -> no-op.
            if Wtot <= 0.0:
                return
            dep_w = melt / Wtot

        # Routing can leave nothing depositable (e.g. no resolved outlet); guard.
        Wdep = MPI.COMM_WORLD.allreduce(float(np.sum(dep_w[owned])), op=MPI.SUM)
        if Wdep <= 0.0:
            return

        dz_ero = -Eg * self.dt                            # bed lowering (m, ≤0)

        if self.stratNb > 0:
            # Stratigraphy on: route the till through the stratigraphic pile so
            # the abraded material leaves the layers it came from and the
            # moraine is recorded as a fresh deposit (split coarse/fine under
            # dual lithology). See _glacialTillStrata.
            self._glacialTillStrata(dz_ero, dep_w, owned)
        else:
            # Bulk bed transport: erode at abrasion cells, deposit the till per
            # the deposition weight. Σ(dz·area) == 0 (rock moved, not created).
            dz_dep = (Vtot * dep_w) / self.larea
            dz = dz_ero + dz_dep
            self.tmpL.setArray(dz)
            self.dm.localToGlobal(self.tmpL, self.tmp)
            self.cumED.axpy(1.0, self.tmp)
            self.hGlobal.axpy(1.0, self.tmp)
            self.dm.globalToLocal(self.cumED, self.cumEDLocal)
            self.dm.globalToLocal(self.hGlobal, self.hLocal)
            self._tillEroded += Vtot
            self._tillDeposited += MPI.COMM_WORLD.allreduce(
                float(np.sum((dz_dep * self.larea)[owned])), op=MPI.SUM
            )

        if MPIrank == 0 and self.verbose:
            print(
                "Glacial Till (%0.02f seconds)" % (process_time() - ti),
                flush=True,
            )

        return

    def _glacialTillStrata(self, dz_ero, dep_w, owned):
        r"""
        Stratigraphic / dual-lithology coupling for glacial till.

        The abraded bed lowering ``dz_ero`` (m, ≤0) is removed from the
        stratigraphic pile with :meth:`erodeStrat`, which returns the
        uncompacted solid removed per fraction (``thCoarse`` / ``thFine``,
        m/yr). That solid is the till: its total uncompacted volume is
        redistributed to the ablation cells weighted by the meltwater rate and
        laid down as a fresh moraine layer with :meth:`deposeStrat`, carrying
        the **abraded fine fraction** (a single, ice-mixed composition).

        Solid mass is conserved per fraction by construction: the fine volume
        deposited equals the fine volume eroded (so the dual-lithology
        ``_fineEroded`` / ``_fineDeposited`` budget stays balanced), and the
        coarse volume likewise. The bed bulks up by the porosity contrast
        between the uncompacted till and the compacted source rock.

        ``erodeStrat`` / ``deposeStrat`` overwrite the transient routing arrays
        (``thCoarse`` / ``thFine`` / ``depoFineFrac``) that the *fluvial*
        ``sedChange`` consumes later this step, so they are saved and restored
        around the till calls.
        """
        # Preserve the fluvial routing state (set by the eroder's erodeStrat,
        # consumed by the later sedChange/_getSedFlux). glacialTill runs
        # between the two.
        saved_thCoarse = self.thCoarse.copy() if hasattr(self, "thCoarse") else None
        saved_thFine = (
            self.thFine.copy()
            if self.stratLith and hasattr(self, "thFine")
            else None
        )
        saved_depoFineFrac = self.depoFineFrac.copy()

        # --- Erosion: remove the abraded bulk from the stratigraphic pile. ---
        self.tmpL.setArray(dz_ero)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self.erodeStrat()                                 # sets thCoarse/thFine
        # Apply the bed lowering to the elevation / cumulative-erosion fields.
        self.cumED.axpy(1.0, self.tmp)
        self.hGlobal.axpy(1.0, self.tmp)

        # Uncompacted solid abraded (m^3): the till volume to redeposit.
        coarseV = MPI.COMM_WORLD.allreduce(
            float(np.sum((self.thCoarse * self.dt * self.larea)[owned])), op=MPI.SUM
        )
        if self.stratLith:
            fineV = MPI.COMM_WORLD.allreduce(
                float(np.sum((self.thFine * self.dt * self.larea)[owned])), op=MPI.SUM
            )
        else:
            fineV = 0.0
        Vsolid = coarseV + fineV
        if Vsolid <= 0.0:
            # Nothing actually came out of the pile (e.g. emptied to bedrock
            # floor): no till to deposit. Restore and return.
            self.thCoarse = saved_thCoarse
            if saved_thFine is not None:
                self.thFine = saved_thFine
            self.depoFineFrac = saved_depoFineFrac
            self.dm.globalToLocal(self.cumED, self.cumEDLocal)
            self.dm.globalToLocal(self.hGlobal, self.hLocal)
            self._tillEroded += 0.0
            return

        # Ice-mixed till composition (single fine fraction for the moraine).
        ffrac = fineV / Vsolid if self.stratLith else 0.0

        # --- Deposition: lay the till down per the deposition weight. ---
        dz_dep = (Vsolid * dep_w) / self.larea            # bulk till thickness (m)
        self.depoFineFrac = np.zeros(self.lpoints, dtype=np.float64)
        self.depoFineFrac[dz_dep > 0.0] = ffrac
        # Snapshot the current layer so the bed/cumED update uses the thickness
        # deposeStrat ACTUALLY laid down (it floors sub-1e-4 m deposits), keeping
        # elevation, the stratigraphic pile and the diagnostics consistent.
        h_before = self.stratH[:, self.stratStep].copy()
        self.tmpL.setArray(dz_dep)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self.deposeStrat()                                # adds to stratH/stratHf
        dep = self.stratH[:, self.stratStep] - h_before   # actual deposit (m)
        self.tmpL.setArray(dep)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self.cumED.axpy(1.0, self.tmp)
        self.hGlobal.axpy(1.0, self.tmp)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal)
        self.dm.globalToLocal(self.hGlobal, self.hLocal)

        # Diagnostics: solid abraded vs solid actually re-deposited (the two
        # match to the deposeStrat thickness floor; any sub-floor remainder is
        # dropped, as for every other deposit in the model).
        self._tillEroded += Vsolid
        self._tillDeposited += MPI.COMM_WORLD.allreduce(
            float(np.sum((dep * self.larea)[owned])), op=MPI.SUM
        )

        # Restore the fluvial routing state for sedChange.
        if saved_thCoarse is not None:
            self.thCoarse = saved_thCoarse
        if saved_thFine is not None:
            self.thFine = saved_thFine
        self.depoFineFrac = saved_depoFineFrac

        return

    def _routeTill(self, Vero, Vtot, owned):
        r"""
        Catchment-aware till routing on the ice-surface flow network
        (``till.route: True``).

        The abraded till is transported down the **ice surface**
        :math:`s = z_\mathrm{bed} + H` along steepest descent (the same
        direction the ice flux follows) and **melts out** progressively toward
        each catchment's terminus, building moraine where the ice actually ends
        rather than smearing it across the whole ablation zone. This matters at
        high (sub-km) resolution where individual glacier catchments and termini
        are resolved.

        It is a transport-with-loss on the ice network, solved with the same
        MPI-correct flow-matrix / KSP machinery as the river accumulation:

        .. math::
            L_i = A_i + \sum_{u \to i} w_{ui}\,(1 - f_u)\,L_u,
            \qquad D_i = f_i\,L_i

        where :math:`A_i` is the local abraded volume, :math:`L_i` the till
        load passing through cell :math:`i`, and :math:`f_i \in [0,1]` the
        melt-out fraction — the share of the local ice column lost to ablation
        this step, :math:`f_i = \min(1, \dot{a}_i\,\Delta t / H_i)`. At the ice
        margin (a cell whose steepest-descent receiver is ice-free) :math:`f`
        is forced to 1 so no till leaks onto bare ground, which makes the
        transport exactly conservative: :math:`\sum_i D_i = \sum_i A_i`.

        :arg Vero: per-cell abraded volume (m^3, local array).
        :arg Vtot: total abraded volume over owned nodes (m^3).
        :arg owned: boolean mask of owned (non-ghost) nodes.

        :return: per-cell deposition weight (local array, Σ over owned = 1).
        """
        # Halo-synced ice thickness and meltwater (ghost values are needed for
        # the steepest-descent receivers and the melt-out fraction).
        self.dm.localToGlobal(self.iceHL, self.tmp)
        self.dm.globalToLocal(self.tmp, self.tmpL)
        H = self.tmpL.getArray().copy()
        self.dm.localToGlobal(self.iceMeltL, self.tmp)
        self.dm.globalToLocal(self.tmp, self.tmpL)
        melt = self.tmpL.getArray().copy()

        zbed = self.hLocal.getArray()
        s = zbed + H                                       # ice surface

        # Steepest-descent (single-flow-direction) receivers on the ice surface.
        rcv, _, wght = mfdreceivers(1, self.flowExp, s, self.sealevel, self.gid)
        rcv0 = rcv[:, 0].astype(petsc4py.PETSc.IntType)
        w0 = wght[:, 0].copy()

        thr = 1.0e-2                                       # ice-presence threshold (m)
        ice = H > thr
        # Melt-out fraction f = (ablation thickness this step) / H, clamped.
        meltThick = (melt * self.dt) / np.maximum(self.larea, 1.0e-12)
        f = np.clip(meltThick / np.maximum(H, thr), 0.0, 1.0)
        # Deposit fully (no downstream carry) off-ice and at the ice margin
        # (receiver ice-free), so nothing leaks onto bare ground -> conservative.
        f[~ice] = 1.0
        f[ice & (H[rcv0] <= thr)] = 1.0

        # Downstream-carried weight (1 - f); zero on ghosts/borders.
        outw = w0 * (1.0 - f)
        outw[self.ghostIDs] = 0.0
        if self.flatModel:
            outw[self.idBorders] = 0.0

        # Accumulation matrix (I - W^T) on the ice network, then L = solve(·, A).
        mat = self.iMat.copy()
        indptr = np.arange(0, self.lpoints + 1, dtype=petsc4py.PETSc.IntType)
        nodes = indptr[:-1]
        data = -outw
        data[rcv0 == nodes] = 0.0
        tmpMat = self._matrix_build()
        tmpMat.assemblyBegin()
        tmpMat.setValuesLocalCSR(indptr, rcv0, data)
        tmpMat.assemblyEnd()
        mat.axpy(1.0, tmpMat)
        tmpMat.destroy()
        mat.transpose()

        self.tmpL.setArray(Vero)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self._solve_KSP(True, mat, self.tmp, self.tmp1)
        mat.destroy()
        self.dm.globalToLocal(self.tmp1, self.tmpL)
        L = self.tmpL.getArray().copy()

        # Deposited volume per cell. The transport conserves analytically, but
        # the KSP solve is only accurate to its tolerance, so normalise by the
        # ACTUAL routed total (not Vtot) to make the weight sum to one exactly —
        # mass conservation then matches the melt-weighted path bit-for-bit.
        # The analytical load is non-negative; clamp away sub-tolerance KSP
        # residual noise so the deposition weight can never go slightly negative.
        dep_vol = np.maximum(f * L, 0.0)
        total = MPI.COMM_WORLD.allreduce(
            float(np.sum(dep_vol[owned])), op=MPI.SUM
        )
        if total <= 0.0:
            return np.zeros(self.lpoints, dtype=np.float64)
        return dep_vol / total
