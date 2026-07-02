import os
import petsc4py
import numpy as np

from time import process_time
from mpi4py import MPI

from gospl.tools.constants import ICE_COVER_MIN

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import jacobiancoeff

MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()


class GWMesh(object):
    r"""
    Water table (groundwater) + generic duricrust — opt-in near-surface hydrology
    and chemical armoring. See ``docs/DESIGN_WATERTABLE_DURICRUST.md``.

    An implicit (backward-Euler) Dupuit-Boussinesq head solve on the DMPlex drives
    a generic capillary-fringe duricrust that armors erodibility. Enabled by the
    YAML ``groundwater:`` block (``self.gwOn``); when off, every path here is a
    no-op and goSPL is byte-identical to a run without it.

    **Phase 0 (this file so far):** state allocation only — the persistent Vecs,
    the per-node numpy state and the cached-solver handles are created (gated on
    ``gwOn``); the recharge, head solve, duricrust ODE and K-armoring land in
    later phases. All new Vecs are registered in ``destroy_DMPlex``.
    """

    def __init__(self, *args, **kwargs):
        """
        Initialise the groundwater / duricrust state. Allocated only when
        ``self.gwOn`` (set by ``inputparser._readGroundwater``); otherwise the
        cached-solver handles are still set to ``None`` so ``destroy_DMPlex`` and
        any ``getattr`` guards are safe.
        """

        # Cached elliptic head solver + operator (built lazily in a later phase).
        # Set unconditionally so destroy_DMPlex / guards never hit a missing attr.
        self._gwMat = None
        self._ksp_gw = None

        if getattr(self, "gwOn", False):
            # --- PETSc state (persistent, halo-synced; in destroy_DMPlex) ---
            # Water-table head h (elevation of the saturated surface, m).
            self.headL = self.hLocal.duplicate()
            self.headG = self.hGlobal.duplicate()
            # Duricrust thickness duriH (m).
            self.duriHL = self.hLocal.duplicate()
            self.duriHG = self.hGlobal.duplicate()
            # Net recharge R this step (m/yr) — diagnostic/output.
            self.rechargeL = self.hLocal.duplicate()
            # Seepage return to rivers (m^3/yr) — only used when conserve_baseflow.
            self.baseflowL = self.hLocal.duplicate()

            # Resolve a per-vertex aquifer_base map (done here — needs locIDs).
            # `[file, key]` -> file + ".npz" subset to the local partition; the
            # z_bed depth below the surface (m) then varies in space (regolith /
            # weathering-front depth). Scalar / 'from_soil' left untouched.
            basemap = getattr(self, "_gwAquiferBaseMap", None)
            if basemap is not None:
                data = np.load(basemap[0] + ".npz")
                self.gwAquiferBase = data[basemap[1]][self.locIDs].astype(
                    np.float64
                )

            # Seed the head at the aquifer BASE (a dry start) so recharge fills
            # it UP to the steady water table. Seeding at the surface is wrong:
            # recharge would push h above z, the seepage clip would pin it there,
            # and it could never drain below the surface (spurious full
            # saturation). For a scalar / map aquifer_base, z_bed = z − base;
            # for `from_soil` (Phase 5, a string) fall back to the surface.
            base = self.gwAquiferBase
            if isinstance(base, str):
                self.hGlobal.copy(result=self.headG)
                self.hLocal.copy(result=self.headL)
            else:
                self.headL.setArray(self.hLocal.getArray() - base)
                self.dm.localToGlobal(self.headL, self.headG)
                self.dm.globalToLocal(self.headG, self.headL)
            self.duriHL.set(0.0)
            self.duriHG.set(0.0)
            self.rechargeL.set(0.0)
            self.baseflowL.set(0.0)

            # --- numpy state (rank-local, no halo) ---
            self.wtDepth = np.zeros(self.lpoints, dtype=np.float64)     # z - h (m)
            self.duriF = np.zeros(self.lpoints, dtype=np.float64)       # induration 0..1
            self.duriKarmor = np.ones(self.lpoints, dtype=np.float64)   # K multiplier (<=1)
            self.gwSeepIDs = np.zeros(0, dtype=np.int64)                # Dirichlet seepage nodes

            # Resolve a per-vertex infiltration map (done here — needs locIDs,
            # unavailable at parse time). `[file, key]` -> file + ".npz", subset
            # to the local partition. Scalar `gwInfiltration` is left untouched.
            infilmap = getattr(self, "_gwInfilMap", None)
            if infilmap is not None:
                data = np.load(infilmap[0] + ".npz")
                self.gwInfiltration = data[infilmap[1]][self.locIDs].astype(
                    np.float64
                )

        return

    def updateGroundwater(self):
        """
        Per-step groundwater / duricrust update.

        **Phase 1:** net recharge only — ``R = f_infil * max(0, rain - evap)``
        (m/yr), held at 0 under standing water (marine ``seaID`` or a ponded
        continental lake ``pitIDs>-1 & lFill>hl``), where the water table is
        pinned to the surface. Stored in ``self.rechargeL`` for output. The
        implicit head solve, duricrust ODE and K-armoring are added in later
        phases — nothing consumes the recharge yet.

        No-op when ``gwOn`` is off. Purely local (per-node) — no collective.
        """
        if not getattr(self, "gwOn", False):
            return

        t0 = process_time()
        rain = self.rainVal
        evap = getattr(self, "evapVal", None)
        net = rain if evap is None else (rain - evap)
        R = self.gwInfiltration * np.maximum(0.0, net)

        # No rain-recharge where the surface is not subaerial land:
        #  - standing water (marine `seaID`, or a ponded continental lake) — the
        #    head is pinned to the surface there;
        #  - ice-covered land — precipitation falls as snow/ice and does not
        #    infiltrate the ground (parallels the soil ice-freeze gate). Subglacial
        #    meltwater recharge is a future refinement (the ice model has
        #    `iceMeltRiverL`). See DESIGN_WATERTABLE_DURICRUST.md §3.
        sub = np.zeros(self.lpoints, dtype=bool)
        sub[self.seaID] = True
        pitIDs = getattr(self, "pitIDs", None)
        lFill = getattr(self, "lFill", None)
        if pitIDs is not None and lFill is not None:
            sub |= (pitIDs > -1) & (lFill > self.hLocal.getArray())
        if getattr(self, "iceOn", False):
            iceHL = getattr(self, "iceHL", None)
            if iceHL is not None:
                sub |= iceHL.getArray() > ICE_COVER_MIN
        R = np.where(sub, 0.0, R)

        self.rechargeL.setArray(R)

        # Phase 2: solve the implicit water-table head from this recharge.
        self._solveHead()

        if MPIrank == 0 and self.verbose:
            print(
                "Update Groundwater Table (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )
        return

    def _makeGWKSP(self):
        """
        Cached KSP for the elliptic head solve: fgmres + **hypre BoomerAMG**
        (algebraic multigrid). The head operator ``I + (Δt/S)·L(T)`` is a stiff,
        high-condition 2-D elliptic operator; block-Jacobi/ILU does NOT converge
        on it (it stalls at ``DIVERGED_ITS`` and returns a garbage iterate that
        drifts to the surface) — an elliptic solve needs a multigrid PC. Same
        rationale as the flexure biharmonic wanting a strong solver. Options
        prefix ``gw_`` scopes any override (env ``-gw_pc_type ...``).
        """
        ksp = petsc4py.PETSc.KSP().create(petsc4py.PETSc.COMM_WORLD)
        ksp.setType("fgmres")
        ksp.getPC().setType("hypre")
        ksp.setTolerances(rtol=1.0e-10, max_it=500)
        ksp.setInitialGuessNonzero(True)
        ksp.setOptionsPrefix("gw_")
        ksp.setFromOptions()
        return ksp

    def _solveHead(self):
        r"""
        Implicit (backward-Euler) Dupuit-Boussinesq water-table head solve
        (DESIGN_WATERTABLE_DURICRUST.md §2/§3). One elliptic solve per step:

        .. math::
            (I + \tfrac{\Delta t}{S}\,L(T))\,h = h_{old} + \tfrac{\Delta t}{S}\,R

        with transmissivity :math:`T = K_h\max(h - z_{bed}, b_{min})` (unconfined),
        :math:`L(T)=-\nabla\cdot(T\nabla)` the FV neg-Laplacian (from
        ``jacobiancoeff``, area-normalised — the same operator the marine Picard
        solver uses). Reduces to the steady elliptic limit as :math:`\Delta t/S`
        grows, so it is quasi-steady at century steps.

        - **Unconfined non-linearity** ``T(h)``: ``gwPicardIts`` Picard sweeps
          (lag ``T`` on the previous iterate; operator rebuilt per sweep).
        - **Seepage free boundary** ``h <= z``: Dirichlet ``h = z`` at the base
          seepage set (marine ``seaID`` + ponded lakes + open ``outletIDs``) via
          ``zeroRowsLocal``; after each Picard block, clip ``h = min(h, z)`` and
          add any newly-over-topped cell to the Dirichlet set, up to
          ``gwSeepagePasses`` passes (break on an ``Allreduce``'d new-seepage
          count — collective-safe). ``aquifer_base`` is prescribed here (the
          ``from_soil`` / ``lHbed`` tie is Phase 5).

        Scratch: uses ``self.tmp`` (global rhs); leaves it defined. Writes
        ``self.headL/headG`` and ``self.wtDepth = z - h``.
        """
        z = self.hLocal.getArray()
        # z_bed = surface − aquifer_base (scalar or per-vertex map, numpy-safe).
        zbed = z - self.gwAquiferBase
        bmin = float(self.gwMinSatThick)
        gwdt = self.dt / float(self.gwSpecificYield)          # Δt/S
        R = self.rechargeL.getArray()
        inIDs = self.inIDs
        zeroKp = np.zeros(self.lpoints, dtype=np.float64)
        IntType = petsc4py.PETSc.IntType

        # Base seepage (Dirichlet h = z): standing water + open outlets.
        seep = np.zeros(self.lpoints, dtype=bool)
        seep[self.seaID] = True
        pitIDs = getattr(self, "pitIDs", None)
        lFill = getattr(self, "lFill", None)
        if pitIDs is not None and lFill is not None:
            seep |= (pitIDs > -1) & (lFill > z)
        outletIDs = getattr(self, "outletIDs", None)
        if outletIDs is not None and len(outletIDs) > 0:
            seep[outletIDs] = True

        if self._ksp_gw is None:
            self._ksp_gw = self._makeGWKSP()
        ksp = self._ksp_gw

        self.dm.globalToLocal(self.headG, self.headL)
        hloc = self.headL.getArray().copy()
        hold = hloc.copy()

        npass = 0
        seep_converged = False
        for _ in range(int(self.gwSeepagePasses)):
            npass += 1
            for _ in range(int(self.gwPicardIts)):
                T = self.gwKsat * np.maximum(hloc - zbed, bmin)   # transmissivity (m²/yr)
                coeffs = gwdt * jacobiancoeff(hloc, T, zeroKp)     # (Δt/S)·L(T)
                coeffs[:, 0] += 1.0                                # I + (Δt/S)·L(T)
                M = self._assembleDiffMatCSR(coeffs)
                owned_seep = np.where(seep & (inIDs == 1))[0].astype(IntType)
                M.zeroRowsLocal(owned_seep, diag=1.0)             # collective; h = rhs there

                rhs = hold + gwdt * R
                rhs[seep] = z[seep]                               # pin h = z at seepage
                self.headL.setArray(rhs)
                self.dm.localToGlobal(self.headL, self.tmp)       # rhs (global)
                self.headL.setArray(hloc)
                self.dm.localToGlobal(self.headL, self.headG)     # nonzero guess
                ksp.setOperators(M, M)
                ksp.solve(self.tmp, self.headG)
                M.destroy()
                self.dm.globalToLocal(self.headG, self.headL)
                hloc = self.headL.getArray().copy()

            # Free boundaries: seepage (upper, h <= surface — discover new
            # seepage cells) and the dry-aquifer floor (lower, h >= z_bed — the
            # head cannot drop below the impermeable base).
            over = (hloc > z) & (~seep)
            n_new = MPI.COMM_WORLD.allreduce(
                int((over & (inIDs == 1)).sum()), op=MPI.SUM
            )
            seep |= over
            hloc = np.clip(hloc, zbed, z)
            if n_new == 0:
                seep_converged = True
                break

        if MPIrank == 0 and self.verbose:
            # `reason > 0` = KSP converged; a non-converged elliptic solve (e.g.
            # a too-weak PC) silently corrupts the water table, so surface it.
            reason = int(ksp.getConvergedReason())
            print(
                "[gw] head solve: %d seepage pass(es)%s, last KSP reason %d%s"
                % (
                    npass,
                    "" if seep_converged else " (seepage set still growing)",
                    reason,
                    "" if reason > 0 else " -- WARNING: KSP did not converge",
                ),
                flush=True,
            )

        self.headL.setArray(hloc)
        self.dm.localToGlobal(self.headL, self.headG)
        self.wtDepth = z - hloc
        return
