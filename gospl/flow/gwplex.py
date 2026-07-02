import os
import petsc4py
import numpy as np

from mpi4py import MPI

from gospl.tools.constants import ICE_COVER_MIN

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

            # Seed head to the current surface (a valid starting water table:
            # h = z, i.e. fully saturated / at the surface) so the first solve
            # has a bounded guess.
            self.hGlobal.copy(result=self.headG)
            self.hLocal.copy(result=self.headL)
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
        return
