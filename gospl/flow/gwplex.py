import os
import petsc4py
import numpy as np

from mpi4py import MPI

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

        return
