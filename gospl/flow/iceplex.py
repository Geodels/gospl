import os
import gc
import sys
import scipy
import petsc4py
import numpy as np
import numpy_indexed as npi

from mpi4py import MPI
from time import process_time

from gospl.tools.constants import BOUNDARY_FLOW_SENTINEL, MISSING_DATA_SENTINEL

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import mfdreceivers
    from gospl._fortran import epsfill
    from gospl._fortran import ice_flux

# Ice density (kg/m^3), used for the SIA deformation coefficient and loading.
RHO_ICE = 910.0

MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()


class IceMesh(object):
    """
    This class calculates **ice flow acculation** based on a multiple flow direction paths (MFD).

    Useful links for improvements:
    - `Braun et al. (1999) <https://doi.org/10.3189/172756499781821797>`_
    - `Herman & Braun (2008) <https://doi.org/10.1029/2007JF000807>`_
    - `Tomkin (2009) <https://doi.org/10.1016/j.geomorph.2008.04.021>`_
    - `Egholm (2012) <https://doi.org/10.1016/j.geomorph.2011.12.019>`_
    - `Deal & Prasicek (2020) <https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2020GL089263>`_
    - `Hergarten (2021) <http://dx.doi.org/10.5194/esurf-2021-1>`_
    - `Liebl et al. (2023) <https://gmd.copernicus.org/articles/16/1315/2023/>`_
    """

    def __init__(self, *args, **kwargs):
        """
        Initialisation of the `IceMesh` class.

        This method initializes the ice-related fields (ice height, flow accumulation, and flexural response)
        based on the configuration flags `iceOn` and `flexOn`. Memory is allocated for these fields.
        """

        if self.iceOn:
            self.iceHL = self.hLocal.duplicate()
            self.iceFAG = self.hGlobal.duplicate()
            self.iceFAL = self.hLocal.duplicate()
            # Local meltwater field captured at the end of iceAccumulation
            # and re-injected into the river FA source term in flowplex so
            # glacial discharge is not lost downstream.
            self.iceMeltL = self.hLocal.duplicate()
            if self.flexOn:
                self.iceFlex = self.hLocal.duplicate()
            self.iceHL.set(0.0)
            self.iceMeltL.set(0.0)

        return

    def _matrixIceFlow(self, dir_ice=1):
        """
        Compute Flow Direction Matrix for Ice Flow.

        This function calculates the flow direction matrix for ice flow using a multiple flow direction (MFD) algorithm.
        It fills in elevation data (using `epsfill`) and constructs the matrix for flow direction and weighting.

        Parameters:
        -----------
        dir_ice : int, optional
            Number of flow directions to consider. Defaults to 1.
        """

        # The filled + eps is done on the global grid!
        hl = self.hLocal.getArray().copy()
        minh = self.hGlobal.min()[1] + 0.1
        minh = max(minh, self.sealevel)
        if self.flatModel:
            hl[self.idBorders] = BOUNDARY_FLOW_SENTINEL
        fillz = np.zeros(self.mpoints, dtype=np.float64) + MISSING_DATA_SENTINEL
        fillz[self.locIDs] = hl
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, fillz, op=MPI.MAX)
        if MPIrank == 0:
            fillz = epsfill(minh, fillz)
        # Broadcast filled elevation data across processors
        fillEPS = MPI.COMM_WORLD.bcast(fillz, root=0)
        fillz = fillEPS[self.locIDs]

        # Calculate receivers and weights for the flow direction matrix.
        # `self.gid` is the per-node global ID, passed in for deterministic
        # exact-tie-break on slope (see fortran/functions.F90:mfdreceivers).
        rcv, _, wght = mfdreceivers(
            dir_ice, 1.0, fillz, BOUNDARY_FLOW_SENTINEL, self.gid,
        )

        # Handle borders for flat models
        if self.flatModel:
            rcv[self.idBorders, :] = np.tile(self.idBorders, (dir_ice, 1)).T
            wght[self.idBorders, :] = 0.0

        self.iceRcv = rcv.copy()
        self.iceWght = wght.copy()

        # Create and assemble the flow direction matrix
        self.iceMat = self.iMat.copy()
        indptr = np.arange(0, self.lpoints + 1, dtype=petsc4py.PETSc.IntType)
        nodes = indptr[:-1]

        for k in range(dir_ice):
            # Flow direction matrix for a specific direction
            tmpMat = self._matrix_build()
            data = wght[:, k].copy()
            data[rcv[:, k].astype(petsc4py.PETSc.IntType) == nodes] = 0.0
            tmpMat.assemblyBegin()
            tmpMat.setValuesLocalCSR(
                indptr,
                rcv[:, k].astype(petsc4py.PETSc.IntType),
                data,
            )
            tmpMat.assemblyEnd()

            # Accumulate weights for each direction
            self.iceMat.axpy(-1.0, tmpMat)
            tmpMat.destroy()

        if self.memclear:
            del data, indptr, nodes
            del hl, fillz, fillEPS, rcv, wght
            gc.collect()

        # Transpose the flow matrix for further computations
        self.iceMat.transpose()

        return

    def _iceFlowSIA(self, elaH, iceH, iceT):
        r"""
        Shallow-Ice-Approximation ice-thickness evolution (Phase 1).

        Prognostic: evolves the **persistent** ice thickness ``self.iceHL`` over
        the goSPL step by integrating

        .. math:: \partial H/\partial t = \dot m - \nabla\!\cdot q,\quad q=(u_d+u_s)H,\ s=z_{bed}+H

        with the SIA flux divergence from the Glen's-law ``ice_flux`` kernel
        (deformation `u_d` + sliding `u_s`) on the **total surface** (bed + H).
        Surface mass balance `ṁ` is the ELA model (no `meltfac` — that is an
        MFD numerical sink amplifier, not physical).

        .. important::
            This is the **explicit, CFL-subcycled reference scheme** — the
            validation oracle (see ``docs/DESIGN_ICE_SHEET.md`` §3). It is
            correct but, at goSPL's km / 10²–10⁴ yr scale, can hit
            ``sia_max_substeps``; the **implicit (semi-implicit SNES)** solve is
            the production scheme (Phase 1b). Coefficients (`ad`, `as`, `glen`)
            are uniform (the kernel's signature); spatially-variable sliding is
            a later extension.
        """

        n = self.sia_glen
        # Deformation coefficient ad = 2 A (rho_ice g)^n / (n+2) (SIA flux form).
        ad = 2.0 * self.sia_Aglen * (RHO_ICE * self.gravity) ** n / (n + 2.0)
        as_ = self.sia_slide

        zbed = self.hLocal.getArray().copy()              # bed elevation (fixed this step)
        rainV = self.rainVal.copy()                       # precip rate (m/yr)
        ramp = (zbed - elaH) / (iceH - elaH)
        ramp[ramp > 1.0] = 1.0
        mdot = rainV * ramp                               # surface mass balance (m ice/yr)

        H = self.iceHL.getArray().copy()                  # persistent ice thickness

        # Precompute neighbour offsets/distances once (mesh is fixed).
        if not hasattr(self, "_iceNgbSafe"):
            ngb = self.FVmesh_ngbID
            valid = ngb >= 0
            safe = np.where(valid, ngb, 0)
            d = self.lcoords[safe] - self.lcoords[:, None, :]
            dist = np.sqrt((d * d).sum(axis=2))
            dist[~valid] = np.inf
            finite = dist[valid]
            self._iceNgbSafe = safe
            self._iceNgbDist = dist
            self._iceMinDist2 = float(finite.min() ** 2) if finite.size else 1.0
        safe = self._iceNgbSafe
        dist = self._iceNgbDist
        nm1 = n - 1.0

        # CFL-subcycled explicit integration over the goSPL timestep.
        t_done = 0.0
        sub = 0
        while t_done < self.dt - 1.0e-9:
            # Halo-sync H so ice_flux sees current ghost-node values (MPI).
            self.tmpL.setArray(H)
            self.dm.localToGlobal(self.tmpL, self.tmp)
            self.dm.globalToLocal(self.tmp, self.tmpL)
            H = self.tmpL.getArray().copy()

            div = ice_flux(self.lpoints, H, zbed, ad, as_, n)   # ∇·q per cell

            # Stable substep from the SIA diffusivity D = (ad H^{n+1}+as H^{n-1})·H·|∇s|^{n-1}.
            s = zbed + H
            slope = np.abs(s[:, None] - s[safe]) / dist
            gradS = slope.max(axis=1)
            Dcell = (ad * H ** (n + 1.0) + as_ * np.maximum(H, 0.0) ** nm1) * H * gradS ** nm1
            Dmax = MPI.COMM_WORLD.allreduce(
                float(np.max(Dcell)) if Dcell.size else 0.0, op=MPI.MAX
            )
            if Dmax > 0.0:
                dt_sub = min(
                    self.dt - t_done,
                    self.sia_cfl * self._iceMinDist2 / (2.0 * Dmax),
                )
            else:
                dt_sub = self.dt - t_done

            H = H + dt_sub * (mdot - div)
            H[H < 0.0] = 0.0                              # free boundary: no negative ice
            t_done += dt_sub
            sub += 1
            if sub >= self.sia_max_substeps:
                if MPIrank == 0:
                    print(
                        "SIA ice solve hit max_substeps=%d (%.1f%% of dt done); "
                        "consider the implicit solver."
                        % (self.sia_max_substeps, 100.0 * t_done / self.dt),
                        flush=True,
                    )
                break

        # Glacier terminus clamp (below it, ice is removed).
        H[zbed < iceT] = 0.0

        self.iceHL.setArray(H)

        # Meltwater for the river coupling: ablation water volume where ice is
        # present (m^3/yr), consumed by flowplex.flowAccumulation. No meltfac.
        melt_local = np.maximum(-mdot, 0.0) * self.larea * (H > 1.0e-2)
        self.iceMeltL.setArray(melt_local)

        # The MFD "ice flow accumulation" is not meaningful under SIA; zero it
        # (SIA-mode glacial erosion uses the basal velocity, Phase 3).
        self.iceFAL.set(0.0)
        self.iceFAG.set(0.0)

        return

    def iceAccumulation(self):
        """
        Main Ice Accumulation Calculation.

        This method calculates ice accumulation based on the elevation, equilibrium-line altitude (ELA), and terminus positions. It integrates the flow direction matrix and computes ice flow across the landscape.

        The method also accounts for flexural responses if enabled.
        """

        ti = process_time()

        # Copy ice height to flexural parameter if flexure modeling is active
        if self.flexOn:
            self.iceHL.copy(result=self.iceFlex)

        # Get dynamic properties such as ice cap altitude and equilibrium-line altitude
        iceH = self.iceH(self.tNow)  # Ice Cap Altitude
        elaH = self.elaH(self.tNow)  # Equilibrium-Line Altitude
        iceT = self.iceT(self.tNow)  # Glacier terminus

        # If maximum elevation is below ELA, no ice accumulation occurs
        max_elev = self.hGlobal.max()[1]
        if max_elev < elaH:
            self.iceHL.set(0.)
            self.iceFAL.set(0.)
            self.iceFAG.set(0.)
            self.iceMeltL.set(0.)
            if self.flexOn:
                self.iceFlex.set(0.)
            return

        # Degenerate config: ice-cap altitude must lie above the ELA, otherwise
        # the (hl - elaH) / (iceH - elaH) ratio below produces NaN/Inf.
        if iceH <= elaH:
            self.iceHL.set(0.)
            self.iceFAL.set(0.)
            self.iceFAG.set(0.)
            self.iceMeltL.set(0.)
            if self.flexOn:
                self.iceFlex.set(0.)
            return

        # Dual ice flow models: SIA dynamics (opt-in) vs the MFD routing proxy.
        if self.iceSIA:
            self._iceFlowSIA(elaH, iceH, iceT)
            if self.flexOn and self.tNow == self.tStart:
                self.iceHL.copy(result=self.iceFlex)
            if MPIrank == 0 and self.verbose:
                print(
                    "Glaciers Accumulation - SIA (%0.02f seconds)"
                    % (process_time() - ti),
                    flush=True,
                )
            return

        # Compute the flow direction matrix for ice
        self._matrixIceFlow(self.iceDir)

        # Ice accumulation calculation
        hl = self.hLocal.getArray().copy()
        rainA = self.bL.getArray().copy()
        tmp = (hl - elaH) / (iceH - elaH)
        tmp[tmp > 1.] = 1.0
        tmp[tmp < 0.] *= self.meltfac

        # Calculate accumulation rates
        iceA = np.multiply(rainA, tmp)
        self.tmpL.setArray(iceA)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self._solve_KSP(True, self.iceMat, self.tmp, self.iceFAG)
        self.dm.globalToLocal(self.iceFAG, self.iceFAL)
        ice_clamped = self.iceFAL.getArray()
        ice_clamped[ice_clamped < 0.] = 0.0
        ice_clamped[hl < iceT] = 0.0
        self.iceFAL.setArray(ice_clamped)
        self.dm.localToGlobal(self.iceFAL, self.iceFAG)

        # Capture glacial meltwater for the river FA. Below ELA the local
        # ablation rate is `|hl - elaH|/(iceH - elaH)` of the precipitation;
        # we gate by ice presence (clamped, pre-smoothing iceFAL > 0) so we
        # only release water where ice actually reached. `meltfac` is
        # deliberately NOT applied here because it is a numerical
        # sink-amplifier inside the implicit ice solver, not a physical
        # melt multiplier. This is consumed by flowplex.flowAccumulation
        # and added back to `rainA` before the river FA is solved.
        melt_local = np.zeros(self.lpoints, dtype=np.float64)
        below_ela = hl < elaH
        if below_ela.any():
            ablation = np.zeros(self.lpoints, dtype=np.float64)
            ablation[below_ela] = np.clip(
                -(hl[below_ela] - elaH) / (iceH - elaH), 0.0, 1.0
            )
            has_ice = ice_clamped > 1.0e-8  # TODO-REFACTOR: value matches DISCHARGE_FLOOR but distinct role (ice-presence threshold for melt capture); do not replace
            melt_local = ablation * rainA * has_ice
        self.iceMeltL.setArray(melt_local)

        # Smooth and diffuse glacier flow accumulation.
        # NB: the linear-diffusion smoothing below is NOT mass-conservative
        # (it spreads the field for visualisation and erosion-driver
        # robustness). The output `iceFA` therefore differs slightly from
        # the strict accumulated source. The meltwater field above was
        # captured BEFORE smoothing, so the river re-injection is based on
        # the true (un-smeared) ice presence.
        self.iceFAG.copy(result=self.tmp1)
        smthIce = self._hillSlope(smooth=1)
        smthIce[smthIce < 0.] = 0.
        self.iceFAL.setArray(smthIce)
        self.dm.localToGlobal(self.iceFAL, self.iceFAG)

        # Ice thickness from Bahr-style width-area scaling. Computed
        # whenever ice is on (used by both flexure and the `iceH` output);
        # the flex copy stays gated behind `flexOn`.
        tmp = self.icewe * self.icewf * smthIce**0.3
        tmp[tmp < 1.e-1] = 0.  # TODO-REFACTOR: value matches BEDROCK_EXPOSED but distinct role (minimum ice thickness for visualization); do not replace
        self.tmpL.setArray(tmp)
        self.dm.localToGlobal(self.tmpL, self.tmp1)
        tmp = self._hillSlope(smooth=1)
        self.iceHL.setArray(tmp)

        if self.flexOn:
            # If simulation starts then set the ice flex variable to the initial glacier thickness
            if self.tNow == self.tStart:
                self.iceHL.copy(result=self.iceFlex)

        if MPIrank == 0 and self.verbose:
            print(
                "Glaciers Accumulation (%0.02f seconds)" % (process_time() - ti),
                flush=True,
            )

        return
