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
    from gospl._fortran import ice_velocity

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
            # Basal sliding speed (m/yr) from the SIA solve; the driver for the
            # velocity-based glacial abrasion law (Phase 3). Zero under the MFD
            # proxy. Registered in destroy_DMPlex.
            self.iceUbL = self.hLocal.duplicate()
            self.iceHL.set(0.0)
            self.iceMeltL.set(0.0)
            self.iceUbL.set(0.0)
            # Glacial-till mass-balance diagnostics (m^3, owned-node running
            # totals reduced by the till-conservation test).
            self._tillEroded = 0.0
            self._tillDeposited = 0.0
            # Cached SIA implicit-solver handles (created lazily on first use;
            # destroyed in destroy_DMPlex).
            self._snes_ice = None
            self._snes_ice_f = None
            self._snes_ice_x = None

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
        ramp = (zbed - elaH) / (iceH - elaH)
        ramp[ramp > 1.0] = 1.0
        mdot = self.rainVal * ramp                        # surface mass balance (m ice/yr)
        return n, ad, as_, zbed, mdot

    def _iceSIAFinalize(self, H, zbed, mdot, iceT, as_, n):
        """
        Common post-solve steps for both SIA solvers: terminus clamp, store the
        ice thickness, the basal sliding speed (for Phase-3 abrasion), capture
        ablation meltwater (m^3/yr) for the river coupling, and zero the MFD
        ice-flow field (SIA-mode erosion uses the basal velocity, not it).
        """
        H = np.maximum(H, 0.0)
        H[zbed < iceT] = 0.0
        self.iceHL.setArray(H)
        # Basal sliding speed from the converged thickness (steepest-descent
        # SIA velocity); zero where there is no ice.
        ub = ice_velocity(self.lpoints, H, zbed, as_, n)
        ub[H <= 1.0e-2] = 0.0
        self.iceUbL.setArray(ub)
        melt_local = np.maximum(-mdot, 0.0) * self.larea * (H > 1.0e-2)
        self.iceMeltL.setArray(melt_local)
        self.iceFAL.set(0.0)
        self.iceFAG.set(0.0)
        return

    def _form_residual_ice(self, snes, H, F):
        """
        SNES residual for the implicit SIA thickness solve (mirrors
        ``hillslope._form_residual_nl_hillslope``):

            F(H) = H − H_old − Δt (ṁ − ∇·q(H))

        with `∇·q` from the ``ice_flux`` kernel on the clamped (≥0) thickness so
        the flux stays physical during the nonlinear iterations.
        """
        self.dm.globalToLocal(H, self.tmpL)
        Ha = self.tmpL.getArray()
        Hc = np.maximum(Ha, 0.0)
        div = ice_flux(
            self.lpoints, Hc, self._sia_zbed, self._sia_ad, self._sia_as, self._sia_n
        )
        res = Ha - self._sia_Hold - self.dt * (self._sia_mdot - div)
        F.setArray(res[self.glIDs])
        return

    def _iceFlowSIAImplicit(self, elaH, iceH, iceT):
        r"""
        **Production** SIA ice-thickness solve (Phase 1b): implicit
        (semi-implicit) integration of ``∂H/∂t = ṁ − ∇·q`` via a cached PETSc
        SNES, reusing the non-linear hillslope pattern (``ngmres``, CG/HYPRE
        KSP, Jacobian-free). Unconditionally stable in `dt`, so it takes the
        full goSPL timestep — adequate for goSPL's km / 10²–10⁴ yr scale where
        explicit subcycling is prohibitive. Free boundary `H≥0` by post-clamp;
        validated against the explicit reference (`_iceFlowSIAExplicit`) and the
        analytical SIA dome.
        """
        n, ad, as_, zbed, mdot = self._iceSIAParams(elaH, iceH)
        # Residual context (read by _form_residual_ice).
        self._sia_n = n
        self._sia_ad = ad
        self._sia_as = as_
        self._sia_zbed = zbed
        self._sia_mdot = mdot
        self._sia_Hold = self.iceHL.getArray().copy()

        if self._snes_ice is None:
            snes = petsc4py.PETSc.SNES().create(comm=petsc4py.PETSc.COMM_WORLD)
            snes.setTolerances(rtol=1.0e-6, atol=1.0e-6, max_it=100)
            self._snes_ice_f = self.hGlobal.duplicate()
            snes.setFunction(self._form_residual_ice, self._snes_ice_f)
            snes.setType("ngmres")
            ksp = snes.getKSP()
            ksp.setType("cg")
            pc = ksp.getPC()
            pc.setType("hypre")
            petsc4py.PETSc.Options()["pc_hypre_type"] = "boomeramg"
            pc.setFromOptions()
            ksp.setTolerances(rtol=1.0e-6)
            self._snes_ice_x = self.hGlobal.duplicate()
            self._snes_ice = snes

        snes = self._snes_ice
        x = self._snes_ice_x
        # Initial guess = previous ice thickness.
        self.tmpL.setArray(self._sia_Hold)
        self.dm.localToGlobal(self.tmpL, x)
        snes.solve(None, x)
        r = snes.getConvergedReason()
        if r < 0 and MPIrank == 0:
            print(
                "SIA implicit ice SNES failed to converge after %d iterations "
                "(reason %d)" % (snes.getIterationNumber(), r),
                flush=True,
            )

        self.dm.globalToLocal(x, self.tmpL)
        self._iceSIAFinalize(self.tmpL.getArray().copy(), zbed, mdot, iceT, as_, n)

        return

    def _iceFlowSIAExplicit(self, elaH, iceH, iceT):
        r"""
        Explicit, CFL-subcycled SIA reference scheme — the **validation oracle**
        (see ``docs/DESIGN_ICE_SHEET.md`` §3). Correct but, at goSPL's km /
        10²–10⁴ yr scale, can hit ``sia_max_substeps``; use the implicit solver
        for production. Integrates ``∂H/∂t = ṁ − ∇·q`` with `ice_flux` on the
        total surface (bed + H), substepping under the SIA diffusivity CFL.
        """

        n, ad, as_, zbed, mdot = self._iceSIAParams(elaH, iceH)
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

        self._iceSIAFinalize(H, zbed, mdot, iceT, as_, n)

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
            if self.sia_solver == "explicit":
                self._iceFlowSIAExplicit(elaH, iceH, iceT)
            else:
                self._iceFlowSIAImplicit(elaH, iceH, iceT)
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

    def glacialTill(self):
        r"""
        Glacial till — production, transport and deposition (Phase 4, SIA only).

        Glacial abrasion (:math:`E_g = K_g |u_b|^l`) erodes rock under sliding
        ice; the abraded material becomes **till** carried by the ice and
        deposited where the ice **melts out** — the ablation zone / terminal
        moraine. The bed is therefore **lowered** under fast ice and **raised**
        in the ablation zone, conserving rock volume exactly (the net bed change
        is zero by construction): a transport, not a source/sink.

        Conservation is structural: the total abraded volume ``Vtot`` is
        redistributed to the ablation cells weighted by the meltwater rate
        (``iceMeltL``), so ``Σ deposited == Vtot``. Guarded by a dedicated
        till-conservation test (the dual-lithology lesson: a volume the total
        budget can't see needs its own guard).

        Active only with ``flow_model: sia``, ``till.on`` and ``Kg > 0``; a
        no-op otherwise, so existing behaviour is unchanged.

        Simplifications (v1, documented): till deposits across the ablation zone
        weighted by melt rather than routed per ice-catchment to each terminus;
        and it operates on the bulk bed (cumED/elevation), not yet the
        stratigraphic pile — couple to stratigraphy / dual lithology later.
        """
        if not (self.iceOn and self.iceSIA and self.ice_till_on) or self.ice_Kg <= 0.0:
            return

        ti = process_time()
        ub = self.iceUbL.getArray()
        Eg = self.ice_Kg * np.power(np.maximum(ub, 0.0), self.ice_abr_l)  # m/yr
        Eg[self.seaID] = 0.0
        Vero = Eg * self.dt * self.larea                  # abraded volume (m^3) per cell
        owned = self.inIDs == 1
        Vtot = MPI.COMM_WORLD.allreduce(float(np.sum(Vero[owned])), op=MPI.SUM)

        melt = self.iceMeltL.getArray()                   # ablation water volume (m^3/yr)
        Wtot = MPI.COMM_WORLD.allreduce(float(np.sum(melt[owned])), op=MPI.SUM)

        # No abrasion or no ablation zone to receive the till -> nothing happens
        # (skipping the erosion too keeps it conservative: no orphaned till).
        if Vtot <= 0.0 or Wtot <= 0.0:
            return

        # Bed change (m): erode at abrasion cells, deposit the till where ice
        # melts, weighted by melt. Σ(dz·area) == 0 (rock moved, not created).
        dz = -Eg * self.dt
        dz_dep = (Vtot * melt / Wtot) / self.larea
        dz = dz + dz_dep

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
