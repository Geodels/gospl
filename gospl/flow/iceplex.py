import os
import petsc4py
import numpy as np

from mpi4py import MPI
from time import process_time

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import ice_flux
    from gospl._fortran import ice_velocity

# Ice density (kg/m^3), used for the SIA deformation coefficient and loading.
RHO_ICE = 910.0

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
    The thickness is integrated **implicitly** with a cached Jacobian-free PETSc
    ``SNES`` (unconditionally stable, so it takes the full goSPL timestep —
    adequate at the km / :math:`10^2`–:math:`10^4` yr scale where explicit
    sub-cycling is prohibitive).

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
        Common post-solve steps for the SIA solve: terminus clamp, store the
        ice thickness and the basal sliding speed (the abrasion driver), and
        capture ablation meltwater (m^3/yr) for the river coupling.
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

    def _iceFlowSIA(self, elaH, iceH, iceT):
        r"""
        SIA ice-thickness solve: implicit (semi-implicit) integration of
        ``∂H/∂t = ṁ − ∇·q`` via a cached PETSc SNES, reusing the non-linear
        hillslope pattern (``ngmres``, CG/HYPRE KSP, Jacobian-free).
        Unconditionally stable in `dt`, so it takes the full goSPL timestep —
        adequate for goSPL's km / 10²–10⁴ yr scale where explicit subcycling is
        prohibitive. Free boundary `H≥0` by post-clamp; validated against the
        analytical SIA dome (ice-volume conservation under zero mass balance).
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

    def iceAccumulation(self):
        """
        Main Ice Accumulation Calculation.

        This method evolves the ice thickness with the SIA solve
        (:meth:`_iceFlowSIA`) using the dynamic ice-cap altitude, the
        equilibrium-line altitude (ELA) and the glacier terminus. It also
        snapshots the ice load for the flexural isostasy when enabled.
        """

        ti = process_time()

        # Snapshot the current ice load for flexure (applyFlexure loads the
        # increment iceHL - iceFlex), taken BEFORE the thickness is updated.
        if self.flexOn:
            self.iceHL.copy(result=self.iceFlex)

        # Get dynamic properties such as ice cap altitude and equilibrium-line altitude
        iceH = self.iceH(self.tNow)  # Ice Cap Altitude
        elaH = self.elaH(self.tNow)  # Equilibrium-Line Altitude
        iceT = self.iceT(self.tNow)  # Glacier terminus

        # No ice when the whole surface is below the ELA, or for a degenerate
        # config where the ice-cap altitude is not above the ELA (the
        # (z - elaH) / (iceH - elaH) mass-balance ramp would be NaN/Inf).
        max_elev = self.hGlobal.max()[1]
        if max_elev < elaH or iceH <= elaH:
            self.iceHL.set(0.)
            self.iceUbL.set(0.)
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
        in the ablation zone, conserving rock volume exactly (the net bed change
        is zero by construction): a transport, not a source/sink.

        Conservation is structural: the total abraded volume ``Vtot`` is
        redistributed to the ablation cells weighted by the meltwater rate
        (``iceMeltL``), so ``Σ deposited == Vtot``. Guarded by a dedicated
        till-conservation test (the dual-lithology lesson: a volume the total
        budget can't see needs its own guard).

        Active only with ``till.on`` and ``Kg > 0``; a no-op otherwise.

        Simplification (documented): till deposits across the ablation zone
        weighted by melt rather than routed per ice-catchment to each terminus.
        When stratigraphy is on the till is layered into the stratigraphic
        record (and split into the coarse/fine lithology fractions under
        dual-lithology); otherwise it moves the bulk bed only.
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
