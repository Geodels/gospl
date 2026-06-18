import os
import gc
import sys
import petsc4py
import numpy as np
import pandas as pd
import pyshtools as pysh

from mpi4py import MPI
from scipy import spatial
from time import process_time

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import jacobiancoeff
    from gospl._fortran import getfacevelocity
    from gospl._fortran import advecupwind

MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()
MPIsize = petsc4py.PETSc.COMM_WORLD.Get_size()
MPIcomm = petsc4py.PETSc.COMM_WORLD


class GridProcess(object):
    """
    When running goSPL on a 2D grid (*i.e.* not a global simulation), this class defines two processes solved directly on the unstructured mesh:

    1. **Flexural isostasy**: it allows to compute isostatic deflections of Earth's lithosphere with uniform or non-uniform flexural rigidity. Evolving surface loads are defined from erosion/deposition values associated to modelled surface processes.
    2. **Orographic rain**: it accounts for change in rainfall patterns associated to change in topography. The orographic precipitation function is based on the `Smith & Barstad (2004) <https://journals.ametsoc.org/view/journals/atsc/61/12/1520-0469_2004_061_1377_altoop_2.0.co_2.xml>`_ linear model, solved **directly on the unstructured mesh in parallel** as two steady advection–relaxation equations (cloud water → hydrometeors → precipitation). The stratified airflow (mountain-wave) term of the original spectral solution is dropped, so the windward/lee rain-shadow is retained while no regular grid or FFT is required.

    For global simulation, the library `pyshtools <https://shtools.github.io/SHTOOLS/>`_ provides a framework to estimate global-scale flexural isostasy. 

    """

    def __init__(self):
        """
        Initialisation of the `gridProcess` class.
        """

        self.flexIDs = None
        self.flex = None

        if self.flexOn:
            self.rho_water = 1030.0
            self.localFlex = np.zeros(self.lpoints)
            self.flexTe_dh = None
            # Warm-start cache for the varying-Te iterative SH solve: the
            # previous step's converged DH-grid deflection seeds the next solve
            # (loads evolve slowly ⇒ far fewer Anderson iterations). Rank-0 only.
            self._flex_w_prev = None
            # Cached state for the parallel FV biharmonic flexure solver
            # (method == 'fem', flat models). Lm is geometry-only; the operator A
            # and its KSP are rebuilt only when the Te interval advances (never,
            # for constant Te) — see _cmptFlexFEM / _buildFlexFEM.
            self._flexLm = None
            self._flexA = None
            self._flexKSP = None
            self._flexTeNb = -2
            self._flexDirNodes = np.empty(0, dtype=int)
            # Step counter driving the flexure interval (`flex_interval`): the
            # load reference (hOldFlex) is snapshotted at interval starts and
            # flexure applied at interval ends, so the load accumulates in
            # between. -1 so the first step (incremented to 0) is an interval
            # start; interval=1 reproduces every-step behaviour exactly.
            self.flexCount = -1
        if self.oroOn:
            # Orographic precipitation is solved on the mesh (no regular grid).
            # The advection-relaxation operators depend only on the uniform,
            # constant wind and the fixed mesh, so they are assembled once and
            # cached; each step only the surface load (source) and the solves
            # change. Cloud-water / hydrometeor solutions are warm-started.
            self._oroKSP = None
            self._oroAc = None        # cloud-water operator  (v·∇ + 1/τc)
            self._oroAf = None        # hydrometeor operator  (v·∇ + 1/τf)
            self._oroAdvDiag = None   # diagonal of the upwind advection operator
            self._oroLcoeff = None    # off-diagonal advection coefficients
            self._oroQc = None        # cached cloud-water solution (warm start)
            self._oroQs = None        # cached hydrometeor solution (warm start)

        if self.flexOn:
            # Both flexure methods solve directly on the mesh ('fem' = parallel
            # FV biharmonic on the flat DMPlex, 'global' = spherical-harmonic).
            if MPIrank == 0 and self.flex_method == 'global':
                self._buildDHGrid()

        return

    def _updateTe(self):
        """
        Finds current elastic thickness map for the considered time interval.
        """

        nb = self.teNb
        if nb < len(self.tedata) - 1:
            if self.tedata.at[nb + 1, "start"] <= self.tNow:
                nb += 1

        if nb > self.teNb or nb == -1:
            if nb == -1:
                nb = 0
            self.teNb = nb
            if self.flex_method == 'global':
                loadData = np.load(self.tedata.at[nb, "tMap"])
                self.flexTe = loadData[self.tedata.at[nb, "tKey"]]
                del loadData
                if MPIrank == 0:
                    self.flexTe_dh = self._unstr2dh(self.flexTe)
            else:
                # Parallel FV biharmonic (flat, 'fem'): Te is needed per LOCAL
                # mesh node (no regular grid). Uniform value or per-vertex map
                # indexed by this rank's owned+ghost nodes.
                if self.tedata["tUni"][nb] == 0.:
                    loadData = np.load(self.tedata.at[nb, "tMap"])
                    teVal = loadData[self.tedata.at[nb, "tKey"]]
                    del loadData
                    self.flexTe_local = teVal[self.locIDs]
                else:
                    self.flexTe_local = self.tedata.at[nb, "tUni"] * np.ones(self.lpoints)

        return

    def _buildDHGrid(self):
        """
        Build the Driscoll-Healy (DH2, shape N x 2N) target grid for global
        flexure and precompute the interpolation weights linking the
        unstructured spherical mesh (``self.mCoords``) and the DH grid.

        Forward direction (mesh -> DH): 3-D inverse-distance weighting via
        ``scipy.spatial.cKDTree`` on the existing Cartesian coordinates,
        which avoids polar / dateline singularities.

        Backward direction (DH -> mesh): bilinear interpolation with
        precomputed integer indices and fractional offsets. The DH grid is
        implicitly extended at call time with a south-pole row and a
        wrap-around longitude column so every mesh point sits inside a
        2 x 2 stencil.

        Because the unstructured mesh is fixed for a given simulation,
        weights are built once and reused for every flexure step.
        """

        # ---- 1. DH2 grid coordinates --------------------------------------
        N = int(round(180.0 / self.flex_res_deg))
        if N % 2:
            N += 1
        self.dh_N = N
        self.dh_dlat = 180.0 / N
        self.dh_dlon = 360.0 / (2 * N)
        self.dh_lat = 90.0 - self.dh_dlat * np.arange(N)
        self.dh_lon = self.dh_dlon * np.arange(2 * N)

        # ---- 2. Forward weights: unstructured -> DH -----------------------
        lon2d, lat2d = np.meshgrid(np.deg2rad(self.dh_lon),
                                   np.deg2rad(self.dh_lat))
        cosLat = np.cos(lat2d)
        gridXYZ = np.stack([self.radius * cosLat * np.cos(lon2d),
                            self.radius * cosLat * np.sin(lon2d),
                            self.radius * np.sin(lat2d)],
                           axis=-1).reshape(-1, 3)

        tree = spatial.cKDTree(self.mCoords[:, :3], leafsize=10)
        distances, self.dhIDs = tree.query(gridXYZ, k=self.rgrd_interp)
        self.dhWeights = np.divide(
            1.0, distances ** 2,
            out=np.zeros_like(distances), where=distances != 0,
        )
        self.dhSumWeights = np.sum(self.dhWeights, axis=1)
        self.dhOnIDs = np.where(self.dhSumWeights == 0)[0]
        self.dhSumWeights[self.dhSumWeights == 0] = 1.0e-4

        # ---- 3. Backward indices: padded DH -> unstructured (bilinear) ----
        # Padded DH grid (built on demand inside `_dh2unstr`):
        #   lat: dh_lat with -90 appended  -> shape (N+1,), step  dh_dlat
        #   lon: dh_lon with 360 appended  -> shape (2N+1,), step dh_dlon
        norm = np.linalg.norm(self.mCoords[:, :3], axis=1)
        uLat = np.rad2deg(np.arcsin(np.clip(self.mCoords[:, 2] / norm,
                                            -1.0, 1.0)))
        uLon = np.rad2deg(np.arctan2(self.mCoords[:, 1],
                                     self.mCoords[:, 0])) % 360.0

        fracLat = (90.0 - uLat) / self.dh_dlat
        self.uIdxLat = fracLat.astype(np.intp)
        self.uFracLat = fracLat - self.uIdxLat
        mask = self.uIdxLat >= N
        self.uIdxLat[mask] = N - 1
        self.uFracLat[mask] = 1.0

        fracLon = uLon / self.dh_dlon
        self.uIdxLon = fracLon.astype(np.intp)
        self.uFracLon = fracLon - self.uIdxLon
        mask = self.uIdxLon >= 2 * N
        self.uIdxLon[mask] = 2 * N - 1
        self.uFracLon[mask] = 1.0

        # ---- 4. Constant elastic-operator eigenvalues ---------------------
        lmax = N // 2 - 1
        ll = np.arange(lmax + 1)
        self.dh_P_l = ((ll * (ll + 1)) ** 2 - 4.0 * ll * (ll + 1)) \
            / self.radius ** 4

        del tree, distances, gridXYZ, lon2d, lat2d, cosLat
        del norm, uLat, uLon, fracLat, fracLon, mask
        gc.collect()

        return

    def _unstr2dh(self, field):
        """
        Inverse-distance interpolation of a 1-D unstructured field defined at
        ``self.mCoords`` onto the DH2 grid (shape ``(dh_N, 2*dh_N)``).
        """

        g = np.sum(self.dhWeights * field[self.dhIDs], axis=1) \
            / self.dhSumWeights
        if len(self.dhOnIDs) > 0:
            g[self.dhOnIDs] = field[self.dhIDs[self.dhOnIDs, 0]]
        return g.reshape(self.dh_N, 2 * self.dh_N)

    def _dh2unstr(self, field_dh):
        """
        Bilinear interpolation of a DH2-grid array back to the unstructured
        mesh. The grid is padded at call time with a south-pole row (mean of
        the southernmost DH row) and a wrap-around longitude column.
        """

        N = self.dh_N
        south = field_dh[-1:, :].mean(axis=1, keepdims=True)
        south = np.broadcast_to(south, (1, 2 * N)).copy()
        padded = np.concatenate([field_dh, south], axis=0)
        padded = np.concatenate([padded, padded[:, :1]], axis=1)

        i, j = self.uIdxLat, self.uIdxLon
        fy, fx = self.uFracLat, self.uFracLon
        return (
            (1.0 - fy) * (1.0 - fx) * padded[i,     j    ] +
            (1.0 - fy) *        fx  * padded[i,     j + 1] +
                   fy  * (1.0 - fx) * padded[i + 1, j    ] +
                   fy  *        fx  * padded[i + 1, j + 1]
        )

    def _cmptFlexGlobal(self, erodep_dh, te_dh, rho_infill=0.0, max_iter=50,
                        flex_tol=5.0e-4, relax=1.0, anderson_depth=5):
        r"""
        Spectral thin-elastic-shell flexure solve on the precomputed
        Driscoll-Healy (DH2) grid.

        Parameters
        ----------
        erodep_dh : np.ndarray, shape ``(dh_N, 2*dh_N)``
            Erosion(-) / deposition(+) thickness in metres on the DH2 grid.
            Produced by :meth:`_unstr2dh` from the unstructured mesh field.
        te_dh : float or np.ndarray of shape ``(dh_N, 2*dh_N)``
            Elastic thickness in metres. Scalar -> constant-Te spectral
            solve. Array (same grid as ``erodep_dh``) -> iterative
            varying-Te solve.
        rho_infill : float
            Density (kg/m^3) of the material filling the flexural moat
            (0 = air, 1030 = sea water).
        max_iter, flex_tol, relax, anderson_depth : int, float, float, int
            Picard / Anderson iteration controls for the varying-Te branch.
            ``relax < 1`` damps the update (useful for strong Te contrasts).

        Returns
        -------
        np.ndarray, shape ``(dh_N, 2*dh_N)``
            Flexural deflection in metres on the DH2 grid. Sign:
            **negative = down** (subsidence under deposition; rebound gives
            positive w under erosion).
        """

        # ---- 1. Te field / reference rigidity ------------------------------
        if isinstance(te_dh, np.ndarray):
            # Maximum Te as reference guarantees dD = D(x) - D0 <= 0 and
            # ||dD/D0||_inf < 1, required for Picard contraction. Mean-Te
            # diverges when oceanic Te << continental Te.
            Te0 = float(np.max(te_dh))
            varying_te = True
        else:
            Te0 = float(te_dh)
            te_dh = None
            varying_te = False

        if Te0 <= 0:
            raise ValueError("elastic thickness must be positive")

        D0 = self.young * Te0 ** 3 / (12.0 * (1.0 - self.nu ** 2))
        if varying_te:
            D_field = self.young * te_dh ** 3 / (12.0 * (1.0 - self.nu ** 2))
            dD = D_field - D0
        drho = self.flex_rhoa - rho_infill

        # ---- 2. spectral filter (P_l precomputed in _buildDHGrid) ----------
        P_l = self.dh_P_l
        filt = 1.0 / (drho * self.gravity + D0 * P_l)
        filt[:2] = 0.0   # drop degree 0 (mean) and 1 (centre-of-mass drift)

        def _sh_solve(q_field):
            qlm = pysh.SHGrid.from_array(q_field, grid='DH').expand()
            qlm.coeffs[:] *= filt[None, :, None]
            return qlm.expand(grid='DH2', extend=False).to_array()

        def _elastic_op(w_field):
            wlm = pysh.SHGrid.from_array(w_field, grid='DH').expand()
            wlm.coeffs[:] *= P_l[None, :, None]
            return wlm.expand(grid='DH2', extend=False).to_array()

        # ---- 3. surface load -----------------------------------------------
        q_dh = erodep_dh * self.flex_rhos * self.gravity   # Pa, signed

        # ---- 4. solve ------------------------------------------------------
        if not varying_te:
            w_dh = _sh_solve(q_dh)
            if self.verbose:
                print(f"[shflex] constant Te = {Te0:.1f} m -> single-pass solve")
        else:
            # Warm-start: the previous step's converged deflection is usually
            # the closest available guess (loads evolve slowly between steps),
            # so it cuts Anderson iterations vs restarting from the constant-Te0
            # solution. The iteration is a contraction, so the converged result
            # is unchanged regardless of the initial guess. Fall back to the
            # constant-Te0 solution on the first step / if the grid changed.
            w_prev = getattr(self, "_flex_w_prev", None)
            if w_prev is not None and w_prev.shape == q_dh.shape:
                w_dh = w_prev
            else:
                w_dh = _sh_solve(q_dh)
            F_hist = []   # F(w_k)        - Anderson history
            g_hist = []   # F(w_k) - w_k  - residuals
            rel = np.nan
            for it in range(max_iter):
                Mw    = _elastic_op(w_dh)
                q_eff = q_dh - dD * Mw
                Fw    = _sh_solve(q_eff)
                g     = Fw - w_dh

                rel = np.linalg.norm(g) / max(np.linalg.norm(Fw), 1.0e-30)
                if self.verbose:
                    print(f"[shflex] iter {it:3d}  |dw|/|w| = {rel:.2e}")
                if rel < flex_tol:
                    w_dh = Fw
                    break

                F_hist.append(Fw)
                g_hist.append(g)
                if len(F_hist) > anderson_depth + 1:
                    F_hist.pop(0)
                    g_hist.pop(0)

                if anderson_depth > 0 and len(g_hist) >= 2:
                    # Anderson type-II: find gamma minimising ||g_n - dG gamma||,
                    # then w_{n+1} = F(w_n) - dF gamma.
                    dG = np.stack([g_hist[i + 1] - g_hist[i]
                                   for i in range(len(g_hist) - 1)])
                    dF = np.stack([F_hist[i + 1] - F_hist[i]
                                   for i in range(len(F_hist) - 1)])
                    dG_flat = dG.reshape(dG.shape[0], -1).T
                    dF_flat = dF.reshape(dF.shape[0], -1).T
                    gamma, *_ = np.linalg.lstsq(dG_flat,
                                                g_hist[-1].ravel(),
                                                rcond=None)
                    w_new = Fw - (dF_flat @ gamma).reshape(Fw.shape)
                else:
                    w_new = Fw

                if relax != 1.0:
                    w_new = relax * w_new + (1.0 - relax) * w_dh
                w_dh = w_new
            else:
                if self.verbose:
                    print(f"[shflex] warning: did not converge in {max_iter} "
                          f"iters (rel={rel:.2e})")

            # Cache the converged deflection to warm-start the next step.
            self._flex_w_prev = w_dh.copy()

        return -w_dh

    def _cmptFlexFEM(self, dzLocal):
        r"""
        **Parallel finite-volume biharmonic flexure for FLAT (2D) models**
        (``flexure: method: 'fem'`` — the flat-model flexure solver).

        Solves the plate equation with a single PETSc solve **directly on the
        DMPlex**: fully parallel, no gather-to-root, no regular grid, and
        **varying elastic thickness in one linear solve** (no iteration over the
        rigidity contrast). It replaced the former serial ``'FD'`` and ``'FFT'``
        solvers, which gathered the load to one rank, solved on a regular grid,
        and interpolated back.

        Solves the thin-plate-on-elastic-foundation equation with a spatially
        varying rigidity :math:`D(x)=E\,T_e(x)^3/[12(1-\nu^2)]`

        .. math::
          \nabla^2\!\big(D\,\nabla^2 w\big) + \Delta\rho\,g\,w = q,
          \qquad q = \rho_s\,g\,\mathrm{(erodep)}

        Discretised with the cached finite-volume negative-Laplacian ``Lm``
        (:math:`= -\nabla^2`, the hillslope stencil) applied twice, this is the
        **single-field** SPD-structured system

        .. math::
          \big[\,L_m\,\mathrm{diag}(D)\,L_m + \Delta\rho g\,I\,\big]\,w = q,

        solved with GMRES + GAMG (``Lm`` is row-area-normalised, hence slightly
        non-symmetric, so GMRES rather than CG). The elastic foundation
        :math:`\Delta\rho g\,I` (a positive diagonal) makes the operator
        well-posed — it pins the absolute deflection, so there is **no nullspace
        and no mean removal** (unlike the closed-sphere spectral solver). The
        solver defaults to a cached DIRECT factorisation — serial PETSc LU,
        parallel MUMPS — reused every step (only the RHS changes); it is
        options-prefixed (``flexfem_``) so an iterative Krylov method can be
        requested for meshes too large to factorise.

        When this solver replaced the former regular-grid solver it agreed with it
        on a flat mesh where the deflection decays inside the domain: correlation
        0.998 (natural BC), 0.9996 (clamped). Where the flexural wavelength
        approaches the domain size the two boundary discretisations differ by ~10 %.

        .. note::
            Boundary conditions (per side, from ``flex_bcN/S/E/W``):

            * ``0Slope0Shear`` and ``Mirror`` — the natural zero-flux FV boundary
              (``w'=0, w'''=0``; a thin plate's Mirror *is* 0Slope0Shear), imposed
              for free by the FV Laplacian — no modification.
            * ``0Displacement0Slope`` (clamped) — pin ``w=0`` (Dirichlet
              ``zeroRows`` on that side's nodes); the natural ``w'=0`` from the
              inner FV Laplacian completes the clamp.

            Sides map geographically. ``0Moment0Shear`` and ``Periodic`` are not
            implemented.

            The operator ``A`` and its KSP are CACHED and reused across steps:
            the geometry (``Lm``) is fixed and, for constant or piecewise-constant
            ``Te``, ``A`` only changes when the ``Te`` interval advances. Each step
            then re-uses the factorisation / preconditioner and only the RHS
            changes — so after the one-off setup a step costs a back-substitution
            (serial) or a few warm Krylov iterations (parallel). Serial runs use a
            direct LU (fast on small meshes); parallel runs use
            GMRES+GAMG. The KSP is options-prefixed (``flexfem_``) to override.

        :arg dzLocal: local (lpoints) elevation change = surface load thickness (m).

        :return: local (lpoints) flexural deflection (m), same sign convention as
            ``_cmptFlexGlobal`` (negative = subsidence under deposition).
        """

        # ---- rigidity D(x); track the Te interval so the cached operator is
        #      rebuilt only when Te actually changes (never, for constant Te). --
        if self.tedata is not None:
            self._updateTe()
            Te = self.flexTe_local
            teNb = self.teNb
        else:
            Te = self.flex_eet * np.ones(self.lpoints)
            teNb = -1
        rebuild = (self._flexKSP is None) or (teNb != self._flexTeNb)

        if rebuild:
            self._buildFlexFEM(Te)
            self._flexTeNb = teNb

        # ---- RHS q = ρ_s g·(erodep); pin w=0 at any clamped boundary nodes ----
        qloc = dzLocal * self.flex_rhos * self.gravity
        if len(self._flexDirNodes) > 0:
            qloc[self._flexDirNodes] = 0.0
        self.tmpL.setArray(qloc)
        self.dm.localToGlobal(self.tmpL, self.tmp)       # self.tmp = q (global)

        sol = self.hGlobal.duplicate()
        self._flexKSP.solve(self.tmp, sol)
        r = self._flexKSP.getConvergedReason()
        if r < 0 and MPIrank == 0:
            print(
                "[femflex] KSP failed to converge (reason %d) after %d its"
                % (r, self._flexKSP.getIterationNumber()),
                flush=True,
            )

        self.dm.globalToLocal(sol, self.tmpL)
        w = self.tmpL.getArray().copy()
        sol.destroy()

        return -w

    def _buildFlexFEM(self, Te):
        r"""
        Assemble (and cache) the flat-model FV biharmonic operator
        ``A = Lm·diag(D)·Lm + Δρg·I`` with its boundary conditions, and set up the
        reusable KSP. Called from :meth:`_cmptFlexFEM` on the first flexure step
        and whenever the elastic-thickness interval advances; the result is reused
        for every step in between (only the RHS changes).
        """

        D = self.young * Te ** 3 / (12.0 * (1.0 - self.nu ** 2))
        kfound = self.flex_rhoa * self.gravity          # Δρ·g (rho_infill = 0)

        # FV negative-Laplacian Lm (= -∇²); geometry only ⇒ cache.
        if self._flexLm is None:
            coeffs = jacobiancoeff(
                self.hLocal.getArray(),
                np.ones(self.lpoints),
                np.zeros(self.lpoints),
            )
            self._flexLm = self._assembleDiffMatCSR(coeffs)
        Lm = self._flexLm

        # A = Lm·diag(D)·Lm + Δρg·I (biharmonic + elastic foundation).
        Dvec = self.hGlobal.duplicate()
        self.tmpL.setArray(D)
        self.dm.localToGlobal(self.tmpL, Dvec)
        LmD = Lm.copy()
        LmD.diagonalScale(R=Dvec)                        # Lm·diag(D)
        A = LmD.matMult(Lm)                              # Lm·diag(D)·Lm
        A.shift(kfound)                                  # + Δρg·I
        LmD.destroy()
        Dvec.destroy()

        # Per-side BCs. 0Slope0Shear and Mirror are the natural zero-flux FV
        # boundary (w'=0, w'''=0; a thin plate's Mirror IS 0Slope0Shear) — no
        # modification. 0Displacement0Slope (clamped) pins w=0 (Dirichlet
        # zeroRows); the natural w'=0 from the inner FV Laplacian completes the
        # clamp. Sides map geographically (N=ymax, E=xmax, S=ymin, W=xmin). The clamped node
        # set is cached for the per-step RHS. `flex_bc*` are identical on every
        # rank, so the zeroRows guard is collective-safe.
        sides = (
            (self.southPts, self.flex_bcS),
            (self.eastPts, self.flex_bcE),
            (self.northPts, self.flex_bcN),
            (self.westPts, self.flex_bcW),
        )
        if any(bc == "0Displacement0Slope" for _, bc in sides):
            dbc = [
                pts for pts, bc in sides
                if bc == "0Displacement0Slope" and pts is not None and len(pts) > 0
            ]
            dnodes = (
                np.unique(np.concatenate(dbc)) if dbc else np.empty(0, dtype=int)
            )
            owned = dnodes[self.inIDs[dnodes] == 1]
            grows = (
                Lm.getLGMap()[0].apply(owned.astype(petsc4py.PETSc.IntType))
                if len(owned) else np.empty(0, dtype=petsc4py.PETSc.IntType)
            )
            A.zeroRows(grows, diag=1.0)                  # collective; clamps w=0
            self._flexDirNodes = dnodes
        else:
            self._flexDirNodes = np.empty(0, dtype=int)

        # Reusable KSP — a cached DIRECT factorisation, reused every step (only
        # the RHS changes), so a step costs a back-substitution. Serial uses
        # PETSc's built-in LU; parallel uses MUMPS (the GMRES+GAMG alternative
        # does not converge on the stiff biharmonic at realistic Te, hitting the
        # iteration cap — a direct solve is both faster and robust here, and the
        # one-off factorisation amortises over the run). Options-prefixed
        # (`flexfem_`) so an iterative solver can be requested for meshes too
        # large to factorise.
        if self._flexKSP is not None:
            self._flexKSP.destroy()
        if self._flexA is not None:
            self._flexA.destroy()
        self._flexA = A
        ksp = petsc4py.PETSc.KSP().create(petsc4py.PETSc.COMM_WORLD)
        ksp.setOptionsPrefix("flexfem_")
        ksp.setOperators(A)
        ksp.setType("preonly")
        pc = ksp.getPC()
        pc.setType("lu")
        if MPIsize > 1:
            pc.setFactorSolverType("mumps")              # parallel direct solve
        ksp.setFromOptions()
        ksp.setUp()
        pc.setReusePreconditioner(True)                  # keep the factor across steps
        self._flexKSP = ksp

        return

    def applyFlexure(self):
        r"""
        This function computes the flexural isostasy equilibrium based on topographic change. It is a simple routine that accounts for flexural isostatic rebound associated with erosional loading/unloading.

        The function takes an initial (at time t) and final topography (at time t+Dt) (i.e. before and after erosion/deposition) and returns a corrected final topography that includes the effect of erosional/depositional unloading/loading. It uses a spectral method to solve the bi-harmonic equation governing the bending/flexure of a thin elastic plate floating on an inviscid fluid (the asthenosphere).

        .. math::

          D (d^4 w / d^4 x ) + \Delta \rho g w = q

        where :math:`D` is the flexural rigidity,  :math:`w` is vertical deflection of the plate, :math:`q` is the applied surface load, and :math:`\Delta \rho = \rho_m − \rho_f` is the density of the mantle minus the density of the infilling material.
        """

        t0 = process_time()

        # Get elevations from time of equilibrium and after erosion deposition.
        # Build the load locally (lpoints); the flexure solve is serial on rank
        # 0, so the global field is assembled there from the owned nodes only.
        hl = self.hLocal.getArray().copy()
        local_dZ = hl - self.hOldFlex.getArray()

        # If glaciers exist then add corresponding equivalent sediment thickness
        if self.iceOn:
            dIce = self.iceHL.getArray() - self.iceFlex.getArray()
            local_dZ += dIce * 910.0 / self.flex_rhos  # 910 kg/m3 ice density

        # Parallel mixed-FV biharmonic flexure for FLAT models (opt-in
        # `method: 'fem'`): solved directly on the DMPlex with no gather-to-root,
        # no regular grid. Returns the local deflection; apply and exit.
        if self.flex_method == 'fem':
            tmpFlex = self._cmptFlexFEM(local_dZ)
            self.localFlex += tmpFlex
            self.hLocal.setArray(hl + tmpFlex)
            self.dm.localToGlobal(self.hLocal, self.hGlobal)
            if MPIrank == 0 and self.verbose:
                print(
                    "Compute Flexural Isostasy [fem] (%0.02f seconds)"
                    % (process_time() - t0),
                    flush=True,
                )
            return

        # Global (spherical-harmonic) flexure: gather the load to rank 0, solve
        # on the Driscoll-Healy grid, scatter back. (The flat 'fem' method has
        # already returned above; the former 'FD'/'FFT' solvers were removed.)
        dZ = self._gatherGlobalOnRoot(local_dZ)
        flexZ = None
        if self.tedata is not None:
            self._updateTe()
        if MPIrank == 0:
            erodep_dh = self._unstr2dh(dZ)
            te_dh = self.flexTe_dh if self.tedata is not None \
                else float(self.flex_eet)
            wflex_dh = self._cmptFlexGlobal(
                erodep_dh, te_dh,
                max_iter=self.flex_max_iter,
                flex_tol=self.flex_tol,
                relax=self.flex_relax,
            )
            flexZ = self._dh2unstr(wflex_dh)

        # Send flexural response globally
        flexZ = MPI.COMM_WORLD.bcast(flexZ, root=0)
        tmpFlex = flexZ[self.locIDs]

        self.localFlex += tmpFlex

        # Update elevation
        self.hLocal.setArray(hl + tmpFlex)
        self.dm.localToGlobal(self.hLocal, self.hGlobal)

        if MPIrank == 0 and self.verbose:
            print(
                "Compute Flexural Isostasy (%0.02f seconds)" % (process_time() - t0)
            )

        return

    def _buildOroMat(self, diag, offcoeff, dirichlet=False):
        """
        Assemble a PETSc operator on the DMPlex from a per-row diagonal and the
        upwind off-diagonal coefficients returned by ``advecupwind`` (one column
        per neighbour, ``FVmesh_ngbID`` giving the neighbour indices). This is the
        same incremental local-CSR assembly used for the tectonic advection
        matrix.

        :arg diag: diagonal entries (length ``lpoints``)
        :arg offcoeff: off-diagonal coefficients, shape ``(lpoints, maxnb)``
        :arg dirichlet: when True, the (non-cyclic) domain edge rows are replaced
            by an identity row so the precipitation tracers are pinned to zero
            there (clean inflow / zero-padding equivalent).

        :return: assembled PETSc matrix
        """

        d = diag.copy()
        off = offcoeff.copy()
        if dirichlet:
            d[self.advectBorders] = 1.0
            off[self.advectBorders, :] = 0.0

        mat = self._matrix_build_diag(d)
        indptr = np.arange(0, self.lpoints + 1, dtype=petsc4py.PETSc.IntType)
        for k in range(0, self.maxnb):
            tmpMat = self._matrix_build()
            indices = self.FVmesh_ngbID[:, k].copy()
            data = off[:, k]
            ids = np.nonzero(data == 0.0)
            indices[ids] = ids
            tmpMat.assemblyBegin()
            tmpMat.setValuesLocalCSR(
                indptr,
                indices.astype(petsc4py.PETSc.IntType),
                data,
            )
            tmpMat.assemblyEnd()
            mat.axpy(1.0, tmpMat)
            tmpMat.destroy()

        return mat

    def _oroSolve(self, matrix, rhs, sol):
        """
        Solve one orographic advection-relaxation system with a dedicated
        (cached) GMRES + block-Jacobi KSP. The operators are non-symmetric but
        diagonally-dominant M-matrices, so this converges quickly; the solution
        is warm-started from the previous step.
        """

        if self._oroKSP is None:
            ksp = petsc4py.PETSc.KSP().create(petsc4py.PETSc.COMM_WORLD)
            ksp.setType("gmres")
            ksp.getPC().setType("bjacobi")
            ksp.setTolerances(rtol=1.0e-8)
            ksp.setInitialGuessNonzero(True)
            ksp.setOptionsPrefix("oro_")
            petsc4py.PETSc.Options()["oro_sub_pc_factor_shift_type"] = "nonzero"
            ksp.setFromOptions()
            self._oroKSP = ksp

        ksp = self._oroKSP
        ksp.setOperators(matrix, matrix)
        ksp.solve(rhs, sol)
        petsc4py.PETSc.garbage_cleanup()

        return sol

    def cptOrography(self):
        """
        Linear Theory of Orographic Precipitation following `Smith & Barstad (2004) <https://journals.ametsoc.org/view/journals/atsc/61/12/1520-0469_2004_061_1377_altoop_2.0.co_2.xml>`_, solved **directly on the unstructured mesh in parallel**.

        The model is cast as two steady, vertically-integrated advection-relaxation
        equations for the cloud-water density :math:`q_c` and the hydrometeor
        density :math:`q_s`, advected by a uniform wind :math:`\\mathbf{v}`:

        .. math::

            (\\mathbf{v}\\cdot\\nabla + 1/\\tau_c)\\,q_c = C_w\\,\\mathbf{v}\\cdot\\nabla h, \\qquad
            (\\mathbf{v}\\cdot\\nabla + 1/\\tau_f)\\,q_s = q_c/\\tau_c,

        with the precipitation rate :math:`P = q_s/\\tau_f`. The source is the
        terrain-forced uplift :math:`C_w\\,\\mathbf{v}\\cdot\\nabla h`: positive
        (condensation) on windward slopes, negative (evaporation) on lee slopes,
        which together with the downwind advection of :math:`q_c, q_s` produces
        the windward-wet / lee-dry rain shadow. The elevation in the source is
        clamped to sea level, so submarine bathymetry produces no orographic
        forcing (the airflow follows the flat sea surface, not the seafloor).

        In Fourier space these equations recover the Smith-Barstad transfer
        function with the stratified mountain-wave term :math:`(1-i\\,h_w m)` set
        to one. That airflow term is therefore dropped (and the parameters
        ``nm``, ``hw`` and ``latitude`` are inert); everything else is identical.
        The advection operator is the first-order upwind finite-volume scheme on
        the Voronoi mesh, so the operators depend only on the (uniform, constant)
        wind and are assembled once and cached.

        .. note::

            Please refer to the original manuscript of Smith and Barstad (2004) to understand the model physics and limitations.

        Common bounds:

         - precip_base : 0-10 [mm/h]
         - precip_min : 0.001 - 1  [mm/h]
         - conv_time :  200-2000 [s]
         - fall_time :  200-2000 [s]
         - cw :  0.001-0.02 [kg/m^3]
         - rainfall_frequency : 0.1 - 24 [number of storms of 1 hour duration per day]
        """

        t0 = process_time()

        # Ensure the local elevation (with halo) is current for the FV operator.
        self.dm.globalToLocal(self.hGlobal, self.hLocal)

        # Assemble and cache the advection-relaxation operators. They depend only
        # on the uniform wind and the fixed mesh, so this is done once.
        if self._oroAc is None:
            # Uniform wind in the flat (x, y) plane.
            u0 = -np.sin(self.wind_dir * 2 * np.pi / 360) * self.wind_speed
            v0 = np.cos(self.wind_dir * 2 * np.pi / 360) * self.wind_speed
            nodeVel = np.zeros((self.lpoints, 3))
            nodeVel[:, 0] = u0
            nodeVel[:, 1] = v0
            getfacevelocity(self.lpoints, nodeVel)

            # advecupwind(dt=1) returns (I + L) where L is the upwind FV operator
            # for v·∇. So L has diagonal lcoeff[:,0]-1 and off-diagonals lcoeff[:,1:].
            lcoeff = advecupwind(self.lpoints, 1.0)
            self._oroAdvDiag = lcoeff[:, 0] - 1.0
            self._oroLcoeff = lcoeff[:, 1:].copy()
            # Cloud-water (v·∇ + 1/τc) and hydrometeor (v·∇ + 1/τf) operators with
            # zero-Dirichlet domain edges.
            self._oroAc = self._buildOroMat(
                self._oroAdvDiag + 1.0 / self.oro_conv_time,
                self._oroLcoeff, dirichlet=True,
            )
            self._oroAf = self._buildOroMat(
                self._oroAdvDiag + 1.0 / self.oro_fall_time,
                self._oroLcoeff, dirichlet=True,
            )
            self._oroQc = self.hGlobal.duplicate()
            self._oroQc.set(0.0)
            self._oroQs = self.hGlobal.duplicate()
            self._oroQs.set(0.0)

        # Source S = Cw (v·∇h) via the same upwind operator (local matvec using
        # the haloed elevation), zeroed on the Dirichlet edges. The uplift is
        # forced by the surface the airflow follows, which over the ocean is the
        # (flat) sea surface, not the seafloor — so the elevation is clamped to
        # sea level: there is no orographic forcing over submarine bathymetry,
        # only where land rises above the sea.
        hL = np.maximum(self.hLocal.getArray(), self.sealevel)
        src = self._oroAdvDiag * hL
        for k in range(0, self.maxnb):
            src = src + self._oroLcoeff[:, k] * hL[self.FVmesh_ngbID[:, k]]
        src *= self.oro_cw
        src[self.advectBorders] = 0.0
        self.tmpL.setArray(src)
        self.dm.localToGlobal(self.tmpL, self.tmp)

        # Cloud water: Ac q_c = S
        self._oroSolve(self._oroAc, self.tmp, self._oroQc)

        # Hydrometeors: Af q_s = q_c / τc
        self._oroQc.copy(result=self.tmp1)
        self.tmp1.scale(1.0 / self.oro_conv_time)
        self.dm.globalToLocal(self.tmp1, self.tmpL)
        arr = self.tmpL.getArray()
        arr[self.advectBorders] = 0.0
        self.tmpL.setArray(arr)
        self.dm.localToGlobal(self.tmpL, self.tmp1)
        self._oroSolve(self._oroAf, self.tmp1, self._oroQs)

        # Precipitation P = q_s/τf; convert to the goSPL rainfall units (mirrors
        # the previous spectral implementation's post-processing exactly).
        self.dm.globalToLocal(self._oroQs, self.tmpL)
        P = self.tmpL.getArray().copy() / self.oro_fall_time
        P *= 3600.  # mm hr-1
        P += self.oro_precip_base
        # Precipitation rate must exceed the minimum to avoid zero/negative runoff.
        P[P <= self.oro_precip_min] = self.oro_precip_min
        # Conversion from mm/hr to m/yr.
        P *= 0.366 * self.rainfall_frequency

        self.rainVal = P
        self.bL.setArray(self.rainVal * self.larea)
        self.dm.localToGlobal(self.bL, self.bG)

        if MPIrank == 0 and self.verbose:
            print(
                "Compute Orographic Rain (%0.02f seconds)" % (process_time() - t0)
            )

        return
