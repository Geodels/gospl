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
    from gflex.f2d import F2D
    from gospl._fortran import flexure

MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()
MPIsize = petsc4py.PETSc.COMM_WORLD.Get_size()


class GridProcess(object):
    """
    When running goSPL on a 2D grid (*i.e.* not a global simulation), this class defines two processes operating on a regular grid:

    1. **Flexural isostasy**: it allows to compute isostatic deflections of Earth's lithosphere with uniform or non-uniform flexural rigidity. Evolving surface loads are defined from erosion/deposition values associated to modelled surface processes.
    2. **Orographic rain**: it accounts for change in rainfall patterns associated to change in topography. The orographic precipitation function is based on `Smith & Barstad (2004) <https://journals.ametsoc.org/view/journals/atsc/61/12/1520-0469_2004_061_1377_altoop_2.0.co_2.xml>`_ linear model.

    For global simulation, the library `pyshtools <https://shtools.github.io/SHTOOLS/>`_ provides a framework to estimate global-scale flexural isostasy. 

    """

    def __init__(self):
        """
        Initialisation of the `gridProcess` class.
        """

        self.flexIDs = None
        self.flex = None
        self.xIndices = None

        if self.flexOn:
            self.rho_water = 1030.0
            self.localFlex = np.zeros(self.lpoints)
            if self.flex_method == 'FFT':
                self.boundflex = self.boundCond.replace('0', '2')
                self.boundflex = self.boundflex.replace('1', '0')
                self.boundflex = self.boundflex.replace('2', '1')
            self.flexTe_dh = None
        if self.oroOn:
            self.oroEPS = np.finfo(float).eps

        if self.flexOn or self.oroOn:
            # Build regular grid for flexure or orographic precipitation calculation
            if MPIrank == 0 and self.flex_method != 'global':
                self._buildRegGrid()
            if MPIrank == 0 and self.flexOn and self.flex_method == 'global':
                self._buildDHGrid()

        return

    def _regInterp(self, field):
        """
        Perform bilinear interpolation of ``field`` on the regular grid to unstructured 2D mesh.

        :arg field: data to interpolate of size m x n

        :return: ufield ``field`` interpolated to unstructured nodes
        """

        ufield = \
            (1. - self.xFrac) * (1. - self.yFrac) * field[self.yIndices, self.xIndices] + \
            self.xFrac * (1. - self.yFrac) * field[self.yIndices, self.xIndices + 1] + \
            (1. - self.xFrac) * self.yFrac * field[self.yIndices + 1, self.xIndices] + \
            self.xFrac * self.yFrac * field[self.yIndices + 1, self.xIndices + 1]

        return ufield

    def _buildRegGrid(self):
        """
        Builds the regular grid based on nodes coordinates and instantiates two interpolation objects.

        1. The first one uses  `SciPy cKDTree <https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.cKDTree.html>`_ to interpolate values from the unstructured mesh onto the regular grid based on an inverse weighting distance approach.

        2. The second performs a bilinear interpolation from the regular grid to unstructured 2D mesh.

        .. note::
            Here the KDTree is not kept in memory, instead we store the interpolation information, namely the indices of the neighbouring nodes, and the weights of each node in the neighborhood (based on the distance). This implies that the initial distribution of the mesh coordinates remains fixed over the simulation time window.
        """

        # Build regular grid for flexure and orographic precipitation calculation
        xmin = self.mCoords[:, 0].min()
        xmax = self.mCoords[:, 0].max()
        ymin = self.mCoords[:, 1].min()
        ymax = self.mCoords[:, 1].max()

        newx = np.arange(xmin, xmax + self.reg_dx, self.reg_dx)
        newy = np.arange(ymin, ymax + self.reg_dx, self.reg_dx)
        rx, ry = np.meshgrid(newx, newy)
        rPts = np.stack((rx.ravel(), ry.ravel())).T
        xmin, xmax = newx[0], newx[-1]
        ymin, ymax = newy[0], newy[-1]
        self.reg_xl = xmax - xmin
        self.reg_yl = ymax - ymin

        self.reg_nx = int(self.reg_xl / self.reg_dx + 1)
        self.reg_ny = int(self.reg_yl / self.reg_dx + 1)

        assert np.all(self.mCoords[:, 0] >= newx[0])
        assert np.all(self.mCoords[:, 0] <= newx[-1])
        assert np.all(self.mCoords[:, 1] >= newy[0])
        assert np.all(self.mCoords[:, 1] <= newy[-1])

        self.xFrac = np.interp(self.mCoords[:, 0], newx, np.arange(self.reg_nx))
        self.yFrac = np.interp(self.mCoords[:, 1], newy, np.arange(self.reg_ny))

        self.xIndices = np.array(self.xFrac, dtype=int)
        self.xFrac -= self.xIndices
        self.yIndices = np.array(self.yFrac, dtype=int)
        self.yFrac -= self.yIndices

        mask = self.xIndices == self.reg_nx - 1
        self.xIndices[mask] -= 1
        self.xFrac[mask] += 1.
        mask = self.yIndices == self.reg_ny - 1
        self.yIndices[mask] -= 1
        self.yFrac[mask] += 1.

        treeT = spatial.cKDTree(self.mCoords[:, :2], leafsize=10)
        distances, self.regIDs = treeT.query(rPts, k=self.rgrd_interp)
        # Inverse weighting distance...
        self.regWeights = np.divide(
            1.0, distances ** 2, out=np.zeros_like(distances), where=distances != 0
        )
        self.regSumWeights = np.sum(self.regWeights, axis=1)
        self.regOnIDs = np.where(self.regSumWeights == 0)[0]
        self.regSumWeights[self.regSumWeights == 0] = 1.0e-4

        del treeT, distances, mask, rPts, newx, newy
        gc.collect()

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
            elif self.tedata["tUni"][nb] == 0.:
                loadData = np.load(self.tedata.at[nb, "tMap"])
                teVal = loadData[self.tedata.at[nb, "tKey"]]
                del loadData
                self.flexTe = np.sum(self.regWeights * teVal[self.regIDs][:, :], axis=1) / self.regSumWeights
                if len(self.regOnIDs) > 0:
                    self.flexTe[self.regOnIDs] = teVal[self.regIDs[self.regOnIDs, 0]]
                self.flexTe = self.flexTe.reshape(self.reg_ny, self.reg_nx)
            else:
                self.flexTe = self.tedata.at[nb, "tUni"] * np.ones((self.reg_ny, self.reg_nx))

        return

    def _cptFlex2D(self, dZ):
        """
        Compute the flexural response for 2D cases.
        """

        # Build regular grid for flexure calculation
        if self.xIndices is None:
            self._buildRegGrid()

        # Interpolate values on the flexural regular grid
        regDiff = np.sum(self.regWeights * dZ[self.regIDs][:, :], axis=1) / self.regSumWeights
        if len(self.regOnIDs) > 0:
            regDiff[self.regOnIDs] = dZ[self.regIDs[self.regOnIDs, 0]]
        regDiff = regDiff.reshape(self.reg_ny, self.reg_nx)

        if self.flex_method == 'FFT':
            nFlex = flexure(regDiff, self.reg_ny, self.reg_nx, self.reg_yl, self.reg_xl,
                            self.young, self.nu, self.flex_rhos, self.flex_rhoa,
                            self.flex_eet, self.gravity, int(self.boundflex))

            # Interpolate back to goSPL mesh
            flexZ = self._regInterp(nFlex)

        elif self.flex_method == 'FD':
            simflex = F2D()
            simflex.Quiet = True

            simflex.Method = "FD"
            simflex.PlateSolutionType = "vWC1994"
            simflex.Solver = "direct"

            # gFlex parameters
            simflex.g = self.gravity
            simflex.E = self.young
            simflex.nu = self.nu
            simflex.rho_m = self.flex_rhoa
            simflex.rho_fill = 0.
            simflex.dx = self.reg_dx
            simflex.dy = self.reg_dx

            # Boundary conditions
            simflex.BC_E = self.flex_bcE
            simflex.BC_W = self.flex_bcW
            simflex.BC_S = self.flex_bcN
            simflex.BC_N = self.flex_bcS

            # Assign elastic thickness grid
            if self.tedata is not None:
                self._updateTe()
                simflex.Te = self.flexTe.copy()
            else:
                simflex.Te = self.flex_eet * np.ones(regDiff.shape)

            # Compute loads
            simflex.qs = self.flex_rhos * self.gravity * regDiff

            # Run gFlex
            simflex.initialize()
            simflex.run()
            simflex.finalize()

            # Interpolate back to goSPL mesh
            flexZ = self._regInterp(simflex.w)

            del simflex
            gc.collect()

        return flexZ

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
            # Constant-Te0 solution is much closer to the fixed point than
            # zero, so convergence starts immediately.
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

        return -w_dh

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

        dZ = self._gatherGlobalOnRoot(local_dZ)
        flexZ = None

        if self.flex_method == 'global':
            if self.tedata is not None:
                self._updateTe()

            if MPIrank == 0:
                erodep_dh = self._unstr2dh(dZ)
                te_dh = self.flexTe_dh if self.tedata is not None \
                    else float(self.flex_eet)
                wflex_dh = self._cmptFlexGlobal(erodep_dh, te_dh)
                flexZ = self._dh2unstr(wflex_dh)

        if MPIrank == 0 and self.flex_method != 'global':
            flexZ = self._cptFlex2D(dZ)

        # Send flexural response globally
        flexZ = MPI.COMM_WORLD.bcast(flexZ, root=0)
        tmpFlex = flexZ[self.locIDs]
        if self.flex_method != 'global':
            # Local flexural isostasy
            if self.south == 1:
                tmpFlex[self.southPts] = 0.
            if self.east == 1:
                tmpFlex[self.eastPts] = 0.
            if self.north == 1:
                tmpFlex[self.northPts] = 0.
            if self.west == 1:
                tmpFlex[self.westPts] = 0.

        self.localFlex += tmpFlex

        # Update elevation
        self.hLocal.setArray(hl + tmpFlex)
        self.dm.localToGlobal(self.hLocal, self.hGlobal)

        if MPIrank == 0 and self.verbose:
            print(
                "Compute Flexural Isostasy (%0.02f seconds)" % (process_time() - t0)
            )

        return

    def cptOrography(self):
        """
        Linear Theory of Orographic Precipitation following ` <https://journals.ametsoc.org/view/journals/atsc/61/12/1520-0469_2004_061_1377_altoop_2.0.co_2.xml>`_.

        The model includes airflow dynamics, condensed water advection, and downslope evaporation. It consists of two vertically-integrated steady-state advection equations describing: (i) the cloud water density and (ii) the hydrometeor density. Solving these equations using Fourier transform techniques, derives a single formula relating terrain and precipitation.

        .. note::

            Please refer to the original manuscript of Smith and Barstad (2004) to understand the model physics and limitations.

        Common bounds:

         - latitude : 0-90 [degrees]
         - precip_base : 0-10 [mm/h]
         - precip_min : 0.001 - 1  [mm/h]
         - conv_time :  200-2000 [s]
         - fall_time :  200-2000 [s]
         - nm :  0-0.1 [1/s]
         - hw  : 1000-5000 [m]
         - cw :  0.001-0.02 [kg/m^3]
         - rainfall_frequency : 0.1 - 24 [number of storms of 1 hour duration per day]
        """

        t0 = process_time()

        # Get elevations from the unstructured mesh structure. The orographic
        # precipitation solve is serial on rank 0, so assemble the global
        # elevation there from the owned nodes only (no mpoints array /
        # Allreduce on the other ranks).
        hl = self.hLocal.getArray().copy()
        newZ = self._gatherGlobalOnRoot(hl)

        if MPIrank == 0:
            # Build regular grid for flexure calculation
            if self.xIndices is None:
                self._buildRegGrid()

            # Interpolate values on the regular grid
            regNewZ = np.sum(self.regWeights * newZ[self.regIDs][:, :], axis=1) / self.regSumWeights
            if len(self.regOnIDs) > 0:
                regNewZ[self.regOnIDs] = newZ[self.regIDs[self.regOnIDs, 0]]
            regNewZ = regNewZ.reshape(self.reg_ny, self.reg_nx)

            # Wind components
            u0 = -np.sin(self.wind_dir * 2 * np.pi / 360) * self.wind_speed
            v0 = np.cos(self.wind_dir * 2 * np.pi / 360) * self.wind_speed
            # Coriolis factors
            f_coriolis = 2 * 7.2921e-5 * np.sin(self.wind_latitude * np.pi / 180)

            # Pad raster boundaries prior to FFT
            calc_pad = int(np.ceil(((sum(regNewZ.shape))) / 2) / 100 * 100)
            pad = min([calc_pad, 200])
            h = np.pad(regNewZ, pad, 'constant')
            nx, ny = h.shape

            # FFT
            hhat = np.fft.fft2(h)
            x_n_value = np.fft.fftfreq(ny, (1. / ny))
            y_n_value = np.fft.fftfreq(nx, (1. / nx))
            x_len = nx * self.reg_dx
            y_len = ny * self.reg_dx
            kx_line = 2 * np.pi * x_n_value / x_len
            ky_line = 2 * np.pi * y_n_value / y_len
            kx = np.tile(kx_line, (nx, 1))
            ky = np.tile(ky_line[:, None], (1, ny))

            # Vertical wave number (m)
            sigma = kx * u0 + ky * v0
            mf_num = self.oro_nm ** 2 - sigma ** 2
            mf_den = sigma ** 2 - f_coriolis ** 2

            # Numerical stability
            mf_num[mf_num < 0] = 0.
            mf_den[(mf_den < self.oroEPS) & (mf_den >= 0)] = self.oroEPS
            mf_den[(mf_den > -self.oroEPS) & (mf_den < 0)] = -self.oroEPS
            sign = np.where(sigma >= 0, 1, -1)
            m = sign * np.sqrt(np.abs(mf_num / mf_den * (kx ** 2 + ky ** 2)))

            # Transfer function
            P_karot = ((self.oro_cw * 1j * sigma * hhat) / 
                       ((1 - (self.oro_hw * m * 1j)) * 
                        (1 + (sigma * self.oro_conv_time * 1j)) * 
                        (1 + (sigma * self.oro_fall_time * 1j))))

            # Inverse FFT, de-pad, convert units, add uniform rate
            oroRain = np.fft.ifft2(P_karot)
            if pad > 0:
                oroRain = np.real(oroRain[pad:-pad, pad:-pad])
            else:
                oroRain = np.real(oroRain)
            oroRain *= 3600.  # mm hr-1
            oroRain += self.oro_precip_base
            # Precipitation rate must be a value greater than minimum precipitation/runoff to avoid errors when precip_rate <= 0
            oroRain[oroRain <= self.oro_precip_min] = self.oro_precip_min
            # Conversion from mm/hr to m/yr
            oroRain *= 0.366 * self.rainfall_frequency

            # Interpolate back to goSPL mesh
            oRain = self._regInterp(oroRain)
        else:
            oRain = None

        # Send orographic rain globally
        oroR = MPI.COMM_WORLD.bcast(oRain, root=0)

        # Local orographic rain values
        self.rainVal = oroR[self.locIDs]

        self.bL.setArray(self.rainVal * self.larea)
        self.dm.localToGlobal(self.bL, self.bG)

        if MPIrank == 0 and self.verbose:
            print(
                "Compute Orographic Rain (%0.02f seconds)" % (process_time() - t0)
            )

        return
