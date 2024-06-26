import os
import gc
import sys
import petsc4py
import numpy as np

from mpi4py import MPI
from scipy import spatial
from time import process_time

petsc4py.init(sys.argv)
MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()
MPIsize = petsc4py.PETSc.COMM_WORLD.Get_size()
MPIcomm = MPI.COMM_WORLD

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import flexure

petsc4py.init(sys.argv)
MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()
MPIsize = petsc4py.PETSc.COMM_WORLD.Get_size()
MPIcomm = MPI.COMM_WORLD


class GridProcess(object):
    """
    When running goSPL on a 2D grid (*i.e.* not a global simulation), this class defines two processes operating on a regular grid:

    1. **Flexural isostasy**: it allows to compute isostatic deflections of Earth's lithosphere with uniform or non-uniform flexural rigidity. Evolving surface loads are defined from erosion/deposition values associated to modelled surface processes.
    2. **Orographic rain**: it accounts for change in rainfall patterns associated to change in topography. The orographic precipitation function is based on `Smith & Barstad (2004) <https://journals.ametsoc.org/view/journals/atsc/61/12/1520-0469_2004_061_1377_altoop_2.0.co_2.xml>`_ linear model.
    """

    def __init__(self):
        """
        Initialisation of the `gridProcess` class.
        """

        self.flexIDs = None
        self.flex = None
        self.xIndices = None

        if self.flexOn:
            self.localFlex = np.zeros(self.lpoints)
            self.boundflex = self.boundCond.replace('0', '2')
            self.boundflex = self.boundflex.replace('1', '0')
            self.boundflex = self.boundflex.replace('2', '1')
        if self.oroOn:
            self.oroEPS = np.finfo(float).eps

        if self.flexOn or self.oroOn:
            # Build regular grid for flexure or orographic precipitation calculation
            if MPIrank == 0:
                self._buildRegGrid()

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

        The first one uses  `SciPy cKDTree <https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.cKDTree.html>`_ to interpolate values from the unstructured mesh onto the regular grid based on an inverse weighting distance approach.

        The second performs a bilinear interpolation from the regular grid to unstructured 2D mesh.

        .. note::
            Here that the KDTree is not kept in memory, instead we store the interpolation information, namely the indices of the neighbouring nodes, and the weights of each node in the neighborhood (based on the distance).
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
        distances, self.regIDs = treeT.query(rPts, k=4)
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

    def applyFlexure(self):
        r"""
        This function computes the flexural isostasy equilibrium based on topographic change. It is a simple routine that accounts for flexural isostatic rebound associated with erosional loading/unloading.

        The function takes an initial (at time t) and final topography (at time t+Dt) (i.e. before and after erosion/deposition) and returns a corrected final topography that includes the effect of erosional/depositional unloading/loading. It uses a spectral method to solve the bi-harmonic equation governing the bending/flexure of a thin elastic plate floating on an inviscid fluid (the asthenosphere).

        .. math::

          D (d^4 w / d^4 x ) + \Delta \rho g w = q

        where :math:`D` is the flexural rigidity,  :math:`w` is vertical deflection of the plate, :math:`q` is the applied surface load, and :math:`\Delta \rho = \rho_m − \rho_f` is the density of the mantle minus the density of the infilling material.

        .. warning ::
            This function assumes a value of 1011 Pa for Young's modulus, 0.25 for Poisson's ratio and 9.81 m/s2 for g, the gravitational acceleration.
        """

        t0 = process_time()

        # Get elevations from time of equilibrium and after erosion deposition
        hlflex = self.hOldFlex.getArray().copy()
        oldZ = np.zeros(self.mpoints, dtype=np.float64) - 1.0e8
        oldZ[self.locIDs] = hlflex
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, oldZ, op=MPI.MAX)

        hl = self.hLocal.getArray().copy()
        newZ = np.zeros(self.mpoints, dtype=np.float64) - 1.0e8
        newZ[self.locIDs] = hl
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, newZ, op=MPI.MAX)

        if MPIrank == 0:
            # Build regular grid for flexure calculation
            if self.xIndices is None:
                self._buildRegGrid()

            # Interpolate values on the flexural regular grid
            regOldZ = np.sum(self.regWeights * oldZ[self.regIDs][:, :], axis=1) / self.regSumWeights
            if len(self.regOnIDs) > 0:
                regOldZ[self.regOnIDs] = oldZ[self.regIDs[self.regOnIDs, 0]]
            regOldZ = regOldZ.reshape(self.reg_ny, self.reg_nx)

            regNewZ = np.sum(self.regWeights * newZ[self.regIDs][:, :], axis=1) / self.regSumWeights
            if len(self.regOnIDs) > 0:
                regNewZ[self.regOnIDs] = newZ[self.regIDs[self.regOnIDs, 0]]
            regNewZ = regNewZ.reshape(self.reg_ny, self.reg_nx)

            newh = flexure(regNewZ, regOldZ, self.reg_ny, self.reg_nx,
                           self.reg_yl, self.reg_xl, self.flex_rhos,
                           self.flex_rhoa, self.flex_eet, int(self.boundflex))
            rflexTec = newh - regNewZ

            # Interpolate back to goSPL mesh
            flexTec = self._regInterp(rflexTec)
        else:
            flexTec = None

        # Send flexural response globally
        flex = MPI.COMM_WORLD.bcast(flexTec, root=0)

        # Local flexural isostasy
        tmpFlex = flex[self.locIDs]
        if self.south == 1:
            tmpFlex[self.southPts] = 0.
        if self.east == 1:
            tmpFlex[self.eastPts] = 0.
        if self.north == 1:
            tmpFlex[self.northPts] = 0.
        if self.west == 1:
            tmpFlex[self.westPts] = 0.
        self.localFlex += tmpFlex
        tmp = self.hLocal.getArray().copy()

        # Update elevation
        self.hLocal.setArray(tmp + tmpFlex)
        self.dm.localToGlobal(self.hLocal, self.hGlobal)

        if MPIrank == 0 and self.verbose:
            print(
                "Compute Flexural Isostasy (%0.02f seconds)" % (process_time() - t0)
            )

        return

    def cptOrography(self):
        """
         Linear Theory of Orographic Precipitation following Smith & Barstad (2004).

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

        # Get elevations from the unstructured mesh structure
        hl = self.hLocal.getArray().copy()
        newZ = np.zeros(self.mpoints, dtype=np.float64) - 1.0e8
        newZ[self.locIDs] = hl
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, newZ, op=MPI.MAX)

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
            P_karot = ((self.oro_cw * 1j * sigma * hhat) / ((1 - (self.oro_hw * m * 1j)) * (1 + (sigma * self.oro_conv_time * 1j)) * (1 + (sigma * self.oro_fall_time * 1j))))

            # Inverse FFT, de-pad, convert units, add uniform rate
            oroRain = np.fft.ifft2(P_karot)
            oroRain = np.real(oroRain[pad:-pad, pad:-pad])
            oroRain *= 3600.  # mm hr-1
            oroRain += self.oro_precip_base
            # Precipitation rate must be a value greater than minimum precipitation/runoff to avoid errors when precip_rate <= 0
            oroRain[oroRain <= 0] = self.oro_precip_min
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
