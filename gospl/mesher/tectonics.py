import os
import gc
import sys
import petsc4py
import numpy as np

from mpi4py import MPI
from scipy import spatial
from time import process_time

# from geomstats.geometry.hypersphere import Hypersphere

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import adveciioe
    from gospl._fortran import adveciioe2
    from gospl._fortran import advecupwind
    from gospl._fortran import fitedges
    from gospl._fortran import getrange
    from gospl._fortran import getfacevelocity

petsc4py.init(sys.argv)
MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()
MPIsize = petsc4py.PETSc.COMM_WORLD.Get_size()
MPIcomm = MPI.COMM_WORLD


class Tectonics(object):

    """
    This class defines how spherical mesh surface is changing given a plate reconstruction model.

    The horizontal displacement over time is done through a series of input files that defines the neighboring nodes that will be used for interpolation.

    .. note::

        The interpolation is based on nearest neighbour search based on the tree built from the spherical mesh and is weighted by distance. This is an efficient approach but it might not be the most accurate one.
    """

    def __init__(self):
        """
        The initialisation of the `EarthPlate` class.
        """

        self.hdisp = None
        self.tecNb = -1
        self.paleoZ = None
        self.minZ = None
        self.plateStep = False
        self.fiso = self.hGlobal.duplicate()

        return

    def getTectonics(self):
        """
        Finds the current tectonic regimes (horizontal and vertical) for the considered time interval.

        For horizontal displacements, the mesh variables will have to be first advected over the grid and then reinterpolated on the initial mesh coordinates.

        .. note::
            The approach here does not allow for mesh refinement in zones of convergence and thus can be limiting.

            Yet using a fixed mesh has one main advantage: the mesh and Finite Volume discretisation do not have to be rebuilt each time the mesh is advected.
        """

        if self.tecdata is None:
            return

        nb = self.tecNb
        if nb < len(self.tecdata) - 1:
            if self.tecdata.iloc[nb + 1, 0] < self.tNow + self.dt:
                nb += 1

        if nb > self.tecNb or nb == -1:
            if nb == -1:
                nb = 0

            self.tecNb = nb
            if nb < len(self.tecdata.index) - 1:
                timer = self.tecdata.iloc[nb, 1] - self.tecdata.iloc[nb, 0]
            else:
                timer = self.tEnd - self.tecdata.iloc[nb, 0]

            # Horizontal displacements
            if self.tecdata.iloc[nb, -1] != "empty":
                fname = self.tecdata.iloc[nb, -1][0] + ".npz"
                mdata = np.load(fname)
                key = self.tecdata.iloc[nb, -1][1]
                self.hdisp = mdata[key][self.locIDs, :]
                # In case of advection based on interpolation from plate position
                if self.advscheme == 0:
                    self._readAdvectionData(mdata[key], timer)
                    self.plateStep = True
                    self.plateTimer = self.tecdata.iloc[nb, 1]
                else:
                    # Get the velocity from the input file.
                    nodeVel = np.zeros((self.lpoints, 3))
                    if self.flatModel:
                        nodeVel[:, :2] = self.hdisp[:, :2]
                    else:
                        nodeVel = self.hdisp.copy()
                    # Store velocity on voronoi edges and dot product based on
                    # vecocity vector and face normals.
                    getfacevelocity(self.lpoints, nodeVel)
            else:
                self.hdisp = None

            # Vertical displacements
            if self.tecdata.iloc[nb, 2] != "empty":
                fname = self.tecdata.iloc[nb, 2][0] + ".npz"
                mdata = np.load(fname)
                key = self.tecdata.iloc[nb, 2][1]
                self.upsub = mdata[key][self.locIDs]
            else:
                self.upsub = None

            # Paleo-elevation fitting
            self.minZ = None
            self.paleoZ = None
            if self.tecdata.iloc[nb, 3] != "empty":
                fname = self.tecdata.iloc[nb, 3][0] + ".npz"
                mdata = np.load(fname)
                key = self.tecdata.iloc[nb, 3][1]
                if len(self.tecdata.iloc[nb, 3]) == 3:
                    self.minZ = self.tecdata.iloc[nb, 2][2]
                    self.paleoZ = mdata[key][self.locIDs]

            del mdata

        # Perform advection based on the flow-Implicit/Outflow-Explicit
        # or a `first-order upwind implicitly scheme.
        if self.hdisp is not None and self.advscheme > 0:
            self._varAdvector()

        return

    def _buildAdvecMat(self, iioe, lCoeffs, rCoeffs=None):
        """
        Create the advection matrix.
        """

        if self.flatModel:
            lCoeffs[self.idBorders, 1:] = 0.0
            lCoeffs[self.idBorders, 0] = 1.0
            if iioe:
                rCoeffs[self.idBorders, 1:] = 0.0
                rCoeffs[self.idBorders, 0] = 1.0

        advMat_left = self._matrix_build_diag(lCoeffs[:, 0])
        if iioe:
            advMat_right = self._matrix_build_diag(rCoeffs[:, 0])
        else:
            advMat_right = None

        indptr = np.arange(0, self.lpoints + 1, dtype=petsc4py.PETSc.IntType)

        for k in range(0, self.maxnb):
            tmpMat = self._matrix_build()
            indices = self.FVmesh_ngbID[:, k].copy()
            data = lCoeffs[:, k + 1]
            ids = np.nonzero(data == 0.0)
            indices[ids] = ids
            tmpMat.assemblyBegin()
            tmpMat.setValuesLocalCSR(
                indptr,
                indices.astype(petsc4py.PETSc.IntType),
                data,
            )
            tmpMat.assemblyEnd()
            advMat_left += tmpMat
            tmpMat.destroy()
            if iioe:
                tmpMat = self._matrix_build()
                indices = self.FVmesh_ngbID[:, k].copy()
                data = rCoeffs[:, k + 1]
                ids = np.nonzero(data == 0.0)
                indices[ids] = ids
                tmpMat.assemblyBegin()
                tmpMat.setValuesLocalCSR(
                    indptr,
                    indices.astype(petsc4py.PETSc.IntType),
                    data,
                )
                tmpMat.assemblyEnd()
                advMat_right += tmpMat
                tmpMat.destroy()
        if iioe:
            return advMat_left, advMat_right
        else:
            return advMat_left

    def _varAdvector(self):
        """
        Perform the advection of elevation and erosion by solving the advection equation based on the finite volume space discretization and the semi-implicit discretization in time. The approach is based on a Inflow-Implicit/Outflow-Explicit scheme following the work from `Mikula & Ohlberger, 2014 <https://www.math.sk/mikula/mo-FVCA6.pdf>`_ or a `first-order upwind implicitly scheme <https://www.sciencedirect.com/science/article/pii/S0168927414001032>`_ depending on the user configuration.

        .. note::

            Its basic idea is that outflow from a cell is treated explicitly while inflow is treated implicitly.

            Since the matrix of the system is determined by the inflow fluxes it is an M-matrix yielding favourable solvability and stability properties.

        .. important::

            The method allows large time steps without losing stability and not deteriorating precision.

            It is formally second order accurate in space and time for 1D advection problems with variable velocity and numerical experiments indicates its second order accuracy for smooth solutions in general.

        .. note::

            Velocity at the face is taken to be the linear interpolation for each vertex (in a vertex-centered discretisation the dual of the delaunay triangulation (i.e. the voronoi mesh has its edges on the middle of the nodes edges)

            Similarly we consider that the advected variable at the face is defined by linear interpolation from each connected vertex.
        """

        t0 = process_time()
        iioe = True
        if self.advscheme == 1:
            iioe = False

        # Advection matrix construction
        if iioe:
            if self.advscheme == 3:
                # Minimum and maximum in a local neighborhood
                hL = self.hLocal.getArray().copy()
                hmin, hmax = getrange(self.lpoints, hL)
            nbOut, lCoeffs, rCoeffs = adveciioe(self.lpoints, self.dt)
            advMat_left, advMat_right = self._buildAdvecMat(iioe, lCoeffs, rCoeffs)
        else:
            lCoeffs = advecupwind(self.lpoints, self.dt)
            advMat_left = self._buildAdvecMat(iioe, lCoeffs)

        # Advect elevations
        if iioe:
            # Inflow-Implicit/Outflow-Explicit Scheme 1
            advMat_right.mult(self.hGlobal, self.tmp1)
            self._solve_KSP(True, advMat_left, self.tmp1, self.tmp)
            # Inflow-Implicit/Outflow-Explicit Scheme 2
            if self.advscheme == 3:
                self.dm.globalToLocal(self.tmp, self.tmpL)
                newh = self.tmpL.getArray()
                diffmax = newh - hmax
                diffmax[diffmax < 0] = 0.
                diffmin = newh - hmin
                diffmin[diffmin > 0] = 0.
                diff = np.abs(diffmax) + np.abs(diffmin)
                self.tmpL.setArray(diff)
                self.dm.localToGlobal(self.tmpL, self.tmp)
                excess = self.tmp.sum()
                if excess > 0.:
                    lCoeffs, rCoeffs = adveciioe2(self.lpoints, self.dt, nbOut, hL, hmin, hmax)
                    advMat_left, advMat_right = self._buildAdvecMat(iioe, lCoeffs, rCoeffs)
                    advMat_right.mult(self.hGlobal, self.tmp1)
                    self._solve_KSP(True, advMat_left, self.tmp1, self.tmp)
        else:
            # Upwind scheme with potentially excessive diffusion solved implicitly
            self._solve_KSP(True, advMat_left, self.hGlobal, self.tmp)

        # Update elevations
        self.tmp.copy(result=self.hGlobal)
        self.dm.globalToLocal(self.hGlobal, self.hLocal)
        if self.flatModel:
            hL = self.hLocal.getArray().copy()
            hL[self.idBorders] = -1.e8
            nhL = fitedges(hL)
            self.hLocal.setArray(nhL)
            self.dm.localToGlobal(self.hLocal, self.hGlobal)

        # Advect erosion deposition
        if iioe:
            if self.advscheme == 3:
                self.dm.globalToLocal(self.cumED, self.cumEDLocal)
                # Minimum and maximum in a local neighborhood
                edL = self.cumEDLocal.getArray().copy()
                edmin, edmax = getrange(self.lpoints, edL)
            # Inflow-Implicit/Outflow-Explicit Scheme 1
            advMat_right.mult(self.cumED, self.tmp1)
            self._solve_KSP(True, advMat_left, self.tmp1, self.tmp)
            # Inflow-Implicit/Outflow-Explicit Scheme 2
            if self.advscheme == 3:
                self.dm.globalToLocal(self.tmp, self.tmpL)
                newed = self.tmpL.getArray()
                diffmax = newed - edmax
                diffmax[diffmax < 0] = 0.
                diffmin = newed - edmin
                diffmin[diffmin > 0] = 0.
                dh = np.abs(diffmax) + np.abs(diffmin)
                self.tmpL.setArray(dh)
                self.dm.localToGlobal(self.tmpL, self.tmp)
                excess = self.tmp.sum()
                if excess > 0.:
                    lCoeffs, rCoeffs = adveciioe2(self.lpoints, self.dt, nbOut, edL, edmin, edmax)
                    advMat_left, advMat_right = self._buildAdvecMat(iioe, lCoeffs, rCoeffs)
                    advMat_right.mult(self.cumED, self.tmp1)
                    self._solve_KSP(True, advMat_left, self.tmp1, self.tmp)
        else:
            # Upwind scheme with potentially excessive diffusion solved implicitly
            self._solve_KSP(True, advMat_left, self.cumED, self.tmp)

        # Update erosion deposition
        self.tmp.copy(result=self.cumED)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal)
        if self.flatModel:
            edL = self.cumEDLocal.getArray().copy()
            edL[self.idBorders] = -1.e8
            nedL = fitedges(edL)
            self.cumEDLocal.setArray(nedL)
            self.dm.localToGlobal(self.cumEDLocal, self.cumED)

        # Advect flexural isostasy
        if self.flexOn:
            if iioe:
                if self.advscheme == 3:
                    self.tmpL.setArray(self.localFlex)
                    self.dm.localToGlobal(self.tmpL, self.fiso)
                    # Minimum and maximum in a local neighborhood
                    fimin, fimax = getrange(self.lpoints, self.localFlex)
                # Inflow-Implicit/Outflow-Explicit Scheme 1
                advMat_right.mult(self.fiso, self.tmp1)
                self._solve_KSP(True, advMat_left, self.tmp1, self.tmp)
                # Inflow-Implicit/Outflow-Explicit Scheme 2
                if self.advscheme == 3:
                    self.dm.globalToLocal(self.tmp, self.tmpL)
                    newfi = self.tmpL.getArray()
                    diffmax = newfi - fimax
                    diffmax[diffmax < 0] = 0.
                    diffmin = newfi - fimin
                    diffmin[diffmin > 0] = 0.
                    dh = np.abs(diffmax) + np.abs(diffmin)
                    self.tmpL.setArray(dh)
                    self.dm.localToGlobal(self.tmpL, self.tmp)
                    excess = self.tmp.sum()
                    if excess > 0.:
                        lCoeffs, rCoeffs = adveciioe2(self.lpoints, self.dt, nbOut, self.localFlex, fimin, fimax)
                        advMat_left, advMat_right = self._buildAdvecMat(iioe, lCoeffs, rCoeffs)
                        advMat_right.mult(self.tmp, self.tmp1)
                        self._solve_KSP(True, advMat_left, self.tmp1, self.tmp)
            else:
                self.tmpL.setArray(self.localFlex)
                self.dm.localToGlobal(self.tmpL, self.fiso)
                # Upwind scheme with potentially excessive diffusion solved implicitly
                self._solve_KSP(True, advMat_left, self.fiso, self.tmp)

            # Update flexural isostasy
            self.tmp.copy(result=self.fiso)
            self.dm.globalToLocal(self.fiso, self.tmpL)
            self.localFlex = self.tmpL.getArray().copy()

        if MPIrank == 0 and self.verbose:
            print(
                "Compute Advection Processes (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        # Clean
        if iioe:
            advMat_right.destroy()
            del rCoeffs

        advMat_left.destroy()
        del edL, nedL, hL, nhL, lCoeffs
        # del edmin, edmax, hmin, hmax
        gc.collect()

        return

    def _readAdvectionData(self, hdisp, timer):
        """
        From a tectonic input file reads the horizontal displacements information, containing the displacement rates along each axis.

        .. note::
            A cKDTree is built with the advected coordinates and used to interpolate the mesh variables on the initial local mesh position.

            The interpolation is based on a weighting distance function accounting for the 3 closest advected vertices.

        :arg hdisp: displacement rates in 3D
        :arg timer: displacement time step in years
        """

        t0 = process_time()
        XYZ = self.mCoords + hdisp * timer

        # Build a tree with the advected nodes
        tree = spatial.cKDTree(XYZ, leafsize=10)

        # Query the distances from unstructured nodes
        distances, self.tec_IDs = tree.query(self.mCoords, k=3)

        # Inverse weighting distance...
        self.tec_weights = np.divide(
            1.0, distances, out=np.zeros_like(distances), where=distances != 0
        )
        self.tec_onIDs = np.where(distances[:, 0] == 0)[0]
        self.tec_sumw = np.sum(self.tec_weights, axis=1)

        del tree, distances

        if MPIrank == 0 and self.verbose:
            print(
                "Read advection data (%0.02f seconds)"
                % (process_time() - t0),
                flush=True,
            )

        return

    def _advectPlates(self):
        """
        Advects surface information based on plate evolution and performs interpolation.

        .. important::

            The interpolated values are the elevation, flexural response and cumulative erosion deposition and the interpolation is done on the spherical mesh on a single provessor. In case the stratigraphic information is also recorded then this is also interpolated.
        """

        # Send local elevation globally
        t0 = process_time()
        hl = self.hLocal.getArray().copy()
        gZ = np.zeros(self.mpoints, dtype=np.float64) - 1.0e8
        gZ[self.locIDs] = hl
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, gZ, op=MPI.MAX)

        # Send local erosion deposition globally
        edl = self.cumEDLocal.getArray().copy()
        gED = np.zeros(self.mpoints, dtype=np.float64) - 1.0e10
        gED[self.locIDs] = edl
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, gED, op=MPI.MAX)

        # Send local flexural isostasy globally
        if self.flexOn:
            gFI = np.zeros(self.mpoints, dtype=np.float64) - 1.0e10
            gFI[self.locIDs] = self.localFlex
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, gFI, op=MPI.MAX)

        if MPIrank == 0 and self.verbose:
            print(
                "Transfer local elevation, erosion deposition and flexural information globally (%0.02f seconds)"
                % (process_time() - t0),
                flush=True,
            )

        # Perform interpolation
        t0 = process_time()

        tmp = np.sum(self.tec_weights * gZ[self.tec_IDs], axis=1)
        nelev = np.divide(tmp, self.tec_sumw, out=np.zeros_like(self.tec_sumw), where=self.tec_sumw != 0)
        tmp = np.sum(self.tec_weights * gED[self.tec_IDs], axis=1)
        nerodep = np.divide(tmp, self.tec_sumw, out=np.zeros_like(self.tec_sumw), where=self.tec_sumw != 0)
        if self.flexOn:
            tmp = np.sum(self.tec_weights * gFI[self.tec_IDs], axis=1)
            nfi = np.divide(tmp, self.tec_sumw, out=np.zeros_like(self.tec_sumw), where=self.tec_sumw != 0)

        if len(self.tec_onIDs) > 0:
            nelev[self.tec_onIDs] = gZ[self.tec_IDs[self.tec_onIDs, 0]]
            nerodep[self.tec_onIDs] = gED[self.tec_IDs[self.tec_onIDs, 0]]
            if self.flexOn:
                nfi[self.tec_onIDs] = gFI[self.tec_IDs[self.tec_onIDs, 0]]

        if MPIrank == 0 and self.verbose:
            print(
                "Define local elevation and erosion/deposition after advection (%0.02f seconds)"
                % (process_time() - t0),
                flush=True,
            )

        self.hGlobal.setArray(nelev[self.glbIDs])
        self.dm.globalToLocal(self.hGlobal, self.hLocal)

        self.cumED.setArray(nerodep[self.glbIDs])
        self.dm.globalToLocal(self.cumED, self.cumEDLocal)

        if self.flexOn:
            self.localFlex = nfi[self.locIDs].copy()

        # Update stratigraphic record
        if self.stratNb > 0 and self.stratStep > 0:
            self._advectStrati()

        return

    def _findPts2Reduce(self):
        """
        Finds the nodes part of another partition that are required for interpolation.

        :return: Indices of needed nodes present in other partitions
        """

        # Global neighbours for each local partition
        lgNghbs = self.idNbghs[self.locIDs, :].flatten()

        # Global IDs required locally but part of another partition
        gids = lgNghbs[~np.in1d(lgNghbs, self.locIDs)]

        # Get all points that have to be transferred
        vals = np.zeros(self.mpoints, dtype=int)
        vals[gids] = 1
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, vals, op=MPI.MAX)

        return np.where(vals > 0)[0]

    def _advectStrati(self):
        """
        From advected stratigraphic information, retrieve points that are part of another partition and perform the interpolation.
        """

        # Get lobal point IDs that will need to be transferred globally
        t0 = process_time()
        redIDs = self._findPts2Reduce()
        if MPIrank == 0 and self.verbose:
            print(
                "Send stratigraphy information to local partition (%0.02f seconds)"
                % (process_time() - t0),
                flush=True,
            )

        # Transfer the variables accordingly and perform operation
        t0 = process_time()
        nghs = self.tec_IDs[self.locIDs]
        wgts = self.tec_weights[self.locIDs, :, None]
        sumw = self.tec_sumw[self.locIDs, :, None]
        onID = self.tec_onIDs[self.locIDs]

        self.stratH[:, : self.stratStep] = self._updateStratInfo(
            redIDs, self.stratH[:, : self.stratStep], sumw, onID, wgts, nghs
        )
        self.stratZ[:, : self.stratStep] = self._updateStratInfo(
            redIDs, self.stratZ[:, : self.stratStep], sumw, onID, wgts, nghs
        )
        self.phiS[:, : self.stratStep] = self._updateStratInfo(
            redIDs, self.phiS[:, : self.stratStep], sumw, onID, wgts, nghs
        )
        if MPIrank == 0 and self.verbose:
            print(
                "Define local stratigraphy variables after advection (%0.02f seconds)"
                % (process_time() - t0),
                flush=True,
            )
        del nghs, wgts, sumw, onID

        return

    def _updateStratInfo(self, reduceIDs, variable, sumw, onIDs, weights, nghs):
        """
        Update stratigraphic variables based on interpolated values.

        :arg reduceIDs: indices of nodes information that need to be available globally
        :arg variable: stratigraphic variable to interpolate
        :arg sumw: sum of weights for interpolation
        :arg onIDs: nodes that remain fixed
        :arg weights: weights for interpolation
        :arg nghs: neighbours for interpolation

        :return: nvar interpolated stratigraphic variable
        """

        # Reduce variable so that all values that need to be read from another partition can be
        vals = np.zeros((self.mpoints, self.stratStep), dtype=np.float64) - 1.0e8
        vals[self.locIDs, :] = variable
        redVals = vals[reduceIDs, :]
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, redVals, op=MPI.MAX)
        vals[reduceIDs, :] = redVals

        # Inverse weighting distance...
        tmp = np.sum(weights * vals[nghs, :], axis=1)
        nvar = np.divide(tmp, sumw, where=sumw != 0, out=np.zeros_like(tmp))
        if len(onIDs) > 0:
            nvar[onIDs, :] = vals[nghs[onIDs, 0], :]

        return nvar

    def updatePaleoZ(self):
        """
        Update surface information based on paleo-reconstruction.

        :arg paleoZ: paleo-elevation data
        :arg minZ: elevation threshold used for fitting
        """

        if self.plateStep:
            if self.plateTimer == self.tNow:
                self._advectPlates()

        if self.paleoZ is None:
            return

        t0 = process_time()

        # Send local elevation globally
        hl = self.hLocal.getArray().copy()
        ids = np.where(hl <= self.minZ)[0]
        hl[ids] = self.paleoZ[ids]

        # Fit simulated elevations to paleoelevation ones
        self.hLocal.setArray(hl)
        self.dm.localToGlobal(self.hLocal, self.hGlobal)

        if MPIrank == 0 and self.verbose:
            print(
                "Update model based on paleo-elevations (%0.02f seconds)"
                % (process_time() - t0),
                flush=True,
            )

        return
