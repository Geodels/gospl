import os
import gc
import sys
import petsc4py
import numpy as np
import numpy_indexed as npi

from mpi4py import MPI
from time import process_time

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import fillpit
    from gospl._fortran import setmaxnb
    from gospl._fortran import mfdreceivers
    from gospl._fortran import seaparams
    from gospl._fortran import distsea
    from gospl._fortran import distexcess
    from gospl._fortran import distocean
    from gospl._fortran import mfdrcvs
    from gospl._fortran import sethillslopecoeff

petsc4py.init(sys.argv)
MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()
MPIcomm = petsc4py.PETSc.COMM_WORLD


class SEDMesh(object):
    """
    This class encapsulates all the functions related to sediment transport, production and deposition. `gospl` has the ability to track three types of clastic sediment size and one type of carbonate (still under development). The following processes are considered:

    - inland river deposition in depressions and enclosed sea
    - marine deposition at river mouth
    - hillslope processes in both marine and inland areas

    .. note::
        Most of these functions are ran in parallel using the underlying PETSc library.

    """

    def __init__(self, *args, **kwargs):
        """
        The initialisation of `SEDMesh` class consists in the declaration of several PETSc vectors.
        """

        # Petsc vectors
        self.tmp = self.hGlobal.duplicate()
        self.tmpL = self.hLocal.duplicate()
        self.tmp1 = self.hGlobal.duplicate()
        self.Qs = self.hGlobal.duplicate()
        self.QsL = self.hLocal.duplicate()

        self.vSed = self.hGlobal.duplicate()
        self.vSedLocal = self.hLocal.duplicate()
        if self.stratNb > 0:
            if self.stratF is not None:
                self.vSedf = self.hGlobal.duplicate()
                self.vSedfLocal = self.hLocal.duplicate()
            if self.stratW is not None:
                self.vSedw = self.hGlobal.duplicate()
                self.vSedwLocal = self.hLocal.duplicate()
            if self.carbOn:
                self.vSedc = self.hGlobal.duplicate()
                self.vSedcLocal = self.hLocal.duplicate()

        # Get the maximum number of neighbours on the mesh
        maxnb = np.zeros(1, dtype=np.int64)
        maxnb[0] = setmaxnb(self.lpoints)
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, maxnb, op=MPI.MAX)
        self.maxnb = maxnb[0]

        return

    def getSedFlux(self):
        """
        This function computes sediment flux in cubic metres per year from incoming rivers. Like for the computation of the flow discharge and erosion rates, the sediment flux is solved by an implicit time integration method, the matrix system is the one obtained from the receiver distributions over the unfilled elevation mesh for the flow discharge (`fMat`).

        The PETSc *scalable linear equations solvers* (**KSP**) is used here again with an iterative method obtained from PETSc Richardson solver (`richardson`) with block Jacobian preconditioning (`bjacobi`).
        """

        t0 = process_time()

        # Multi-lithology case
        if self.stratNb > 0:
            # Coarse sediment
            # Get erosion rate (m/yr) to volume
            self.tmpL.setArray(self.thCoarse)
            self.dm.localToGlobal(self.tmpL, self.tmp)
            self.tmp.pointwiseMult(self.tmp, self.areaGlobal)
            # Get the volume of sediment transported in m3 per year
            self._solve_KSP(False, self.fMat, self.tmp, self.vSed)
            # Update local vector
            self.dm.globalToLocal(self.vSed, self.vSedLocal)

            # Fine sediment
            if self.stratF is not None:
                # Get erosion rate (m/yr) to volume
                self.tmpL.setArray(self.thFine)
                self.dm.localToGlobal(self.tmpL, self.tmp)
                self.tmp.pointwiseMult(self.tmp, self.areaGlobal)
                # Get the volume of sediment transported in m3 per year
                self._solve_KSP(False, self.fMat, self.tmp, self.vSedf)
                # Update local vector
                self.dm.globalToLocal(self.vSedf, self.vSedfLocal)

            # Weathered sediment
            if self.stratW is not None:
                # Get erosion rate (m/yr) to volume
                self.tmpL.setArray(self.thClay)
                self.dm.localToGlobal(self.tmpL, self.tmp)
                self.tmp.pointwiseMult(self.tmp, self.areaGlobal)
                # Get the volume of sediment transported in m3 per year
                self._solve_KSP(False, self.fMat, self.tmp, self.vSedw)
                # Update local vector
                self.dm.globalToLocal(self.vSedw, self.vSedwLocal)

            # Carbonate sediment
            if self.carbOn:
                # Get erosion rate (m/yr) to volume
                self.tmpL.setArray(self.thCarb)
                self.dm.localToGlobal(self.tmpL, self.tmp)
                self.tmp.pointwiseMult(self.tmp, self.areaGlobal)
                # Get the volume of sediment transported in m3 per year
                self._solve_KSP(False, self.fMat, self.tmp, self.vSedc)
                # Update local vector
                self.dm.globalToLocal(self.vSedc, self.vSedcLocal)

        else:
            # Get erosion rate (m/yr) to volume
            self.Eb.copy(result=self.tmp)
            self.tmp.pointwiseMult(self.tmp, self.areaGlobal)
            # Get the volume of sediment transported in m3 per year
            self._solve_KSP(False, self.fMat, self.tmp, self.vSed)
            # Update local vector
            self.dm.globalToLocal(self.vSed, self.vSedLocal)

        if MPIrank == 0 and self.verbose:
            print(
                "Update Sediment Load (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        return

    def _matrixDir(self):
        """
        This function defines the transport direction matrix for filled-limited elevations.

        .. note::

            Each matrix is built incrementally looping through the number of flow direction paths defined by the user. It proceeds by assembling a local Compressed Sparse Row (**CSR**) matrix to a global PETSc matrix.

            When setting up transport direction matrix in PETSc, we preallocate the non-zero entries of the matrix before starting filling in the values. Using PETSc sparse matrix storage scheme has the advantage that matrix-vector multiplication is extremely fast.

        The  matrix coefficients consist of weights (comprised between 0 and 1) calculated based on the number of downslope neighbours and proportional to the slope.

        """

        dirMat = self.zMat.copy()
        indptr = np.arange(0, self.lpoints + 1, dtype=petsc4py.PETSc.IntType)
        nodes = indptr[:-1]

        for k in range(0, self.flowDir):
            # Flow direction matrix for a specific direction
            tmpMat = self._matrix_build()
            data = self.lWght[:, k].copy()
            data[self.lRcv[:, k].astype(petsc4py.PETSc.IntType) == nodes] = 0.0
            tmpMat.assemblyBegin()
            tmpMat.setValuesLocalCSR(
                indptr,
                self.lRcv[:, k].astype(petsc4py.PETSc.IntType),
                data,
            )
            tmpMat.assemblyEnd()

            # Add the weights from each direction
            dirMat += tmpMat

            tmpMat.destroy()

        if self.memclear:
            del data, indptr, nodes
            gc.collect()

        # Store flow accumulation matrix
        self.dMat = dirMat.transpose().copy()

        dirMat.destroy()

        return

    def distributeSea(self, Qs):
        """
        From evaluation of enclosed seas (or endorheic basins with elevation below sea-level), this function computes the amount of incoming sediment which are deposited as well as the sediment volume able to flow downstream.

        The function calls the fortran subroutines `distsea` that takes the depressions parameters (sink IDs and overspilling node ID) to compute the different volumes and `distexcess` that moves excess sediment volumes to the overspilling nodes.

        It also calls the following *function*:

        - hierarchicalSink

        :arg Qs: numpy array of incoming sediment volume entering depressions

        :return: Qs (outflowing sediment volume numpy array)
        """

        ndepo = np.zeros(self.mpoints, dtype=np.float64)

        # Distribute sediments to next sink
        if MPIrank == 0:
            Qs, ndepo, outVol = distsea(
                Qs,
                self.gZ,
                self.sealevel,
                self.sGrp,
                self.sPits,
                self.sVol,
                self.sOut,
            )
        ndepo = MPI.COMM_WORLD.bcast(ndepo, root=0)
        if self.flatModel:
            ndepo[self.glBorders] = 0.0
        self.depo += ndepo
        self.gZ += ndepo

        # Was the elevation updated? if so update the downstream parameters
        if (ndepo > 0).any():
            self.hierarchicalSink(False, True, False, self.gZ)

        # Are there excess sea deposits to distribute?
        if MPIrank == 0:
            if (outVol > 0).any():
                Qs, ndepo = distexcess(
                    Qs,
                    self.gZ,
                    self.lFill,
                    self.lPits,
                    self.lVol,
                    outVol,
                    self.sOut,
                    self.lPits[:, 0],
                )
        Qs = MPI.COMM_WORLD.bcast(Qs, root=0)
        ndepo = MPI.COMM_WORLD.bcast(ndepo, root=0)
        if self.flatModel:
            ndepo[self.glBorders] = 0.0
            Qs[self.glBorders] = 0.0
        self.depo += ndepo
        self.gZ += ndepo

        if self.memclear:
            del ndepo
            gc.collect()

        return Qs

    def distributeContinent(self, Qs, sinkFlx):
        """
        This function takes the incoming sediment volume for each continental depression and finds the volume requires to fill the sinks up to overspilling as well as the volume available for distribution in downstream regions.

        The function calls the following *function*:

        - hierarchicalSink

        :arg Qs: numpy array of incoming sediment volume entering depressions
        :arg sinkFlx: numpy array of incoming sediment volume entering open ocean

        :return: Qs, sinkFlx (updated outflowing continental and open ocean sediment volumes numpy arrays)
        """

        self.tmpL.setArray(Qs[self.locIDs])
        self.dm.localToGlobal(self.tmpL, self.tmp)

        # Move to the first downstream nodes from the overspilling ones
        self.dMat.mult(self.tmp, self.tmp1)

        # Downstream distribution of the excess in continental regions
        self.hierarchicalSink(False, False, True, self.gZ)
        self._solve_KSP(False, self.fMat, self.tmp1, self.tmp)

        self.dm.globalToLocal(self.tmp, self.tmpL)
        nQs = np.zeros(self.mpoints)
        nQs[self.locIDs] = self.tmpL.getArray().copy()
        nQs[self.outIDs] = 0.0
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, nQs, op=MPI.MAX)
        if self.flatModel:
            nQs[self.glBorders] = 0.0

        # Update open ocean incoming sediment volume and downstream continental depression ones
        Qs = np.zeros(self.mpoints)
        if sinkFlx is not None:
            sinkFlx[self.sinkID] += nQs[self.sinkID]
        nQs[self.sinkID] = 0.0
        Qs[self.lnoRcv] = nQs[self.lnoRcv]
        nQs[self.lnoRcv] = 0.0

        if sinkFlx is not None:
            # Update sediment volume available in downstream portion of the rivers
            self.tmpL.setArray(nQs[self.locIDs])
            self.vSedLocal.axpy(1.0, self.tmpL)
            self.dm.localToGlobal(self.vSedLocal, self.vSed)
        else:
            # Update water volume available in downstream portion of the rivers
            self.tmpL.setArray(nQs[self.locIDs])
            self.FAL.axpy(1.0, self.tmpL)
            self.dm.localToGlobal(self.FAL, self.FAG)

        # Are the fluxes ending in continental sinks? If so find the overspilling volume if any.
        if (Qs > 0).any():
            ndepo = np.zeros(self.mpoints, dtype=np.float64)
            if MPIrank == 0:
                ids = np.where((Qs > 0.0) & (self.lPits[:, 1] == -1))[0]
                nQs = Qs[ids]
                Qs[ids] = 0.0
                outNb, depFill = self.grpSink.sum(Qs)
                with np.errstate(divide="ignore", over="ignore"):
                    scaleV = np.divide(
                        depFill,
                        self.lVol,
                        out=np.zeros_like(self.lVol),
                        where=self.lVol != 0,
                    )
                scaleV[scaleV > 1.0] = 1.0
                ndepo = (self.lFill - self.gZ) * scaleV[self.lPits[:, 0]]
                excess = np.where(depFill > self.lVol)[0]
                Qs.fill(0.0)
                Qs[self.lOut[excess]] = depFill[excess] - self.lVol[excess]
                Qs[ids] += nQs
                del ids, outNb, depFill, scaleV, excess

            Qs = MPI.COMM_WORLD.bcast(Qs, root=0)
            ndepo = MPI.COMM_WORLD.bcast(ndepo, root=0)
            if self.flatModel:
                ndepo[self.glBorders] = 0.0
                Qs[self.glBorders] = 0.0
            self.depo += ndepo
            self.gZ += ndepo
            if self.memclear:
                del ndepo

        if self.memclear:
            del nQs
            gc.collect()

        if sinkFlx is not None:
            return Qs, sinkFlx
        else:
            return Qs

    def _smoothSurface(self, topo):
        """
        This function uses a diffusion approach to smooth the bathymetry before performing marine deposition. In the marine realm, we assume that local slope plays a small role on sediment distribution and the smoothed surface is used to evaluate marine multiple flow directions.

        :arg topo: numpy array of initial surface to smooth

        :return: topo (smoothed surface numpy array)
        """

        # Diffusion matrix construction
        Cd = np.full(self.lpoints, 1.0e-4, dtype=np.float64)
        seaID = np.where(topo[self.locIDs] < self.sealevel)
        Cd[seaID] = 1.0e5  # this is an hardcoded value for now...

        diffCoeffs = sethillslopecoeff(self.lpoints, Cd * self.dt)
        diff = self._matrix_build_diag(diffCoeffs[:, 0])
        indptr = np.arange(0, self.lpoints + 1, dtype=petsc4py.PETSc.IntType)

        for k in range(0, self.maxnb):
            tmpMat = self._matrix_build()
            indices = self.FVmesh_ngbID[:, k].copy()
            data = diffCoeffs[:, k + 1]
            ids = np.nonzero(data == 0.0)
            indices[ids] = ids
            tmpMat.assemblyBegin()
            tmpMat.setValuesLocalCSR(
                indptr,
                indices.astype(petsc4py.PETSc.IntType),
                data,
            )
            tmpMat.assemblyEnd()
            diff += tmpMat
            tmpMat.destroy()

        # Get smoothed elevation values
        self.tmp1.setArray(topo[self.glbIDs])
        self._solve_KSP(False, diff, self.tmp1, self.tmp)
        self.dm.globalToLocal(self.tmp, self.tmpL)
        hl = self.tmpL.getArray().copy()
        topo.fill(0.0)
        topo[self.locIDs] = hl
        topo[self.outIDs] = -1.0e8
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, topo, op=MPI.MAX)

        diff.destroy()

        if self.memclear:
            del ids, indices, indptr, data, diffCoeffs, Cd, seaID, hl
            gc.collect()

        return topo

    def _marineDeposition(self, h, stype):
        """
        This function computes the maximum elevations of marine deposition based on clinoform slopes and antecedent shelf and continental slope physiographies.

        It also calls the following *private function*:

        - _smoothSurface

        :arg h: numpy array of local elevation values
        :arg stype: sediment type (integer)

        :return: smthZ, hclino (smoothed surface and maximum marine deposition elevations numpy arrays)
        """

        # Marine smoothed bathymetry
        smthZ = np.zeros(self.mpoints)
        smthZ[self.locIDs] = h
        smthZ[self.outIDs] = -1.0e8
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, smthZ, op=MPI.MAX)
        smthZ = self._smoothSurface(smthZ)

        # From the distance to coastline define the upper limit of the shelf to ensure a maximum slope angle
        clinoH = self.sealevel - 1.0 - self.coastDist * self.slps[stype] * 1.0e-5

        # Get the maximum marine deposition heights
        depo = smthZ[self.locIDs] + 500.0 - h
        shelfIDs = smthZ[self.locIDs] + 500.0 > clinoH
        depo[shelfIDs] = clinoH[shelfIDs] - h[shelfIDs]
        depo[depo < 0.0] = 0.0

        # From the marine deposition heights define the marine deposition elevations
        hclino = np.zeros(self.mpoints, dtype=np.float64) + 1.0e8
        hclino[self.locIDs] = h + depo
        hclino[self.outIDs] = 1.0e8
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, hclino, op=MPI.MIN)

        if self.memclear:
            del depo, clinoH, shelfIDs
            gc.collect()

        return smthZ, hclino

    def _diffuseOcean(self, zb, ndepo, stype):
        """
        For sediment reaching the marine realm, this function computes the related marine deposition diffusion. The approach is based on a linear diffusion.

        The diffusion equation is ran for the coarser sediment first and for the finest ones afterwards. This mimicks the standard behaviour observed in stratigraphic architectures where the fine fraction are generally transported over longer distances.

        :arg zb: numpy array of current elevation over which sediment diffusion is performed
        :arg ndepo: numpy array of incoming marine depositional thicknesses
        :arg stype: sediment type (integer)

        :return: ndepo, h (updated elevations and diffused depositions numpy arrays)
        """

        # Get diffusion coefficients based on sediment type
        sedK = self.sedimentK
        if sedK > 0.0:
            if stype == 1:
                sedK = self.sedimentKf
            if stype == 3:
                sedK = self.sedimentKw

        # Top elevations to diffuse
        h = zb + ndepo

        # Diffusion matrix construction
        Cd = np.zeros(self.lpoints, dtype=np.float64)
        Cd[self.seaID] = sedK
        scale = np.divide(ndepo[self.locIDs], ndepo[self.locIDs] + 500.0)
        Cd *= scale

        diffCoeffs = sethillslopecoeff(self.lpoints, Cd * self.dt)
        seaDiff = self._matrix_build_diag(diffCoeffs[:, 0])
        indptr = np.arange(0, self.lpoints + 1, dtype=petsc4py.PETSc.IntType)
        if self.memclear:
            del Cd, scale

        for k in range(0, self.maxnb):
            tmpMat = self._matrix_build()
            indices = self.FVmesh_ngbID[:, k].copy()
            data = diffCoeffs[:, k + 1]
            ids = np.nonzero(data == 0.0)
            indices[ids] = ids
            tmpMat.assemblyBegin()
            tmpMat.setValuesLocalCSR(
                indptr,
                indices.astype(petsc4py.PETSc.IntType),
                data,
            )
            tmpMat.assemblyEnd()
            seaDiff += tmpMat
            tmpMat.destroy()

        if self.memclear:
            del ids, data, indices, indptr, diffCoeffs
            gc.collect()

        # Perform diffusion equation
        self.tmp1.setArray(h[self.glbIDs])
        self._solve_KSP(False, seaDiff, self.tmp1, self.tmp)
        self.dm.globalToLocal(self.tmp, self.tmpL)
        hl = self.tmpL.getArray().copy()
        h = np.zeros(self.mpoints)
        h[self.locIDs] = hl
        h[self.outIDs] = -1.0e8
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, h, op=MPI.MAX)
        ndepo = h - zb
        ndepo[ndepo < 0.0] = 0.0

        seaDiff.destroy()
        if self.flatModel:
            h[self.glBorders] = zb[self.glBorders]

        if self.memclear:
            del hl
            gc.collect()

        return ndepo, h

    def _distributeOcean(self, sinkFlx, stype=0):
        """
        Main entry point to distribute incoming river sediment fluxes in the open ocean and perform diffusion of freshly deposited accumulations.

        The function calls the following *private functions*:

        - _marineDeposition
        - _diffuseOcean

        The diffusion equation is ran for the coarser sediment first and for the finest ones afterwards. This mimicks the standard behaviour observed in stratigraphic architectures where the fine fraction are generally transported over longer distances.

        :arg sinkFlx: numpy array of incoming sediment volume entering open ocean
        :arg stype: sediment type (integer)
        """

        exp = 1.0e-2
        flowDir = 8
        zb = self.gZ.copy()
        smthZ, hclino = self._marineDeposition(zb[self.locIDs], stype)

        maxDist = 2.0e6
        coastDist = np.zeros(self.mpoints, dtype=np.float64)
        coastDist[self.locIDs] = self.coastDist
        coastDist[self.outIDs] = 0.0
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, coastDist, op=MPI.MAX)
        ins = np.where(coastDist > maxDist)[0]
        if self.flatModel:
            ins = self.glBorders.copy()
        gZ = self.gZ.copy()
        smthZ[ins] = -3.0e4
        zsmth = smthZ.copy()
        sinkFlx = sinkFlx / float(self.diffstep)

        ndep = np.zeros(self.mpoints, dtype=np.float64)
        ndepo = np.zeros(self.mpoints, dtype=np.float64)
        for step in range(self.diffstep):
            ndepo.fill(0.0)
            sRcvs, sWghts = mfdrcvs(
                flowDir, exp, self.inIDs, smthZ[self.locIDs], -2.9e4
            )
            if self.flatModel:
                sRcvs[self.idBorders, :] = np.tile(self.idBorders, (flowDir, 1)).T
                sWghts[self.idBorders, :] = 0.0

            # Distribute open sea sediments
            sWght = np.zeros((self.mpoints, flowDir), dtype=np.float64)
            sWght[self.locIDs, :] = sWghts
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, sWght, op=MPI.MAX)

            sRcv = -np.ones((self.mpoints, flowDir), dtype=int)
            val = self.locIDs[sRcvs]
            ins = np.where(sRcvs == -1)[0]
            val[ins] = -1
            sRcv[self.locIDs, :] = val
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, sRcv, op=MPI.MAX)

            if MPIrank == 0:
                sid = np.argsort(smthZ)[::-1]
                topID = np.where(sinkFlx[sid] > 0)[0][0]
                vdepth = hclino - gZ
                vdepth[vdepth < 0.0] = 0.0
                ndepo = distocean(
                    flowDir, sid[topID:], sinkFlx, sRcv, sWght, self.marea, vdepth
                )
            ndepo = MPI.COMM_WORLD.bcast(ndepo, root=0)
            if self.flatModel:
                ndepo[self.glBorders] = 0.0
            ndep += ndepo
            ndep, gZ = self._diffuseOcean(zb, ndep, stype)
            if self.flatModel:
                ndep[self.glBorders] = 0.0
            smthZ = zsmth + ndep

        self.depo += ndep  # gZ - self.oldh
        self.gZ = gZ.copy()

        if self.memclear:
            del hclino, coastDist, gZ, smthZ, zsmth
            del sinkFlx, ndep, ndepo, sRcvs, sWghts, val
            if MPIrank == 0:
                del sid, vdepth
            gc.collect()

        return

    def _currentElevation(self):
        """
        Get current elevation.
        """

        hl = self.hLocal.getArray().copy()
        gZ = np.zeros(self.mpoints)
        gZ[self.locIDs] = hl
        gZ[self.outIDs] = -1.0e8
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, gZ, op=MPI.MAX)
        self.gZ = gZ.copy()
        self.oldh = gZ.copy()

        if self.memclear:
            del hl
            gc.collect()

        return gZ

    def _fill2D(self, minz, elev, fmax=1.0e6):
        """
        Fill pit for the 2D case.
        """

        nelev = elev.copy()
        nelev[self.glBorders] = minz - 1.0
        fill, pits = fillpit(minz, nelev, fmax)
        fill[self.glBorders] = elev[self.glBorders]
        pits[self.glBorders] = -1

        if self.memclear:
            del nelev
            gc.collect()

        return fill, pits

    def _nullitfyBorders(self, rcv, wght):
        """
        Set receivers and weights to zero for border points when running a 2D case.
        """

        rcv[self.idBorders, :] = np.tile(self.idBorders, (self.flowDir, 1)).T
        wght[self.idBorders, :] = 0.0

        return rcv, wght

    def _getDepressions(self):
        """
        Get continental depressions parameters.
        """

        sum_weight = np.sum(self.lWght, axis=1)
        tmpw = np.zeros(self.mpoints, dtype=np.float64)
        tmpw[self.locIDs] = sum_weight
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, tmpw, op=MPI.MAX)
        self.lnoRcv = tmpw == 0.0
        if self.memclear:
            del tmpw

        if MPIrank == 0:
            # Compute continental pit volumes
            self.grpSink = npi.group_by(self.lPits[:, 0])
            self.lGrp = self.grpSink.unique
            pitNb, self.lVol = self.grpSink.sum((self.lFill - self.gZ) * self.marea)
            _, outids, _ = np.intersect1d(self.lPits[:, 0], pitNb, return_indices=True)
            self.lOut = self.lPits[outids, 1]
            if self.memclear:
                del pitNb, outids

        return

    def _getCloseSea(self):
        """
        Get closed sea parameters.
        """

        # Set all nodes below sea-level as sinks
        self.sinkID = np.where(self.sFill < self.sealevel)[0]

        if MPIrank == 0:
            # Find closed seas
            closeIDs = np.where(
                (self.sFill > self.sealevel) & (self.gZ < self.sealevel)
            )[0]
            grpSink = npi.group_by(self.sPits[closeIDs, 0])
            self.sGrp = grpSink.unique
            self.sIns, self.sOut, self.sVol = seaparams(
                self.gZ, self.sealevel, self.sGrp, self.sPits
            )
            if self.memclear:
                del closeIDs, grpSink

        return

    def hierarchicalSink(self, offshore=False, land=False, limited=False, gZ=None):
        """
        This function defines the depression properties used to evaluate downstream sediment transport. To do so, the following three surfaces could be considered:

        1. a completely filled surface starting from the deep offshore regions
        2. a completely filled surface starting from the sea-level position
        3. a filled-limited surface starting from the sea-level position and limited by user defined value (to distribute sediment flux)

        In addition, we extract the volume of all depressions based on current elevation, depressionless one and voronoi cell areas. It also stores the spillover vertices indices for each of the depression. It is ran over the global mesh.

        .. note::

            This function uses the **numpy-indexed** library which contains functionality for indexed operations on numpy ndarrays and provides efficient vectorized functionality such as grouping and set operations.

        :arg offshore: boolean to compute sinks and associated parameters for surface 1
        :arg land: boolean to compute sinks and associated parameters for surface 2
        :arg limited: boolean to compute sinks and associated parameters for surface 3
        """

        # Get current state global elevations...
        if gZ is None:
            gZ = self._currentElevation()

        mingZ = -4.0e3
        if self.flatModel:
            mingZ = gZ.min()

        # Perform pit filling on process rank 0
        if MPIrank == 0:
            if offshore:
                if self.flatModel:
                    sFill, self.sPits = self._fill2D(mingZ, gZ, 1.0e6)
                else:
                    sFill, self.sPits = fillpit(-4.0e3, gZ, 1.0e6)
            if land:
                if self.flatModel:
                    lFill, self.lPits = self._fill2D(self.sealevel, gZ, 1.0e6)
                else:
                    lFill, self.lPits = fillpit(self.sealevel, gZ, 1.0e6)
            if limited:
                fill = self.sedfill
                if self.flowDist:
                    fill = self.waterfill
                if self.flatModel:
                    lFill, self.lPits = self._fill2D(self.sealevel, gZ, fill)
                else:
                    lFill, self.lPits = fillpit(self.sealevel, gZ, fill)
        else:
            if offshore:
                sFill = np.zeros(self.mpoints, dtype=np.float64)
            if land or limited:
                lFill = np.zeros(self.mpoints, dtype=np.float64)

        if offshore:
            self.sFill = MPI.COMM_WORLD.bcast(sFill, root=0)
            self._getCloseSea()

            # Define multiple flow directions for filled elevation
            self.sRcv, _, self.sWght = mfdreceivers(
                self.flowDir,
                self.flowExp,
                self.inIDs,
                self.sFill[self.locIDs],
                mingZ,
            )
            if self.flatModel:
                self.sRcv, self.sWght = self._nullitfyBorders(self.sRcv, self.sWght)

        if land or limited:
            self.lFill = MPI.COMM_WORLD.bcast(lFill, root=0)
            # Define multiple flow directions for filled elevation
            self.lRcv, _, self.lWght = mfdreceivers(
                self.flowDir,
                self.flowExp,
                self.inIDs,
                self.lFill[self.locIDs],
                self.sealevel,
            )
            if self.flatModel:
                self.lRcv, self.lWght = self._nullitfyBorders(self.lRcv, self.lWght)
            if land:
                self._matrixDir()
            else:
                self.matrixFlow(False)
            # Get depression parameters
            self._getDepressions()

        return

    def _getSedVol(self, stype):
        """
        Pick the relevant PETSc array for the specified sediment type.
        """

        if stype == 0:
            self.vSedLocal.copy(result=self.QsL)
        elif stype == 1:
            self.vSedfLocal.copy(result=self.QsL)
        elif stype == 2:
            self.vSedcLocal.copy(result=self.QsL)
        elif stype == 3:
            self.vSedwLocal.copy(result=self.QsL)
        self.dm.localToGlobal(self.QsL, self.Qs)

        return

    def _distributeSediment(self, stype=0):
        """
        Main entry point to distribute rivers' sediment fluxes over the surface. It accounts for both continental and marine deposition based on previously eroded sediments from rivers. The algorithm conserves sediment volumes and accounts for downstream deposition in both continental sink (endorheic basins), enclosed sea or lakes and open ocean environments.

        .. note::

            The deposition in depressions (both marine and continental) is performed based on the differences between sink volumes and incoming sediments ones and scale deposition accordingly. For the open seas, deposited sediment thicknesses are defined based on clinorform slopes and shelf and continental slopes physiographies before being diffused based on user diffusion coefficient.

        The function calls the following *private and non-private functions*:

        - hierarchicalSink
        - _getSedVol
        - distributeSea
        - distributeContinent
        - _distributeOcean

        :arg stype: sediment type (integer)
        """

        t0 = process_time()
        # Get downstream sinks in marine and continental regions
        self.hierarchicalSink(True, True, False, None)

        # Define the volume to distribute
        self._getSedVol(stype)

        # Get the volumetric sediment rate (m3/yr) to distribute during the time step and convert in volume (m3) for considered timestep
        lQs = self.QsL.getArray().copy()
        nQs = np.zeros(self.mpoints, dtype=np.float64)
        nQs[self.locIDs] = lQs * self.dt
        nQs[self.outIDs] = 0.0
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, nQs, op=MPI.MAX)
        Qs = np.zeros(self.mpoints, dtype=np.float64)
        if self.flatModel:
            Qs[self.glBorders] = 0.0
        Qs[self.pitPts] = nQs[self.pitPts]

        # Update all flux in the main sink
        sinkFlx = np.zeros(self.mpoints, dtype=np.float64)
        self.depo = np.zeros(self.mpoints, dtype=np.float64)
        sinkFlx[self.sinkID] += Qs[self.sinkID]
        Qs[self.sinkID] = 0.0

        step = 0
        while (Qs > 0).any():

            # Check if sediments are in closed sea
            inQs = np.where(Qs > 0)[0]
            insea = np.zeros(1, dtype=np.int64)
            if MPIrank == 0:
                if (self.sIns[inQs] > -1).any():
                    insea[0] = 1
            MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, insea, op=MPI.MAX)
            if self.memclear:
                del inQs
            if insea[0] > 0:
                # Distribute sediments in closed sea
                Qs = self.distributeSea(Qs)
                # Remove flux that flow into open sea
                sinkFlx[self.sinkID] += Qs[self.sinkID]
                Qs[self.sinkID] = 0.0
                if self.flatModel:
                    Qs[self.glBorders] = 0.0

            # At this point excess sediments are all in continental regions
            if (Qs > 0).any():
                Qs, sinkFlx = self.distributeContinent(Qs, sinkFlx)
                if self.flatModel:
                    Qs[self.glBorders] = 0.0

            step += 1
            if step > 20:
                if MPIrank == 0:
                    print("Forced Distribution to Ocean")
                sinkFlx += Qs
                break

        if MPIrank == 0 and self.verbose:
            print("Distribution step nb:", step)
            print(
                "Distribute Continental Sediments (%0.02f seconds)"
                % (process_time() - t0),
                flush=True,
            )

        # Distribute ocean sediments
        if (sinkFlx > 0).any():
            t0 = process_time()
            if self.flatModel:
                sinkFlx[self.glBorders] = 0.0
            self._distributeOcean(sinkFlx, stype)
            if MPIrank == 0 and self.verbose:
                print(
                    "Distribute Ocean Sediments (%0.02f seconds)"
                    % (process_time() - t0),
                    flush=True,
                )

        # Update cumed and elev
        self.tmpL.setArray(self.depo[self.locIDs])
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self.cumED.axpy(1.0, self.tmp)
        self.hGlobal.axpy(1.0, self.tmp)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal)
        self.dm.globalToLocal(self.hGlobal, self.hLocal)

        # Update stratigraphic layer parameters
        if self.stratNb > 0:
            self.deposeStrat(stype)
            self.elevStrat()

        return

    def sedChange(self):
        """
        This function is the main entry point to perform both continental and marine river-induced deposition. It calls the private function `_distributeSediment`.
        """

        if self.stratNb > 0:
            self._distributeSediment(stype=0)
            if self.carbOn:
                self._distributeSediment(stype=2)
            if self.stratF is not None:
                self._distributeSediment(stype=1)
            if self.stratW is not None:
                self._distributeSediment(stype=3)
        else:
            self._distributeSediment(stype=0)

        return

    def getHillslope(self):
        r"""
        This function computes hillslope processes. The code assumes that gravity is the main driver for transport and states that the flux of sediment is proportional to the gradient of topography.

        As a result, we use a linear diffusion law commonly referred to as **soil creep**:

        .. math::
          \frac{\partial z}{\partial t}= \kappa_{D} \nabla^2 z

        in which :math:`\kappa_{D}` is the diffusion coefficient and can be defined with different values for the marine and land environments (set with `hillslopeKa` and `hillslopeKm` in the YAML input file). It encapsulates, in a simple formulation, processes operating on superficial sedimentary layers. Main controls on variations of :math:`\kappa_{D}` include substrate, lithology, soil depth, climate and biological activity.
        """

        # Compute Hillslope Diffusion Law
        h = self.hLocal.getArray().copy()
        self.seaID = np.where(h <= self.sealevel)[0]
        self._hillSlope()

        # Update layer elevation
        if self.stratNb > 0:
            self.elevStrat()
        if self.memclear:
            del h
            gc.collect()

        return

    def _growCarbonates(self):
        """
        When carbonates is turned on update carbonate thicknesses based on fuzzy logic controls. The carbonate thicknesses created are uncompacted ones.

        .. warning::

            This function is a place order and will be updated in a future version of `gospl`.
        """

        # Limit fuzzy controllers to possible value range
        hl = self.sealevel - self.hLocal.getArray().copy()

        # TO DO: this will need to be updated with additional controls...
        carbH = np.zeros(self.lpoints, dtype=np.float64)
        validIDs = np.ones(self.lpoints, dtype=np.int64)
        for k in range(len(self.carbCtrl.controlKeys)):
            # Accomodation control
            if self.carbCtrl.controlKeys[k] == "depth":
                ids = np.where(
                    np.logical_or(
                        hl <= self.carbCtrl.controlBounds[k][0],
                        hl >= self.carbCtrl.controlBounds[k][1],
                    )
                )[0]
                validIDs[ids] = 0

        # Compute carbonate growth rates in metres per thousand years
        # and associated thicknesses
        ids = np.where(validIDs > 0)[0]
        for k in range(len(ids)):
            self.carbCtrl.carbControlSystem.input["depth"] = hl[ids[k]]
            # self.carbCtrl.carbControlSystem.input["temperature"] = 20.0
            self.carbCtrl.carbControlSystem.compute()
            growthRate = self.carbCtrl.carbControlSystem.output["growth"]
            carbH[ids[k]] = min(hl[ids[k]] - 1.0e-6, growthRate * self.dt * 1.0e-3)

        # Update top layers based on carbonate growth
        self.tmpL.setArray(carbH)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self.cumED.axpy(1.0, self.tmp)
        self.hGlobal.axpy(1.0, self.tmp)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal)
        self.dm.globalToLocal(self.hGlobal, self.hLocal)

        if self.stratNb > 0:
            # Update carbonate reef content in the stratigraphic layer
            self.deposeStrat(2)

        if self.memclear:
            del ids, validIDs, carbH, hl
            gc.collect()

        return

    def _hillSlope(self):
        r"""
        This function computes hillslope using a linear diffusion law commonly referred to as **soil creep**:

        .. math::
          \frac{\partial z}{\partial t}= \kappa_{D} \nabla^2 z

        in which :math:`\kappa_{D}` is the diffusion coefficient and can be defined with different values for the marine and land environments (set with `hillslopeKa` and `hillslopeKm` in the YAML input file).

        .. note::
            The hillslope processes in `gospl` are considered to be happening at the same rate for coarse and fine sediment sizes.

        """

        if self.Cda == 0.0 and self.Cdm == 0.0:
            return

        t0 = process_time()

        # Diffusion matrix construction
        Cd = np.full(self.lpoints, self.Cda, dtype=np.float64)
        Cd[self.seaID] = self.Cdm

        diffCoeffs = sethillslopecoeff(self.lpoints, Cd * self.dt)
        diffMat = self._matrix_build_diag(diffCoeffs[:, 0])
        indptr = np.arange(0, self.lpoints + 1, dtype=petsc4py.PETSc.IntType)

        for k in range(0, self.maxnb):
            tmpMat = self._matrix_build()
            indices = self.FVmesh_ngbID[:, k].copy()
            data = diffCoeffs[:, k + 1]
            ids = np.nonzero(data == 0.0)
            indices[ids] = ids
            tmpMat.assemblyBegin()
            tmpMat.setValuesLocalCSR(
                indptr,
                indices.astype(petsc4py.PETSc.IntType),
                data,
            )
            tmpMat.assemblyEnd()
            diffMat += tmpMat
            tmpMat.destroy()

        if self.memclear:
            del ids, indices, indptr, diffCoeffs, Cd
            gc.collect()

        # Get elevation values for considered time step
        self.hGlobal.copy(result=self.hOld)
        self._solve_KSP(True, diffMat, self.hOld, self.hGlobal)
        diffMat.destroy()

        # Update cumulative erosion/deposition and elevation
        self.tmp.waxpy(-1.0, self.hOld, self.hGlobal)
        self.cumED.axpy(1.0, self.tmp)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal)
        self.dm.globalToLocal(self.hGlobal, self.hLocal)

        if self.stratNb > 0:
            self.erodeStrat()
            if self.stratW is not None:
                self.deposeStrat(3)
            else:
                self.deposeStrat(0)

        if MPIrank == 0 and self.verbose:
            print(
                "Compute Hillslope Processes (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        return
