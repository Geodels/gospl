import os
import gc
import sys
import petsc4py
import numpy as np
import numpy_indexed as npi

from mpi4py import MPI
from time import process_time

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import setMaxNb
    from gospl._fortran import marineCoeff
    from gospl._fortran import setHillslopeCoeff

petsc4py.init(sys.argv)
MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()
MPIcomm = petsc4py.PETSc.COMM_WORLD


class SEDMesh(object):
    """
    Performing surface evolution induced by considered processes:

     - hillslope processes
     - rivers using stream power law
    """

    def __init__(self, *args, **kwargs):

        # Petsc vectors
        self.tmp = self.hGlobal.duplicate()
        self.tmpL = self.hLocal.duplicate()
        self.Qs = self.hGlobal.duplicate()
        self.QsL = self.hLocal.duplicate()
        self.vSed = self.hGlobal.duplicate()
        self.vSedLocal = self.hLocal.duplicate()

        self.maxnb = setMaxNb(self.npoints)
        self.scaleIDs = np.zeros(self.npoints)
        self.scaleIDs[self.lIDs] = 1.0

        return

    def sedChange(self):
        """
        Perform erosion deposition changes.

        This function contains methods for the following operations:

         - stream induced deposition diffusion
         - hillslope diffusion
        """

        # Compute Continental Sediment Deposition
        self.seaQs = np.zeros(self.npoints, dtype=np.float64)
        self.vSedLocal.copy(result=self.QsL)
        self.dm.localToGlobal(self.QsL, self.Qs)

        # Define points in the marine environment for filled surface
        self.fillSeaID = np.where(self.hFill <= self.sealevel)[0]

        perc = 1.0
        iters = 0
        fill = False
        totQs = self.Qs.sum()
        while self.Qs.sum() > 0.0 and iters < 100:
            self._continentalDeposition()
            if perc < 1.0e-2:
                self.Qs.set(0.0)
                self.QsL.set(0.0)
            perc = self.Qs.sum() / totQs
            if self.Qs.sum() > 0.0:
                if perc >= 1.0e-2:
                    self.flowAccumulation(filled=False)
                else:
                    fill = True
                self._moveFluxes(filled=fill)
                iters += 1
                if iters == 100 and MPIrank == 0:
                    print(
                        "Continental sediment transport not converging; decrease time step",
                        flush=True,
                    )
            if MPIrank == 0 and self.verbose:
                print(
                    "Remaining percentage to transport: %0.01f " % (perc * 100.0),
                    flush=True,
                )

        # Compute Marine Sediment Deposition
        self.QsL.setArray(self.seaQs)
        self.dm.localToGlobal(self.QsL, self.Qs)
        self._marineDeposition()

        # Compute Hillslope Diffusion Law
        h = self.hLocal.getArray().copy()
        self.seaID = np.where(h <= self.sealevel)[0]
        self._hillSlope()

        del h
        gc.collect()

        return

    def _moveFluxes(self, filled=False):
        """
        Move continental sediment fluxes downstream based on depressions.
        """

        t0 = process_time()
        self.Qs.copy(result=self.tmp)

        # Transport sediment volume in m3 per year
        if filled:
            self._solve_KSP(False, self.fillMat, self.tmp, self.Qs)
        else:
            self._solve_KSP(False, self.wMat, self.tmp, self.Qs)

        # Update local vector
        self.dm.globalToLocal(self.Qs, self.QsL, 1)

        if MPIrank == 0 and self.verbose:
            print(
                "Move Sediment Downstream (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        return

    def getSedFlux(self):
        """
        Compute sediment flux in cubic metres per year.
        """

        t0 = process_time()

        # Get erosion rate (m/yr) to volume
        self.Eb.copy(result=self.tmp)
        self.tmp.pointwiseMult(self.tmp, self.areaGlobal)

        # Get the volume of sediment transported in m3 per year
        if self.tNow == self.rStart:
            self._solve_KSP(False, self.wMat, self.tmp, self.vSed)
        else:
            self._solve_KSP(False, self.wMat, self.tmp, self.vSed)

        # Update local vector
        self.dm.globalToLocal(self.vSed, self.vSedLocal, 1)
        if MPIrank == 0 and self.verbose:
            print(
                "Update Sediment Load (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        return

    def _continentalDeposition(self):
        """
        Perform continental deposition from incoming river flux.
        """

        t0 = process_time()

        # Get the marine volumetric sediment rate (m3 / yr) to distribute
        # during the time step
        tmp = self.QsL.getArray().copy()

        # Convert in volume (m3) for considered timestep
        # We only consider the nodes that are their own receivers
        Qs = np.zeros(self.npoints, dtype=np.float64)
        Qs[self.pitPts] = tmp[self.pitPts] * self.dt

        # To avoid counting twice the volume on the partition
        # boundaries we removed ghost nodes volumes
        Qs = np.multiply(Qs, self.scaleIDs)

        # Do not take ocean nodes they will be updated later
        self.seaQs[self.fillSeaID] += Qs[self.fillSeaID] / self.dt
        Qs[self.fillSeaID] = 0.0

        # Get the sediment volume available for transport and deposition
        # globally
        gQ = Qs[self.lgIDs]
        gQ[self.outIDs] = 0.0
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, gQ, op=MPI.MAX)

        # Find sediment load originally in a depression which are now
        # able to flow downslope
        ids = np.where(np.logical_and(self.pits[:, 0] == -1, gQ > 0))[0]
        newgQ = np.zeros(self.gpoints, dtype=np.float64)
        newgQ[ids] = gQ[ids]
        gQ[ids] = 0.0

        # Get the cumulative volume for each depression
        groupPits = npi.group_by(self.pits[:, 0])
        outNb, depFill = groupPits.sum(gQ)

        # Add the excess volume for each pit that needs to be distributed
        excess = np.where(depFill - self.pitVol > 0.0)[0]
        newgQ[self.outFlows[excess]] += depFill[excess] - self.pitVol[excess]
        newQs = newgQ[self.glIDs]

        # Scale the deposition based on available volume
        with np.errstate(divide="ignore", over="ignore"):
            scaleV = np.divide(
                depFill,
                self.pitVol,
                out=np.zeros_like(self.pitVol),
                where=self.pitVol != 0,
            )
        scaleV[scaleV > 1.0] = 1.0

        # Update available volume on each pit
        self.pitVol -= depFill
        self.pitVol[self.pitVol < 0.0] = 0.0

        # Deposit sediment in depressions based on the volumetric scaling
        h = self.hLocal.getArray().copy()
        dep = (self.hFill - h) * scaleV[self.pits[self.glIDs, 0]]

        # Update PETSc vectors
        self.QsL.setArray(newQs / self.dt)
        self.dm.localToGlobal(self.QsL, self.Qs)

        self.tmpL.setArray(dep)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self.cumED.axpy(1.0, self.tmp)
        self.hGlobal.axpy(1.0, self.tmp)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal, 1)
        self.dm.globalToLocal(self.hGlobal, self.hLocal, 1)

        if self.stratNb > 0:
            self._deposeStrat()

        del newQs, excess, newgQ, dep, h, scaleV, groupPits
        del gQ, ids, Qs, tmp, outNb, depFill
        gc.collect()

        if MPIrank == 0 and self.verbose:
            print(
                "Compute Continental Deposition (%0.02f seconds)"
                % (process_time() - t0),
                flush=True,
            )

        return

    def _marineDeposition(self):
        """
        Perform marine deposition from incoming river continental flux.
        """

        t0 = process_time()

        if self.sedimentK == 0.0:
            return

        # Get the marine volumetric sediment rate (m3 / yr) to diffuse
        # during the time step as suspended material...
        tmp = self.QsL.getArray().copy()
        Qs = np.zeros(self.npoints, dtype=np.float64)

        # Convert in volume (m3) for considered timestep
        Qs[self.fillSeaID] = tmp[self.fillSeaID] * self.dt

        # Diffusion matrix construction
        Cd = np.zeros(self.npoints, dtype=np.float64)
        Cd[self.fillSeaID] = self.sedimentK

        # From the distance to coastline define the upper limit
        # of the shelf to ensure a maximum slope angle
        if self.vtkMesh is not None:
            toplimit = self.sealevel - self.coastDist * 1.0e-4
            ids = self.coastDist < 2.0 * self.edgeMax
            toplimit[ids] = self.sealevel - self.coastDist[ids] * 1.0e-5
            ids = self.coastDist < self.edgeMax
            toplimit[ids] = self.sealevel
        else:
            toplimit = np.full((self.npoints), 0.9 * self.sealevel)

        # Define maximum deposition thicknesses and initialise
        # cumulative deposits
        h0 = self.hLocal.getArray().copy()
        maxDep = toplimit - h0
        maxDep[maxDep < 0.0] = 0.0
        cumDep = np.zeros(self.npoints, dtype=np.float64)

        # Build suspended sediment volume per unit area (m) vector
        self.tmpL.setArray(Qs)
        self.dm.localToGlobal(self.tmpL, self.Qs)
        maxSedVol = self.Qs.sum()
        self.Qs.pointwiseDivide(self.Qs, self.areaGlobal)

        diffCoeffs = marineCoeff(self.npoints, Cd * self.dt)
        self.Diff = self._matrix_build_diag(diffCoeffs[:, 0])
        indptr = np.arange(0, self.npoints + 1, dtype=petsc4py.PETSc.IntType)

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
                petsc4py.PETSc.InsertMode.INSERT_VALUES,
            )
            tmpMat.assemblyEnd()
            self.Diff += tmpMat
            tmpMat.destroy()
        del ids, indices, indptr, diffCoeffs, Cd
        gc.collect()

        iters = 0
        remainPerc = 1.0
        while iters < 100 and remainPerc > max(0.1, self.frac_fine):

            # Get erosion values for considered time step
            self._solve_KSP(False, self.Diff, self.Qs, self.tmp)

            # Find overfilled nodes
            self.dm.globalToLocal(self.tmp, self.tmpL, 1)
            dH = self.tmpL.getArray().copy()
            dH[dH < 0] = 0.0
            overDep = dH - maxDep
            overDep[overDep < 0] = 0.0
            overIDs = np.where(dH > maxDep)[0]

            # Update space both for cumulative and available depths
            cumDep += dH
            cumDep[overIDs] = toplimit[overIDs] - h0[overIDs]
            cumDep[cumDep < 0] = 0.0
            maxDep -= dH
            maxDep[maxDep < 0] = 0.0

            # Update sediment to diffuse
            Qs.fill(0.0)
            Qs[overIDs] = overDep[overIDs]

            # Update PETSc vector
            self.tmpL.setArray(Qs)
            self.dm.localToGlobal(self.tmpL, self.Qs)

            self.Qs.pointwiseMult(self.Qs, self.areaGlobal)
            sedVol = self.Qs.sum()
            remainPerc = sedVol / maxSedVol
            self.Qs.pointwiseDivide(self.Qs, self.areaGlobal)

            iters += 1
            if iters > 100 and MPIrank == 0:
                print(
                    "Sediment marine diffusion not converging; decrease time step",
                    flush=True,
                )

            if MPIrank == 0 and self.verbose:
                print(
                    "Remaining percentage to diffuse: %0.01f " % (remainPerc * 100.0),
                    flush=True,
                )

        # Update cumulative erosion/deposition and elevation
        cumDep[cumDep < 0] = 0.0
        self.tmpL.setArray(cumDep)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self.cumED.axpy(1.0, self.tmp)
        self.hGlobal.axpy(1.0, self.tmp)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal, 1)
        self.dm.globalToLocal(self.hGlobal, self.hLocal, 1)

        if self.stratNb > 0:
            self._deposeStrat()

        del h0, cumDep, dH, overDep, maxDep, Qs
        gc.collect()

        if MPIrank == 0 and self.verbose:
            print(
                "Compute Marine Sediment Diffusion (%0.02f seconds)"
                % (process_time() - t0),
                flush=True,
            )

        return

    def _hillSlope(self):
        """
        Perform hillslope diffusion.
        """

        if self.Cda == 0.0 and self.Cdm == 0.0:
            return

        t0 = process_time()

        # Diffusion matrix construction
        Cd = np.full(self.npoints, self.Cda, dtype=np.float64)
        Cd[self.seaID] = self.Cdm

        diffCoeffs = setHillslopeCoeff(self.npoints, Cd * self.dt)
        self.Diff = self._matrix_build_diag(diffCoeffs[:, 0])
        indptr = np.arange(0, self.npoints + 1, dtype=petsc4py.PETSc.IntType)

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
                petsc4py.PETSc.InsertMode.INSERT_VALUES,
            )
            tmpMat.assemblyEnd()
            self.Diff += tmpMat
            tmpMat.destroy()
        del ids, indices, indptr, diffCoeffs, Cd
        gc.collect()

        # Get elevation values for considered time step
        self.hGlobal.copy(result=self.hOld)
        self._solve_KSP(True, self.Diff, self.hOld, self.hGlobal)

        # Update cumulative erosion/deposition and elevation
        self.tmp.waxpy(-1.0, self.hOld, self.hGlobal)
        self.cumED.axpy(1.0, self.tmp)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal, 1)
        self.dm.globalToLocal(self.hGlobal, self.hLocal, 1)

        if self.stratNb > 0:
            self._deposeStrat()
            self.erodeStrat()
            self._elevStrat()

        if MPIrank == 0 and self.verbose:
            print(
                "Compute Hillslope Processes (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        return

    def _deposeStrat(self):
        """
        Add deposition on top of existing stratigraphic layer.
        """

        self.dm.globalToLocal(self.tmp, self.tmpL)
        depo = self.tmpL.getArray().copy()
        depo[depo < 0] = 0.0
        self.stratH[:, self.stratStep] += depo
        del depo
        gc.collect()

        return

    def erodeStrat(self):
        """
        Remove thickness from the stratigraphic pile.
        """

        self.dm.globalToLocal(self.tmp, self.tmpL)
        ero = self.tmpL.getArray().copy()
        ero[ero > 0] = 0.0

        # Nodes experiencing erosion
        nids = np.where(ero < 0)[0]
        if len(nids) == 0:
            return

        # Cumulative thickness for each node
        self.stratH[nids, 0] += 1.0e6
        cumThick = np.cumsum(self.stratH[nids, self.stratStep :: -1], axis=1)[:, ::-1]

        # Find nodes with no remaining stratigraphic thicknesses
        ids = np.where(-ero[nids] >= cumThick[:, 0])[0]
        self.stratH[nids[ids], : self.stratStep + 1] = 0.0

        # Erode remaining stratal layers
        if len(ids) < len(nids):
            ero[nids[ids]] = 0.0

            # Clear all stratigraphy points which are eroded
            cumThick[cumThick < -ero[nids].reshape((len(nids), 1))] = 0
            mask = (cumThick > 0).astype(int) == 0
            tmp = self.stratH[nids, : self.stratStep + 1]
            tmp[mask] = 0
            self.stratH[[nids], : self.stratStep + 1] = tmp

            # Update thickness of top stratigraphic layer
            eroIDs = np.bincount(np.nonzero(cumThick)[0]) - 1
            eroVal = cumThick[np.arange(len(nids)), eroIDs] + ero[nids]
            eroVal[ids] = 0.0
            self.stratH[nids, eroIDs] = eroVal

        self.stratH[nids, 0] -= 1.0e6
        self.stratH[self.stratH < 0] = 0.0

        return

    def _elevStrat(self):
        """
        Update current stratigraphic layer elevation.
        """

        self.stratZ[:, self.stratStep] = self.hLocal.getArray()

        return
