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
    This class encapsulates all the functions related to sediment transport, production and deposition.
    `gospl` has the ability to track two types of clastic sediment size and one type of carbonate (still
    under development). The following processes are considered:

    - inland river deposition in depressions
    - marine deposition at river mouth
    - hillslope processes in both marine and inland areas
    - sediment compaction as stratigraphic layers geometry and properties change

    .. note::
        All these functions are ran in parallel using the underlying PETSc library.

    """

    def __init__(self, *args, **kwargs):
        """
        The initialisation of `SEDMesh` class consists in the declaration of several PETSc vectors.
        """

        # Petsc vectors
        self.tmp = self.hGlobal.duplicate()
        self.tmpL = self.hLocal.duplicate()
        self.Qs = self.hGlobal.duplicate()
        self.QsL = self.hLocal.duplicate()
        self.vSed = self.hGlobal.duplicate()
        self.vSedLocal = self.hLocal.duplicate()
        if self.stratNb > 0:
            self.vSedf = self.hGlobal.duplicate()
            self.vSedfLocal = self.hLocal.duplicate()
            if self.carbOn:
                self.vSedc = self.hGlobal.duplicate()
                self.vSedcLocal = self.hLocal.duplicate()
        maxnb = np.zeros(1, dtype=np.int)
        maxnb[0] = setMaxNb(self.npoints)
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, maxnb, op=MPI.MAX)
        self.maxnb = maxnb[0]
        self.scaleIDs = np.zeros(self.npoints)
        self.scaleIDs[self.lIDs] = 1.0

        return

    def getSedFlux(self):
        """
        This function computes sediment flux in cubic metres per year from incoming rivers. Like for
        the computation of the flow discharge and erosion rates, the sediment flux is solved by an
        implicit time integration method, the matrix system is the one obtained from the receiver
        distributions over the unfilled elevation mesh for the flow discharge (`wMat`). The PETSc
        *scalable linear equations solvers* (**KSP**) is used here again with an iterative method
        obtained from PETSc Richardson solver (`richardson`) with block Jacobian preconditioning
        (`bjacobi`).
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
            self._solve_KSP(False, self.wMat, self.tmp, self.vSed)
            # Update local vector
            self.dm.globalToLocal(self.vSed, self.vSedLocal, 1)

            # Fine sediment
            # Get erosion rate (m/yr) to volume
            self.tmpL.setArray(self.thFine)
            self.dm.localToGlobal(self.tmpL, self.tmp)
            self.tmp.pointwiseMult(self.tmp, self.areaGlobal)
            # Get the volume of sediment transported in m3 per year
            self._solve_KSP(False, self.wMat, self.tmp, self.vSedf)
            # Update local vector
            self.dm.globalToLocal(self.vSedf, self.vSedfLocal, 1)

            # Carbonate sediment
            if self.carbOn:
                # Get erosion rate (m/yr) to volume
                self.tmpL.setArray(self.thCarb)
                self.dm.localToGlobal(self.tmpL, self.tmp)
                self.tmp.pointwiseMult(self.tmp, self.areaGlobal)
                # Get the volume of sediment transported in m3 per year
                self._solve_KSP(False, self.wMat, self.tmp, self.vSedc)
                # Update local vector
                self.dm.globalToLocal(self.vSedc, self.vSedcLocal, 1)

        else:
            # Get erosion rate (m/yr) to volume
            self.Eb.copy(result=self.tmp)
            self.tmp.pointwiseMult(self.tmp, self.areaGlobal)
            # Get the volume of sediment transported in m3 per year
            self._solve_KSP(False, self.wMat, self.tmp, self.vSed)
            # Update local vector
            self.dm.globalToLocal(self.vSed, self.vSedLocal, 1)

        if MPIrank == 0 and self.verbose:
            print(
                "Update Sediment Load (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        return

    def sedChange(self):
        """
        This function is the main entry point to perform both continental and marine river-induced
        deposition. It calls the private function `_sedChange`.
        """

        # Define points in the marine environment for filled surface
        self.fillSeaID = np.where(self.hFill <= self.sealevel)[0]

        if self.stratNb > 0:
            self._sedChange(stype=0)
            if self.carbOn:
                self._sedChange(stype=2)
            self._sedChange(stype=1)
        else:
            self._sedChange(stype=0)

        return

    def _sedChange(self, stype):
        """
        Deposition in depressions and the marine environments.

        This function contains methods for the following operations:

        - stream induced deposition in depression
        - marine river sediments diffusion
        - carbonate production from fuzzy logic

        .. note::

            When dealing with multiple depressions, the calculation of sediment flux
            over a landscape needs to be done carefully. The approach proposed here
            consists in iteratively filling the depressions as the sediment flux is
            transported downstream. To do so it *(1)* records the volume change in depressions
            if any, *(2)* updates elevation changes and *(3)* recomputes the flow discharge
            based on the `FAMesh` class functions.

        Once river sediments have reached the marine environment, a marine sediment-transport
        equations taking into account the tracking of the two different grain sizes is performed using a linear diffusion equation characterized by distinct transport coefficients.

        """

        # Compute continental sediment deposition
        self.seaQs = np.zeros(self.npoints, dtype=np.float64)
        if stype == 0:
            self.vSedLocal.copy(result=self.QsL)
        elif stype == 1:
            self.vSedfLocal.copy(result=self.QsL)
        elif stype == 2:
            self.vSedcLocal.copy(result=self.QsL)
        self.dm.localToGlobal(self.QsL, self.Qs)

        perc = 1.0
        minperc = 1.0e-3
        iters = 0
        fill = False
        totQs = self.Qs.sum()
        while self.Qs.sum() > 0.0 and iters < 100:
            self._continentalDeposition(stype, filled=fill)
            if perc < minperc:
                self.Qs.set(0.0)
                self.QsL.set(0.0)
            perc = self.Qs.sum() / totQs
            if self.Qs.sum() > 0.0:
                if perc >= minperc and iters < 99:
                    self.flowAccumulation(filled=False, limit=False)
                else:
                    fill = True
                self._moveFluxes(filled=fill)
                iters += 1
            if MPIrank == 0 and self.verbose:
                print(
                    "Remaining percentage to transport: %0.01f %d"
                    % (perc * 100.0, iters),
                    flush=True,
                )

        # Compute Marine Sediment Deposition
        self.QsL.setArray(self.seaQs)
        self.dm.localToGlobal(self.QsL, self.Qs)
        sedK = self.sedimentK
        if sedK > 0.0:
            if stype == 1:
                sedK = self.sedimentKf
            self._marineDeposition(stype, sedK)

        # Compute Fuzzy Logic Carbonate Growth
        if self.carbOn:
            self._growCarbonates()

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
            self._elevStrat()

        del h
        gc.collect()

        return

    def _growCarbonates(self):
        """
        When carbonates is turned on update carbonate thicknesses based on fuzzy logic controls.
        The carbonate thicknesses created are uncompacted ones.

        .. warning::

            This function is a place order and will be updated in a future version of `gospl`.
        """

        # Limit fuzzy controllers to possible value range
        hl = self.sealevel - self.hLocal.getArray().copy()

        # TO DO: this will need to be updated with additional controls...
        carbH = np.zeros(self.npoints, dtype=np.float64)
        validIDs = np.ones(self.npoints, dtype=np.int64)
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
        self.dm.globalToLocal(self.cumED, self.cumEDLocal, 1)
        self.dm.globalToLocal(self.hGlobal, self.hLocal, 1)

        if self.stratNb > 0:
            # Update carbonate reef content in the stratigraphic layer
            self._deposeStrat(2)

        del ids, validIDs, carbH, hl
        gc.collect()

        return

    def _moveFluxes(self, filled=False):
        """
        This function updates downstream continental sediment fluxes accounting for changes in
        elevation induced by deposition.

        As mentionned above, the PETSc *scalable linear equations solvers* (**KSP**) is used here
        for obtaining the solution.

        :arg filled: boolean to choose between unfilled and filled surface
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

    def _continentalDeposition(self, stype, filled=False):
        """
        Accounting for available volume for deposition in each depression, this function records
        the amount of sediment deposited inland and update the stratigraphic record accordingly.

        For each sediment type, freshly deposits are given the surface porosities given in the
        input file.

        If a depression is overfilled the excess sediment is added to the overspill node of the
        considered depression and is subsequently transported downstream in a subsequent iteration.

        During a given time step, the process described above will be repeated iteratively until
        all sediment transported are either deposited in depressions or are reaching the shoreline.

        :arg stype: sediment type (integer)
        :arg filled: boolean to choose between unfilled and filled surface
        """

        t0 = process_time()

        # Get the volumetric sediment rate (m3 / yr) to distribute
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

        # In case we are using the filled topography all continental rivers
        # drain into the ocean, therefore we can leave the function.
        if filled:
            return

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
            self._deposeStrat(stype)

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

    def _marineDeposition(self, stype, sedK):
        """
        For sediment reaching the coastline, this function computes the related marine
        deposition. The approach is based on a linear diffusion which is applied iteratively
        over a given time step until all the sediment have been diffused.

        The diffusion equation is ran for the coarser sediment first and for the finest one
        afterwards. This mimicks the standard behaviour observed in stratigraphic architectures
        where the fine fraction are generally transported over longer distances.

        The diffusion process stops when the sediment entering the ocean at river mouths are
        distributed and the accumulations on these specific nodes remains below water depth.

        :arg stype: sediment type (integer)
        :arg sedK: sediment diffusion coefficient for river transported sediment in the marine realm
        """

        t0 = process_time()

        # Get the marine volumetric sediment rate (m3 / yr) to diffuse
        # during the time step as suspended material...
        tmp = self.QsL.getArray().copy()
        Qs = np.zeros(self.npoints, dtype=np.float64)

        # Convert in volume (m3) for considered timestep
        Qs[self.fillSeaID] = tmp[self.fillSeaID] * self.dt

        # Diffusion matrix construction
        Cd = np.zeros(self.npoints, dtype=np.float64)
        Cd[self.fillSeaID] = sedK

        # From the distance to coastline define the upper limit
        # of the shelf to ensure a maximum slope angle
        if self.vtkMesh is not None:
            if stype == 1:
                toplimit = np.full((self.npoints), 0.9 * self.sealevel)
            else:
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
        while (
            iters < 500 and remainPerc > max(0.05, self.frac_fine) and maxSedVol > 0.0
        ):

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
            if iters >= 500 and MPIrank == 0:
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
            self._deposeStrat(stype)

        del h0, cumDep, maxDep, Qs
        if maxSedVol > 0.0:
            del dH, overDep
        gc.collect()

        if MPIrank == 0 and self.verbose:
            print(
                "Compute Marine Sediment Diffusion (%0.02f seconds)"
                % (process_time() - t0),
                flush=True,
            )

        return

    def _hillSlope(self):
        r"""
        This function computes hillslope using a linear diffusion law commonly referred to as **soil creep**:

        .. math::
          \frac{\partial z}{\partial t}= \kappa_{D} \nabla^2 z

        in which :math:`\kappa_{D}` is the diffusion coefficient and can be defined with different values for the marine and land environments (set with `hillslopeKa` and `hillslopeKm` in the YAML input file).

        .. note::
            The hillslope processes in `gospl` are considered to be happening at the same rate for coarse
            and fine sediment sizes.

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
            self.erodeStrat()
            self._deposeStrat(1)

        if MPIrank == 0 and self.verbose:
            print(
                "Compute Hillslope Processes (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        return

    def _deposeStrat(self, stype):
        """
        Add deposition on top of an existing stratigraphic layer. The following variables will be recorded:

        - thickness of each stratigrapic layer `stratH` accounting for both
          erosion & deposition events.
        - proportion of fine sediment `stratF` contains in each stratigraphic layer.
        - porosity of coarse sediment `phiS` in each stratigraphic layer computed at
          center of each layer.
        - porosity of fine sediment `phiF` in each stratigraphic layer computed at
          center of each layer.
        - proportion of carbonate sediment `stratC` contains in each stratigraphic layer
          if the carbonate module is turned on.
        - porosity of carbonate sediment `phiC` in each stratigraphic layer computed at
          center of each layer when the carbonate module is turned on.

        :arg stype: sediment type (integer)
        """

        self.dm.globalToLocal(self.tmp, self.tmpL)
        depo = self.tmpL.getArray().copy()
        depo[depo < 0] = 0.0
        fineH = self.stratH[:, self.stratStep] * self.stratF[:, self.stratStep]
        if self.carbOn:
            carbH = self.stratH[:, self.stratStep] * self.stratC[:, self.stratStep]
        self.stratH[:, self.stratStep] += depo
        ids = np.where(depo > 0)[0]

        if stype == 0:
            self.phiS[ids, self.stratStep] = self.phi0s
            self.stratF[ids, self.stratStep] = (
                fineH[ids] / self.stratH[ids, self.stratStep]
            )
            if self.carbOn:
                self.stratC[ids, self.stratStep] = (
                    carbH[ids] / self.stratH[ids, self.stratStep]
                )
                del carbH

        elif stype == 1:
            fineH[ids] += depo[ids]
            self.stratF[ids, self.stratStep] = (
                fineH[ids] / self.stratH[ids, self.stratStep]
            )
            self.phiF[ids, self.stratStep] = self.phi0f
            if self.carbOn:
                self.stratC[ids, self.stratStep] = (
                    carbH[ids] / self.stratH[ids, self.stratStep]
                )
                del carbH

        elif stype == 2:
            self.stratF[ids, self.stratStep] = (
                fineH[ids] / self.stratH[ids, self.stratStep]
            )
            carbH[ids] += depo[ids]
            self.stratC[ids, self.stratStep] = (
                carbH[ids] / self.stratH[ids, self.stratStep]
            )
            self.phiC[ids, self.stratStep] = self.phi0c
            del carbH

        # Cleaning arrays
        del depo, fineH, ids
        gc.collect()

        return

    def erodeStrat(self):
        """
        This function removes eroded sediment thicknesses from the stratigraphic pile.
        The function takes into account the porosity values of considered lithologies in
        each eroded stratigraphic layers.

        It follows the following assumptions:

        - Eroded thicknesses from stream power law and hillslope diffusion are considered
          to encompass both the solid and void phase.
        - Only the solid phase will be moved dowstream by surface processes.
        - The corresponding deposit thicknesses for those freshly eroded sediments correspond to
          uncompacted thicknesses based on the porosity at surface given from the input file.
        """

        self.dm.globalToLocal(self.tmp, self.tmpL)
        ero = self.tmpL.getArray().copy()
        ero[ero > 0] = 0.0

        # Nodes experiencing erosion
        nids = np.where(ero < 0)[0]
        if len(nids) == 0:
            del ero, nids
            gc.collect()
            return

        # Cumulative thickness for each node
        self.stratH[nids, 0] += 1.0e6
        cumThick = np.cumsum(self.stratH[nids, self.stratStep :: -1], axis=1)[:, ::-1]
        boolMask = cumThick < -ero[nids].reshape((len(nids), 1))
        mask = boolMask.astype(int)

        # Get fine sediment eroded from river incision
        thickF = (
            self.stratH[nids, 0 : self.stratStep + 1]
            * self.stratF[nids, 0 : self.stratStep + 1]
        )
        # From fine thickness extract the solid phase that is eroded
        thFine = thickF * (1.0 - self.phiF[nids, 0 : self.stratStep + 1])
        thFine = np.sum((thFine * mask), axis=1)

        # Get carbonate sediment eroded from river incision
        if self.carbOn:
            thickC = (
                self.stratH[nids, 0 : self.stratStep + 1]
                * self.stratC[nids, 0 : self.stratStep + 1]
            )
            # From carbonate thickness extract the solid phase that is eroded
            thCarb = thickC * (1.0 - self.phiC[nids, 0 : self.stratStep + 1])
            thCarb = np.sum((thCarb * mask), axis=1)
            # From sand thickness extract the solid phase that is eroded
            thickS = self.stratH[nids, 0 : self.stratStep + 1] - thickF - thickC
            thCoarse = thickS * (1.0 - self.phiS[nids, 0 : self.stratStep + 1])
            thCoarse = np.sum((thCoarse * mask), axis=1)
        else:
            # From sand thickness extract the solid phase that is eroded
            thickS = self.stratH[nids, 0 : self.stratStep + 1] - thickF
            thCoarse = thickS * (1.0 - self.phiS[nids, 0 : self.stratStep + 1])
            thCoarse = np.sum((thCoarse * mask), axis=1)

        # Clear all stratigraphy points which are eroded
        cumThick[boolMask] = 0.0
        tmp = self.stratH[nids, : self.stratStep + 1]
        tmp[boolMask] = 0
        self.stratH[nids, : self.stratStep + 1] = tmp

        # Erode remaining stratal layers
        # Get non-zero top layer number
        eroLayNb = np.bincount(np.nonzero(cumThick)[0]) - 1
        eroVal = cumThick[np.arange(len(nids)), eroLayNb] + ero[nids]

        # Get thickness of each sediment type eroded in the remaining layer
        self.thFine = np.zeros(self.npoints)
        # From fine thickness extract the solid phase that is eroded from this last layer
        tmp = (self.stratH[nids, eroLayNb] - eroVal) * self.stratF[nids, eroLayNb]
        thFine += tmp * (1.0 - self.phiF[nids, eroLayNb])
        # Define the uncompacted fine thickness that will be deposited dowstream
        self.thFine[nids] = thFine / (1.0 - self.phi0f)
        self.thFine[self.thFine < 0.0] = 0.0

        self.thCoarse = np.zeros(self.npoints)
        if self.carbOn:
            # From carb thickness extract the solid phase that is eroded from this last layer
            self.thCarb = np.zeros(self.npoints)
            tmp = (self.stratH[nids, eroLayNb] - eroVal) * self.stratC[nids, eroLayNb]
            thCarb += tmp * (1.0 - self.phiC[nids, eroLayNb])
            # Define the uncompacted carbonate thickness that will be deposited dowstream
            self.thCarb[nids] = thCarb / (1.0 - self.phi0c)
            self.thCarb[self.thCarb < 0.0] = 0.0
            # From sand thickness extract the solid phase that is eroded from this last layer
            tmp = self.stratH[nids, eroLayNb] - eroVal
            tmp *= 1.0 - self.stratC[nids, eroLayNb] - self.stratF[nids, eroLayNb]
            # Define the uncompacted sand thickness that will be deposited dowstream
            thCoarse += tmp * (1.0 - self.phiS[nids, eroLayNb])
            self.thCoarse[nids] = thCoarse / (1.0 - self.phi0s)
            self.thCoarse[self.thCoarse < 0.0] = 0.0
        else:
            # From sand thickness extract the solid phase that is eroded from this last layer
            tmp = self.stratH[nids, eroLayNb] - eroVal
            tmp *= 1.0 - self.stratF[nids, eroLayNb]
            # Define the uncompacted sand thickness that will be deposited dowstream
            thCoarse += tmp * (1.0 - self.phiS[nids, eroLayNb])
            self.thCoarse[nids] = thCoarse / (1.0 - self.phi0s)
            self.thCoarse[self.thCoarse < 0.0] = 0.0

        # Update thickness of top stratigraphic layer
        self.stratH[nids, eroLayNb] = eroVal
        self.stratH[nids, 0] -= 1.0e6
        self.stratH[self.stratH <= 0] = 0.0
        self.stratF[self.stratH <= 0] = 0.0
        self.phiS[self.stratH <= 0] = 0.0
        self.phiF[self.stratH <= 0] = 0.0
        if self.carbOn:
            self.stratC[self.stratH <= 0] = 0.0
            self.phiC[self.stratH <= 0] = 0.0

        self.thCoarse /= self.dt
        self.thFine /= self.dt
        if self.carbOn:
            self.thCarb /= self.dt

        del ero, nids, cumThick, boolMask, mask, tmp, eroLayNb, eroVal
        del thFine, thickF, thickS, thCoarse
        if self.carbOn:
            del thickC, thCarb

        gc.collect()

        return

    def _elevStrat(self):
        """
        This function updates the current stratigraphic layer elevation.
        """

        self.stratZ[:, self.stratStep] = self.hLocal.getArray()

        return

    def getCompaction(self):
        """
        This function computes the changes in sedimentary layers porosity and thicknesses due to
        compaction.

        .. note::

            We assume simple depth-porosiy relationships for each sediment type available in each
            layers.
        """

        t0 = process_time()
        topZ = np.vstack(self.hLocal.getArray())
        totH = np.sum(self.stratH[:, : self.stratStep + 1], axis=1)

        # Height of the sediment column above the center of each layer is given by
        cumZ = -np.cumsum(self.stratH[:, self.stratStep :: -1], axis=1) + topZ
        elev = np.append(topZ, cumZ[:, :-1], axis=1)
        zlay = np.fliplr(elev - np.fliplr(self.stratH[:, : self.stratStep + 1] / 2.0))

        # Compute lithologies porosities for each depositional layers
        # Get depth below basement
        depth = zlay - topZ

        # Now using depth-porosity relationships we compute the porosities
        # We assume that porosity cannot increase after unloading
        phiS = self.phi0s * np.exp(depth / self.z0s)
        phiS = np.minimum(phiS, self.phiS[:, : self.stratStep + 1])
        phiF = self.phi0f * np.exp(depth / self.z0f)
        phiF = np.minimum(phiF, self.phiF[:, : self.stratStep + 1])
        if self.carbOn:
            phiC = self.phi0c * np.exp(depth / self.z0c)
            phiC = np.minimum(phiC, self.phiC[:, : self.stratStep + 1])

        # Compute the solid phase in each layers
        tmpF = (
            self.stratH[:, : self.stratStep + 1] * self.stratF[:, : self.stratStep + 1]
        )
        tmpF *= 1.0 - self.phiF[:, : self.stratStep + 1]

        if self.carbOn:
            tmpC = (
                self.stratH[:, : self.stratStep + 1]
                * self.stratC[:, : self.stratStep + 1]
            )
            tmpC *= 1.0 - self.phiC[:, : self.stratStep + 1]
            tmpS = (
                self.stratC[:, : self.stratStep + 1]
                + self.stratF[:, : self.stratStep + 1]
            )
            tmpS = self.stratH[:, : self.stratStep + 1] * (1.0 - tmpS)
            tmpS *= 1.0 - self.phiS[:, : self.stratStep + 1]
            solidPhase = tmpC + tmpS + tmpF
        else:
            tmpS = self.stratH[:, : self.stratStep + 1] * (
                1.0 - self.stratF[:, : self.stratStep + 1]
            )
            tmpS *= 1.0 - self.phiS[:, : self.stratStep + 1]
            solidPhase = tmpS + tmpF

        # Get new layer thickness after porosity change
        tmpF = self.stratF[:, : self.stratStep + 1] * (
            1.0 - phiF[:, : self.stratStep + 1]
        )
        if self.carbOn:
            tmpC = self.stratC[:, : self.stratStep + 1] * (
                1.0 - phiC[:, : self.stratStep + 1]
            )
            tmpS = (
                1.0
                - self.stratF[:, : self.stratStep + 1]
                - self.stratC[:, : self.stratStep + 1]
            )
            tmpS *= 1.0 - phiS[:, : self.stratStep + 1]
            tot = tmpS + tmpC + tmpF
        else:
            tmpS = (1.0 - self.stratF[:, : self.stratStep + 1]) * (
                1.0 - phiS[:, : self.stratStep + 1]
            )
            tot = tmpS + tmpF
        ids = np.where(tot > 0.0)
        newH = np.zeros(tot.shape)
        newH[ids] = solidPhase[ids] / tot[ids]
        newH[newH <= 0] = 0.0
        phiS[newH <= 0] = 0.0
        phiF[newH <= 0] = 0.0
        if self.carbOn:
            phiC[newH <= 0] = 0.0

        # Update porosities in each sedimentary layer
        self.phiF[:, : self.stratStep + 1] = phiF
        self.phiS[:, : self.stratStep + 1] = phiS
        if self.carbOn:
            self.phiC[:, : self.stratStep + 1] = phiC

        # Get the total thickness changes induced by compaction and update the elevation accordingly
        dz = totH - np.sum(newH, axis=1)
        dz[dz <= 0] = 0.0
        self.hLocal.setArray(topZ.flatten() - dz.flatten())
        self.dm.localToGlobal(self.hLocal, self.hGlobal)

        # Update each layer thicknesses
        self.stratH[:, : self.stratStep + 1] = newH
        del dz, newH, totH, topZ, phiS, phiF, solidPhase
        del ids, tmpF, tmpS, tot, depth, zlay, cumZ, elev
        if self.carbOn:
            del phiC, tmpC

        gc.collect()

        if MPIrank == 0 and self.verbose:
            print(
                "Compute Lithology Porosity Values (%0.02f seconds)"
                % (process_time() - t0),
                flush=True,
            )

        return
