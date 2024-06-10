import os
import gc
import sys
import petsc4py
import numpy as np
import numpy_indexed as npi

from mpi4py import MPI
from time import process_time

if "READTHEDOCS" not in os.environ:
    from gospl._fortran import donorslist
    from gospl._fortran import donorsmax
    from gospl._fortran import mfdreceivers

petsc4py.init(sys.argv)
MPIrank = petsc4py.PETSc.COMM_WORLD.Get_rank()
MPIcomm = petsc4py.PETSc.COMM_WORLD


class FAMesh(object):
    """
    This class calculates **drainage area** in an implicit, iterative manner using PETSc solvers. It accounts  for multiple flow direction paths (SFD to MFD) based on user input declaration.

    .. note::

        The class follows the parallel approach described in `Richardson et al., 2014 <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2013WR014326>`_ where the iterative nature of the computational algorithms used to solve the linear system creates the possibility of accelerating the solution by providing an initial guess.

    For drainage computation, the class uses a depression-less surface and computes river incision expressed using a **stream power formulation** function of river discharge and slope.

    """

    def __init__(self, *args, **kwargs):
        """
        The initialisation of `FAMesh` class consists in the declaration of PETSc vectors and matrices.
        """

        # KSP solver parameters
        self.rtol = 1.0e-10

        # Identity matrix construction
        self.II = np.arange(0, self.lpoints + 1, dtype=petsc4py.PETSc.IntType)
        self.JJ = np.arange(0, self.lpoints, dtype=petsc4py.PETSc.IntType)
        self.iMat = self._matrix_build_diag(np.ones(self.lpoints))

        # Petsc vectors
        self.fillFAL = self.hLocal.duplicate()
        self.FAG = self.hGlobal.duplicate()
        self.FAL = self.hLocal.duplicate()
        self.hOld = self.hGlobal.duplicate()
        self.hOldLocal = self.hLocal.duplicate()
        self.hOldFlex = self.hLocal.duplicate()
        self.Eb = self.hGlobal.duplicate()
        self.stepED = self.hGlobal.duplicate()
        self.EbLocal = self.hLocal.duplicate()
        self.EbLocal.set(0.0)

        if self.iceOn:
            self.iceFAG = self.hGlobal.duplicate()
            self.iceFAL = self.hLocal.duplicate()

        return

    def _matrix_build(self, nnz=(1, 1)):
        """
        Creates a sparse PETSc matrix.

        .. note::

            To achieve good performance during matrix assembly, the function preallocates the matrix storage by setting the array nnz.

        :arg nnz: array containing the number of nonzeros in the various rows

        :return: sparse PETSc matrix
        """

        matrix = petsc4py.PETSc.Mat().create(comm=MPIcomm)
        matrix.setType("aij")
        matrix.setSizes(self.sizes)
        matrix.setLGMap(self.lgmap_row, self.lgmap_col)
        matrix.setFromOptions()
        matrix.setPreallocationNNZ(nnz)

        return matrix

    def _matrix_build_diag(self, V, nnz=(1, 1)):
        """
        Builds a PETSc diagonal matrix based on a given array `V`

        :arg V: diagonal data array
        :arg nnz: array containing the number of nonzero blocks

        :return: sparse PETSc matrix
        """

        matrix = self._matrix_build()

        # Define diagonal matrix
        matrix.assemblyBegin()
        matrix.setValuesLocalCSR(
            self.II,
            self.JJ,
            V,
        )
        matrix.assemblyEnd()

        return matrix

    def _make_reasons(self, reasons):
        """
        Provides reasons for PETSc error...
        """

        return dict(
            [(getattr(reasons, r), r) for r in dir(reasons) if not r.startswith("_")]
        )

    def _solve_KSP2(self, matrix, vector1, vector2):
        ksp = petsc4py.PETSc.KSP().create(petsc4py.PETSc.COMM_WORLD)
        ksp.setInitialGuessNonzero(True)
        ksp.setOperators(matrix, matrix)
        ksp.setType("fgmres")
        pc = ksp.getPC()
        pc.setType("asm")
        ksp.setTolerances(rtol=1.0e-6, divtol=1.e20)
        ksp.solve(vector1, vector2)
        r = ksp.getConvergedReason()
        if r < 0:
            KSPReasons = self._make_reasons(petsc4py.PETSc.KSP.ConvergedReason())
            if MPIrank == 0:
                print(
                    "LinearSolver failed to converge after iterations",
                    ksp.getIterationNumber(),
                    flush=True,
                )
                print("with reason: ", KSPReasons[r], flush=True)
            vector2.set(0.0)
            ksp.destroy()
            # raise RuntimeError("LinearSolver failed to converge!")
        else:
            ksp.destroy()

        return vector2

    def _solve_KSP(self, guess, matrix, vector1, vector2):
        """

        PETSc *scalable linear equations solvers* (**KSP**) component provides Krylov subspace iterative method and a preconditioner. Here, flow accumulation solution is obtained using PETSc Richardson solver (`richardson`) with block Jacobian preconditioning (`bjacobi`).

        .. note::

            The solver choice was made based on the convergence results from `Richardson et al. (2014) <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2013WR014326>`_ but can be changed if better solver and preconditioner combinations are found.

        Using such iterative method allows for an initial guess to be provided. When this initial guess is close to the solution, the number of iterations required for convergence dramatically decreases. Here the flow discharge solution from previous time step can be passed as an initial `guess` to the solver as discharge often exhibits little change between successive time intervals.

        :arg guess: Boolean specifying if the iterative KSP solver initial guess is nonzero (when provided it corresponds to the previous flow discharge values)
        :arg matrix: PETSc sparse matrix used by the KSP solver composed of diagonal terms set to unity (identity matrix) and off-diagonal terms (weights between 0 and 1). The weights are calculated based on the number of downslope neighbours (based on the chosen number of flow direction directions) and are proportional to the slope.
        :arg vector1: PETSc vector corresponding to the local volume of water available for runoff during a given time step (*e.g.* voronoi area times local precipitation rate)
        :arg vector2: PETSc vector corresponding to the unknown flow discharge values

        :return: vector2 PETSc vector of the new flow discharge values
        """

        ksp = petsc4py.PETSc.KSP().create(petsc4py.PETSc.COMM_WORLD)
        if guess:
            ksp.setInitialGuessNonzero(guess)
        ksp.setOperators(matrix, matrix)
        ksp.setType("richardson")
        pc = ksp.getPC()
        pc.setType("bjacobi")
        ksp.setTolerances(rtol=self.rtol)
        ksp.solve(vector1, vector2)
        r = ksp.getConvergedReason()
        if r < 0:
            ksp.destroy()
            vector2 = self._solve_KSP2(matrix, vector1, vector2)
        else:
            ksp.destroy()

        return vector2
    
    def matrixFlow(self, flowdir, dep=None):
        """
        This function defines the flow direction matrices.

        .. note::

            The matrix is built incrementally looping through the number of flow direction paths defined by the user. It proceeds by assembling a local Compressed Sparse Row (**CSR**) matrix to a global PETSc matrix.

            When setting up the flow matrix in PETSc, we preallocate the non-zero entries of the matrix before starting filling in the values. Using PETSc sparse matrix storage scheme has the advantage that matrix-vector multiplication is extremely fast.

        The  matrix coefficients consist of weights (comprised between 0 and 1) calculated based on the number of downslope neighbours and proportional to the slope.

        :arg flow: boolean to compute matrix for either downstream water or sediment transport

        """

        flowMat = self.iMat.copy()
        indptr = np.arange(0, self.lpoints + 1, dtype=petsc4py.PETSc.IntType)
        nodes = indptr[:-1]
        if dep is None:
            wght = self.wghtVal
        else:
            wght = np.multiply(self.wghtVal,dep.reshape((len(dep),1)))
        rcv = self.rcvID
        for k in range(0, flowdir):

            # Flow direction matrix for a specific direction
            tmpMat = self._matrix_build()
            data = -wght[:, k].copy()
            data[rcv[:, k].astype(petsc4py.PETSc.IntType) == nodes] = 0.0
            tmpMat.assemblyBegin()
            tmpMat.setValuesLocalCSR(
                indptr,
                rcv[:, k].astype(petsc4py.PETSc.IntType),
                data,
            )
            tmpMat.assemblyEnd()
            # Add the weights from each direction
            flowMat += tmpMat
            tmpMat.destroy()
            
        if self.memclear:
            del data, indptr, nodes
            gc.collect()

        # Store flow accumulation matrix
        if dep is None:
            self.fMat = flowMat.transpose().copy()
        else:
            self.fDepMat = flowMat.transpose().copy()

        flowMat.destroy()

        return

    def _buildFlowDirection(self, h, down=True):
        """
        This function builds from neighbouring slopes the flow directions. It calls a fortran subroutine that locally computes for each vertice:

        - the indices of receivers (downstream) nodes depending on the desired number of flow directions (SFD to MFD).
        - the distances to the receivers based on mesh resolution.
        - the associated weights calculated based on the number of receivers and proportional to the slope.

        :arg h: elevation numpy array
        """
        
        # Get open marine regions
        self.seaID = np.where(self.lFill <= self.sealevel)[0]

        # Define multiple flow directions for unfilled elevation
        self.donRcvs, self.distRcv, self.wghtVal = mfdreceivers(
            self.flowDir, self.flowExp, h, self.sealevel
        )
        
        self.rcvID = self.donRcvs.copy()
        self.rcvID[self.ghostIDs,:] = -1
        self.distRcv[self.ghostIDs,:] = 0
        self.wghtVal[self.ghostIDs,:] = 0
              
        if down:
            sum_weight = np.sum(self.wghtVal, axis=1)
            ids = (
                (h == self.lFill)
                & (self.pitIDs > -1)
                & (self.flatDirs > -1)
                & (sum_weight == 0.0)
            )
            ids = ids.nonzero()[0]
            self.rcvID[ids, :] = np.tile(ids, (self.flowDir, 1)).T
            self.rcvID[ids, 0] = self.flatDirs[ids]
            self.wghtVal[ids, :] = 0.0
            self.wghtVal[ids, 0] = 1.0

        # Set borders nodes
        if self.flatModel:
            self.rcvID[self.idBorders, :] = np.tile(self.idBorders, (self.flowDir, 1)).T
            self.distRcv[self.idBorders, :] = 0.0
            self.wghtVal[self.idBorders, :] = 0.0

        # Get local nodes with no receivers as boolean array
        sum_weight = np.sum(self.wghtVal, axis=1)
        lsink = sum_weight == 0.0

        # We don't consider open sea nodes and borders as sinks
        lsink[self.idBorders] = False
        lsink[self.seaID] = False
        lsink = lsink.astype(int) * self.inIDs

        self.lsink = lsink == 1

        self.matrixFlow(self.flowDir)

        return

    def _distributeDownstream(self, pitVol, FA, hl, step, ice=False):
        """
        In cases where rivers flow in depressions, they might fill the sink completely and overspill or remain within the depression, forming a lake. This function computes the excess of water (if any) able to flow dowstream.

        .. important::

            The excess water is then added to the downstream flow accumulation (`FA`) and used to estimate rivers' erosion.

        :arg pitVol: volume of depressions
        :arg FA: excess flow accumulation array
        :arg hl: current elevation array
        :arg step: downstream distribution step

        :return: pitVol, excess (updated volume in each depression and boolean set to True is excess flow remains to be distributed)
        """

        excess = False

        # Remove points belonging to other processors
        FA = np.multiply(FA, self.inIDs)

        # Get volume incoming in each depression
        grp = npi.group_by(self.pitIDs[self.lsink])
        uID = grp.unique
        _, vol = grp.sum(FA[self.lsink])
        inV = np.zeros(len(self.pitParams), dtype=np.float64)
        ids = uID > -1
        inV[uID[ids]] = vol[ids]

        # Combine incoming volume globally
        MPI.COMM_WORLD.Allreduce(MPI.IN_PLACE, inV, op=MPI.SUM)
        
        # Get excess volume to distribute downstream
        eV = inV - pitVol
        if (eV > 0.0).any():
            eIDs = eV > 0.0
            pitVol[eIDs] = 0.0
            spillIDs = self.pitInfo[eIDs, 0]
            localSpill = np.where(self.pitInfo[eIDs, 1] == MPIrank)[0]
            localPts = spillIDs[localSpill]
            nFA = np.zeros(self.lpoints, dtype=np.float64)
            nFA[localPts] = eV[eIDs][localSpill]
            ids = np.in1d(self.pitIDs, np.where(eV > 0.0)[0])
            self.waterFilled[ids] = self.lFill[ids]

        # Update unfilled depressions volumes and assign water level in depressions
        if (eV < 0.0).any():
            eIDs = np.where(eV < 0.0)[0]
            pitVol[eIDs] += eV[eIDs]
            nid = np.absolute(self.filled_vol[eIDs] - pitVol[eIDs][:, None]).argmin(
                axis=1
            )
            fill_lvl = self.filled_lvl[eIDs, nid]
            for k in range(len(eIDs)):
                ids = (self.waterFilled <= fill_lvl[k]) & (self.pitIDs == eIDs[k])
                self.waterFilled[ids] = fill_lvl[k]

        # In case there is still remaining water flux to distribute downstream
        if (eV > 1.0e-3).any():
            if step == 100:
                self._buildFlowDirection(self.lFill)
            else:
                self._buildFlowDirection(self.waterFilled)
            self.tmpL.setArray(nFA / self.dt)
            self.dm.localToGlobal(self.tmpL, self.tmp)
            if self.tmp.sum() > self.maxarea[0]:
                excess = True
                self._solve_KSP(True, self.fMat, self.tmp, self.tmp1)
                self.dm.globalToLocal(self.tmp1, self.tmpL)
                nFA = self.tmpL.getArray().copy() * self.dt
                FA = nFA.copy()
                FA[hl < self.waterFilled] = 0.0
                self.tmpL.setArray(FA / self.dt)
                if ice:
                    self.iceFAL.axpy(1.0, self.tmpL)
                else:
                    self.FAL.axpy(1.0, self.tmpL)
            else:
                nFA = None
        else:
            nFA = None

        return excess, pitVol, nFA

    def _iceFlow(self):


        ti = process_time()
        
        # Get depressions information
        ti1 = process_time()
        pitVol = self.pitParams[:, 0].copy()
        hl = self.hLocal.getArray().copy()

        # Solve ice flow accumulation
        iceA = self.bL.getArray().copy()
        iceA[self.seaID] = 0.
        tmp = (hl - self.elaH) / (self.iceH - self.elaH)
        self.iceIDs = tmp > 0
        tmp[tmp>1.] = 1. 
        tmp[tmp<0.] = 0.
        iceA = np.multiply(iceA, tmp)
        self.tmpL.setArray(iceA)
        self.dm.localToGlobal(self.tmpL, self.tmp)
        self._solve_KSP(True, self.fMat, self.tmp, self.iceFAG)
        self.dm.globalToLocal(self.iceFAG, self.iceFAL)

        # # Volume of ice flowing downstream
        self.waterFilled = hl.copy()
        if (pitVol > 0.0).any():
            iFA = self.iceFAL.getArray().copy() * self.dt
            _, _, iFA = self._distributeDownstream(pitVol, iFA, hl, 100, ice=True)
        if MPIrank == 0 and self.verbose:
            print(
                "Downstream ice flow computation (%0.02f seconds)" % (process_time() - ti1),
                flush=True,
            )

        # Smooth the ice flow across cells
        ti1 = process_time()
        self.iceFAL.copy(result=self.tmpL)
        tmp = self.tmpL.getArray().copy()
        self.dm.localToGlobal(self.tmpL, self.tmp1)
        smthIce = self._hillSlope(smooth=2)
        self.tmpL.setArray(smthIce*self.scaleIce)
        self.dm.localToGlobal(self.tmpL, self.tmp1)
        smthIce[~self.iceIDs] = tmp[~self.iceIDs]
        self.iceFAL.setArray(smthIce)
        self.dm.localToGlobal(self.iceFAL, self.iceFAG)
        if MPIrank == 0 and self.verbose:
            print(
                "Glaciers Accumulation (%0.02f seconds)" % (process_time() - ti1),
                flush=True,
            )

        # Update fluvial flow accumulation
        PA = self.FAL.getArray().copy()
        PA[self.iceIDs] -= smthIce[self.iceIDs]
        PA[~self.iceIDs] += smthIce[~self.iceIDs]
        PA[PA<0] = 0.
        self.FAL.setArray(PA)
        self.dm.localToGlobal(self.FAL, self.FAG)

        PA = self.fillFAL.getArray().copy()
        PA[self.iceIDs] -= smthIce[self.iceIDs]
        PA[~self.iceIDs] += smthIce[~self.iceIDs]
        PA[PA<0] = 0.
        self.fillFAL.setArray(PA)

        if MPIrank == 0 and self.verbose:
            print(
                "Compute IceFlow Accumulation (%0.02f seconds)" % (process_time() - ti),
                flush=True,
            )
        
        return

    def flowAccumulation(self):
        """
        This function is the **main entry point** for flow accumulation computation.

        .. note::

            Flow accumulation (`FA`) calculations are a core component of landscape evolution models as they are often used as proxy to estimate flow discharge, sediment load, river width, bedrock erosion, and sediment deposition. Until recently, conventional `FA` algorithms were serial and limited to small spatial problems.

        `gospl` model computes the flow discharge from `FA` and the net precipitation rate using a **parallel implicit drainage area (IDA) method** proposed by `Richardson et al., 2014 <https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2013WR014326>`_ but adapted to unstructured grids.

        It calls the following *private functions*:

        1. _buildFlowDirection
        2. _solve_KSP
        3. _distributeDownstream

        """

        t0 = process_time()

        # Compute depressions information
        self.fillElevation(sed=False)
        pitVol = self.pitParams[:, 0].copy()

        # Build flow direction and downstream matrix
        hl = self.hLocal.getArray().copy()

        self._buildFlowDirection(hl, False)
        self.wghtVali = self.wghtVal.copy()
        self.rcvIDi = self.rcvID.copy()
        self.distRcvi = self.distRcv.copy()
        self.fMati = self.fMat.copy()
        self.lsinki = self.lsink.copy()

        # Solve flow accumulation
        rainA = self.bL.getArray().copy()
        rainA[self.seaID] = 0.
        self.bL.setArray(rainA)
        self.dm.localToGlobal(self.bL, self.bG)
        self._solve_KSP(True, self.fMat, self.bG, self.FAG)
        self.dm.globalToLocal(self.FAG, self.FAL)

        # Volume of water flowing downstream
        self.waterFilled = hl.copy()
        if (pitVol > 0.0).any():
            FA = self.FAL.getArray().copy() * self.dt
            excess = True
            step = 0
            while excess:
                t1 = process_time()
                excess, pitVol, FA = self._distributeDownstream(pitVol, FA, hl, step)
                if MPIrank == 0 and self.verbose:
                    print(
                        "Downstream flow computation step %d (%0.02f seconds)"
                        % (step, process_time() - t1),
                        flush=True,
                    )
                step += 1

            # Get overall water flowing donwstream accounting for filled depressions
            FA = self.FAL.getArray().copy()
            ids = (hl <= self.waterFilled) & (self.pitIDs > -1)
            FA[ids] = 0.0
            self.FAL.setArray(FA)
            self.dm.localToGlobal(self.FAL, self.FAG)
            self.dm.globalToLocal(self.FAG, self.FAL)
            FA[ids] = self.FAG.max()[1] * 0.1
            self.fillFAL.setArray(FA)
        else:
            self.FAL.copy(result=self.fillFAL)

        # Get water level
        self.waterFilled -= hl

        if MPIrank == 0 and self.verbose:
            print(
                "Compute Flow Accumulation (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        # Compute glacier flows
        if self.iceOn:
            self._iceFlow()

        return

    def _upstreamDeposition(self, flowdir):
        """
        Check if over deposition occurs in depressions and/or flat regions and enforce a maximum deposition rate based on upstream elevations computed from the implicit erosion/deposition equation.

        """
        niter = 0
        tolerance = 1e-3 
        equal = False
        h = self.hLocal.getArray().copy()

        # Get donors list based on elevation
        donors = donorslist(flowdir, self.inIDs, self.donRcvs)
        topIDs = np.where(np.max(donors, axis=1)==-1)[0]
        while not equal and niter < 10000:

            Eb = -self.EbLocal.getArray().copy()*self.dt
            elev = (Eb + h).copy()
            elev[topIDs] = h[topIDs]
            maxh = donorsmax(elev,donors)
            maxh[maxh==-1.e8] = elev[maxh==-1.e8]
            self.tmpL.setArray(maxh)
            self.dm.localToGlobal(self.tmpL, self.tmp)
            self.dm.globalToLocal(self.tmp, self.tmpL)
            maxh = self.tmpL.getArray().copy()-1.e-6
            maxh[maxh>elev] = elev[maxh>elev]
            
            self.tmpL.setArray(maxh-elev)
            self.dm.localToGlobal(self.tmpL, self.dh)
            self.dh.abs()
            maxdh = self.dh.max()[1]

            equal = maxdh < tolerance
            self.EbLocal.setArray(-(maxh - h)/self.dt)
            self.dm.localToGlobal(self.EbLocal, self.Eb)
            self.dm.globalToLocal(self.Eb, self.EbLocal)
            niter += 1

        if self.memclear:
            del donors, h, elev, maxh, Eb
            gc.collect()

        return

    def _getEroDepRate(self):
        r"""
        This function computes erosion deposition rates in metres per year. This is done on the filled elevation. We use the filled-limited elevation to ensure that erosion/deposition is not going to be underestimated by small depressions which are likely to be filled (either by sediments or water) during a single time step.

        The simplest law to simulate fluvial incision is based on the detachment-limited stream power law, in which erosion rate  depends on drainage area :math:`A`, net precipitation :math:`P` and local slope :math:`S` and takes the form:

        .. math::

          E = − \kappa P^d (PA)^m S^n

        :math:`\kappa` is a dimensional coefficient describing the erodibility of the channel bed as a function of rock strength, bed roughness and climate, :math:`d`, :math:`m` and :math:`n` are dimensionless positive constants.

        Default formulation assumes :math:`d = 0`, :math:`m = 0.5` and :math:`n = 1`. The precipitation exponent :math:`d` allows for representation of climate-dependent chemical weathering of river bed across non-uniform rainfall.

        .. important::

            In `gospl`, the coefficient `n` is fixed and the only variables that the user can tune are the coefficients `m`, `d` and the erodibility :math:`\kappa`.

        The erosion rate is solved by an implicit time integration method, the matrix system is based on the receiver distributions and is assembled from local Compressed Sparse Row (**CSR**) matrices into a global PETSc matrix. The PETSc *scalable linear equations solvers* (**KSP**) is used with both an iterative method and a preconditioner and erosion rate solution is obtained using PETSc Richardson solver (`richardson`) with block Jacobian preconditioning (`bjacobi`).

        Once the erosion rate solution has been obtained, local sediment flux depends on upstream fluxes, local eroded flux and local deposition flux. We assume that local deposition depends on the user-defined forced deposition value (:math:`fDep`), the local water flux and the cell area. Here again the sediment flux is determined implicitly and corresponding deposition flux are calculated subsequently once the local total flux is known.      
        """

        t0 = process_time()
        hOldArray = self.hLocal.getArray().copy()
        self.oldH = hOldArray.copy()
        hOldArray[self.seaID] = self.sealevel
        if self.flexOn:
            self.hLocal.copy(result=self.hOldFlex)

        # Upstream-averaged mean annual precipitation rate based on drainage area
        PA = self.FAL.getArray()

        # Incorporate the effect of local mean annual precipitation rate on erodibility
        if self.sedfacVal is not None:
            Kbr = self.K * self.sedfacVal * (self.rainVal ** self.coeffd)
        else:
            Kbr = self.K * (self.rainVal ** self.coeffd)
        Kbr *= self.dt * (PA ** self.spl_m)
        Kbr[self.seaID] = 0.0

        if self.iceOn:
            GA = self.iceFAL.getArray()
            GA[~self.iceIDs] = 0.
            Kbi = self.Kice * GA

        # Dimensionless depositional coefficient
        if self.fDepa > 0:
            fDep = np.divide(self.fDepa*self.larea, PA, out=np.zeros_like(PA), where=PA != 0)
            if self.dmthd == 1:
                fDep[fDep>0.99] = 0.99
                self.matrixFlow(self.flowDir, 1.-fDep)
            else:
                dMat = self._matrix_build_diag(fDep)
                dMat += self.fMat
    
        # Initialise matrices...
        eMat = self.iMat.copy()
        wght = self.wghtVali.copy()
        indptr = np.arange(0, self.lpoints + 1, dtype=petsc4py.PETSc.IntType)
        nodes = indptr[:-1]
        wght[self.seaID, :] = 0.0

        if self.iceOn:
            gMat = self.iMat.copy()

        # Define erosion coefficients
        for k in range(0, self.flowDir):

            # Define erosion limiter to prevent formation of flat
            dh = hOldArray - hOldArray[self.rcvIDi[:, k]]
            limiter = np.divide(dh, dh + 1.0e-2, out=np.zeros_like(dh), where=dh != 0)

            # Bedrock erosion processes SPL computation (maximum bedrock incision)
            data = np.divide(
                Kbr * limiter,
                self.distRcvi[:, k],
                out=np.zeros_like(PA),
                where=self.distRcvi[:, k] != 0,
            )
            tmpMat = self._matrix_build()
            data = np.multiply(data, -wght[:, k])
            data[self.rcvIDi[:, k].astype(petsc4py.PETSc.IntType) == nodes] = 0.0
            tmpMat.assemblyBegin()
            tmpMat.setValuesLocalCSR(
                indptr,
                self.rcvIDi[:, k].astype(petsc4py.PETSc.IntType),
                data,
            )
            tmpMat.assemblyEnd()
            eMat += tmpMat
            eMat -= self._matrix_build_diag(data)
            tmpMat.destroy()

            if self.iceOn:
                data = np.divide(
                    Kbi * limiter,
                    self.distRcvi[:, k],
                    out=np.zeros_like(PA),
                    where=self.distRcvi[:, k] != 0,
                )
                tmpMat = self._matrix_build()
                data = np.multiply(data, -wght[:, k])
                data[self.rcvIDi[:, k].astype(petsc4py.PETSc.IntType) == nodes] = 0.0
                tmpMat.assemblyBegin()
                tmpMat.setValuesLocalCSR(
                    indptr,
                    self.rcvIDi[:, k].astype(petsc4py.PETSc.IntType),
                    data,
                )
                tmpMat.assemblyEnd()
                gMat += tmpMat
                gMat -= self._matrix_build_diag(data)
                tmpMat.destroy()

        if self.memclear:
            del dh, limiter, wght, data
            gc.collect()

        # Solve SPL erosion implicitly for fluvial erosion
        t1 = process_time()
        self._solve_KSP(True, eMat, self.hOld, self.stepED)
        self.tmp.waxpy(-1.0, self.hOld, self.stepED)
        eMat.destroy()
        if MPIrank == 0 and self.verbose:
            print(
                "Solve SPL erosion (%0.02f seconds)" % (process_time() - t1),
                flush=True,
            )
        # Solve SPL erosion implicitly for glacial erosion
        if self.iceOn:
            t1 = process_time()
            self._solve_KSP(False, gMat, self.hOld, self.stepED)
            self.tmp1.waxpy(-1.0, self.hOld, self.stepED)
            gMat.destroy()
            self.tmp.axpy(1.,self.tmp1)
            if MPIrank == 0 and self.verbose:
                print(
                    "Solve glacial erosion (%0.02f seconds)" % (process_time() - t1),
                    flush=True,
                )

        # Implicit sediment fluxes combining upstream flux, erosion and deposition
        if self.fDepa > 0:
            t1 = process_time()
            self.dm.globalToLocal(self.tmp, self.tmpL)
            QsL = -self.tmpL.getArray()*self.larea
            QsL = np.divide(QsL, self.dt)
            QsL[QsL<1.e-8] = 0.
            self.tmpL.setArray(QsL)
            self.dm.localToGlobal(self.tmpL, self.tmp)
            if self.dmthd == 1:
                self._solve_KSP(False, self.fDepMat, self.tmp, self.tmp1)
                self.fDepMat.destroy()
            else:
                self._solve_KSP(False, dMat, self.tmp, self.tmp1)
                dMat.destroy()
            if MPIrank == 0 and self.verbose:
                print(
                    "Solve sediment fluxes (%0.02f seconds)" % (process_time() - t1),
                    flush=True,
                )
            
            # Extract local sediment deposition fluxes
            self.dm.globalToLocal(self.tmp1, self.tmpL)
            if self.dmthd == 1:
                QsT = self.tmpL.getArray()
                scale = np.divide(fDep, 1.0-fDep, out=np.zeros_like(fDep), where=fDep != 0)
                QsD = (QsT-QsL)*scale
            else:
                QsD = self.tmpL.getArray()*fDep
            self.tmpL.setArray((QsD-QsL)*self.dt/self.larea)
            self.dm.localToGlobal(self.tmpL, self.tmp)

        # Define erosion rate (positive for incision)
        t1 = process_time()
        E = -self.tmp.getArray().copy()
        E = np.divide(E, self.dt)
        self.Eb.setArray(E)
        self.dm.globalToLocal(self.Eb, self.EbLocal)
        E = self.EbLocal.getArray().copy()
        if MPIrank == 0 and self.verbose:
            print(
                "Define erosion rate (%0.02f seconds)" % (process_time() - t1),
                flush=True,
            )
        if self.flatModel:
            E[self.idBorders] = 0.0
        E[self.lsink] = 0.0
        self.EbLocal.setArray(E)
        self.dm.localToGlobal(self.EbLocal, self.Eb)

        # Limiting deposition rates based on upstream conditions
        if self.fDepa > 0:
            self._upstreamDeposition(self.flowDir)

        if self.memclear:
            del E, PA, Kbr, QsL, QsD, fDep
            gc.collect()
        
        if MPIrank == 0 and self.verbose:
            print(
                "Finalise erosion deposition rates (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        return
    
    def erodepSPL(self):
        """
        Modified **stream power law** model used to represent erosion by rivers also taking into account the role played by sediment in modulating erosion and deposition rate.
        
        It calls the private function `_getEroDepRate <https://gospl.readthedocs.io/en/latest/api.html#flow.flowplex.FAMesh._getEroDepRate>`_ described above. Once erosion/deposition rates have been calculated, the function computes local thicknesses for the considered time step and update local elevation and cumulative erosion, deposition values.

        If multiple lithologies are considered, the stratigraphic record is updated based on eroded/deposited thicknesses.

        """

        t0 = process_time()

        # Computes the erosion deposition rates based on flow accumulation
        self.Eb.set(0.0)
        self.hGlobal.copy(result=self.hOld)
        self.dm.globalToLocal(self.hOld, self.hOldLocal)
        self._getEroDepRate()

        # Get erosion / deposition thicknesses
        Eb = self.Eb.getArray().copy()
        self.tmp.setArray(-Eb * self.dt)
        self.cumED.axpy(1.0, self.tmp)
        self.dm.globalToLocal(self.cumED, self.cumEDLocal)
        self.hGlobal.axpy(1.0, self.tmp)
        self.dm.globalToLocal(self.hGlobal, self.hLocal)

        # Update stratigraphic layers
        if self.stratNb > 0:
            self.erodeStrat()
            self.deposeStrat()

        if MPIrank == 0 and self.verbose:
            print(
                "Get Erosion Deposition values (%0.02f seconds)" % (process_time() - t0),
                flush=True,
            )

        if self.memclear:
            del Eb
            gc.collect()

        return