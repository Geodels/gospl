import os
import gc
import glob
import h5py
import shutil
import numpy as np
import petsc4py, sys

petsc4py.init(sys.argv)
from time import clock
from mpi4py import MPI


from .mesher import UnstMesh as _UnstMesh
from .tools import ReadYaml as _ReadYaml
from .tools import WriteMesh as _WriteMesh
from .flow import SPMesh as _SPMesh
from .fit import PFit as _PFit

from petsc4py import PETSc as _PETSc

MPIrank = MPI.COMM_WORLD.Get_rank()


def LandscapeEvolutionModel(filename, *args, **kwargs):
    """
    Instantiates eSCAPE model object and performs surface processes evolution.

    This object contains methods for the following operations:
     - initialisation of eSCAPE mesh based on input file options.
     - computation of surface processes
     - cleaning/destruction of PETSC objects

    Args
        filename : YAML input file
        verbose : True/False
            Output option for model main functions
        showlog : True/False
            Output option for PETSC logging file

    Returns:
        LandscapeEvolutionModel : object
    """

    class LandscapeEvolutionModelClass(
        _ReadYaml, _WriteMesh, _UnstMesh, _SPMesh, _PFit
    ):
        def __init__(self, filename, verbose=True, showlog=False, *args, **kwargs):

            self.showlog = showlog
            if self.showlog:
                self.log = _PETSc.Log()
                self.log.begin()

            self.modelRunTime = clock()
            self.verbose = verbose

            # Read input dataset
            _ReadYaml.__init__(self, filename)

            # Initialise output mesh
            _WriteMesh.__init__(self)

            # Define unstructured mesh
            _UnstMesh.__init__(self)

            # Check if simulations just restarted
            if self.rStep > 0:
                _WriteMesh.readData(self)

            # Get external forces
            _UnstMesh.applyForces(self)
            if self.backward:
                _UnstMesh.applyTectonics(self)

            # Surface processes initialisation
            _SPMesh.__init__(self, *args, **kwargs)

            # Paleotopography convergence initialisation
            if self.paleostep > 0:
                _PFit.__init__(self)

            if MPIrank == 0:
                print(
                    "--- Initialisation Phase (%0.02f seconds)"
                    % (clock() - self.modelRunTime)
                )

            return

        def runProcesses(self):
            """
            Run eSCAPE Earth surface processes.

            This function contains methods for the following operations:
             - calculating flow accumulation
             - erosion/deposition induced by stream power law
             - depression identification and pit filling
             - stream induced deposition diffusion
             - hillslope diffusion
            """

            self.newForcing = True
            steppaleo = 0

            while self.tNow <= self.tEnd:
                tstep = clock()

                if not self.fast:
                    # Compute Flow Accumulation
                    _SPMesh.FlowAccumulation(self)

                # Output time step for first step
                if self.tNow == self.tStart:
                    _WriteMesh.outputMesh(self)
                    self.saveTime += self.tout

                if not self.fast:
                    # Compute Erosion using Stream Power Law
                    _SPMesh.cptErosion(self)

                    # Compute Deposition and Sediment Flux
                    _SPMesh.cptSedFlux(self)

                    # Compute Sediment Deposition
                    _SPMesh.SedimentDeposition(self)

                    # Compute Fresh Sediment Diffusion
                    _SPMesh.SedimentDiffusion(self)

                    # Compute Hillslope Diffusion Law
                    _SPMesh.HillSlope(self)

                # Output stratal evolution
                if self.tNow >= self.saveStrat:
                    _WriteMesh.outputStrat(self)
                    self.saveStrat += self.strat

                # Update Tectonic, Sea-level & Climatic conditions
                if self.backward and self.tNow < self.tEnd:
                    _UnstMesh.applyTectonics(self)

                # Output time step
                if self.tNow >= self.saveTime:
                    _WriteMesh.outputMesh(self)
                    self.saveTime += self.tout

                    # Forcing with backward model
                    if self.forceStep >= 0 and self.newForcing:
                        print("Forcing ", steppaleo)
                        _WriteMesh.forcePaleo(self)
                        steppaleo += 1

                    if steppaleo == 1:
                        steppaleo = 0
                        self.newForcing = False
                        self.forceStep += 1
                    else:
                        self.newForcing = True

                # Update Tectonic, Sea-level & Climatic conditions
                if self.tNow < self.tEnd:
                    _UnstMesh.applyForces(self)
                    if not self.backward:
                        _UnstMesh.applyTectonics(self)

                # Advance time
                self.tNow += self.dt

                if MPIrank == 0:
                    print(
                        "--- Computational Step \
                          (%0.02f seconds)"
                        % (clock() - tstep)
                    )

            return

        def reInitialise(self):
            """
            Reinitialise eSCAPE model for paleo-fitting experiments.
            """

            t0step = clock()

            # Restart time
            self.tNow = self.tStart
            self.step = 0
            self.stratStep = 0
            self.rStart = self.tStart
            self.saveTime = self.tNow
            if self.strat > 0:
                self.saveStrat = self.tNow + self.strat
            else:
                self.saveStrat = self.tEnd + self.tout

            # Forcing functions
            self.rainNb = -1
            self.tecNb = -1
            self.flexNb = -1

            # Getting PETSc vectors values
            loadData = np.load(self.meshFile)
            gZ = loadData["z"]
            self.hLocal.setArray(gZ[self.glIDs])
            self.dm.localToGlobal(self.hLocal, self.hGlobal)
            self.vSed.set(0.0)
            self.vSedLocal.set(0.0)
            self.cumED.set(0.0)
            self.cumEDLocal.set(0.0)

            # Update external forces
            _UnstMesh.applyForces(self)
            _UnstMesh.applyTectonics(self)

            del gZ, loadData
            gc.collect()

            if MPIrank == 0:
                print(
                    "--- Reinitialise Phase \
                      (%0.02f seconds)\n+++"
                    % (clock() - t0step)
                )

            return

        def matchPaleo(self):

            self.tNow = self.tEnd
            self.saveStrat = self.tEnd

            # Output stratal evolution
            if self.strat > 0:
                self.stratStep -= 1
                _WriteMesh.outputStrat(self)

            # Output time step
            self.step -= 1
            self.saveTime = self.tEnd
            _WriteMesh.outputMesh(self)

            return

        def destroy(self):
            """
            Destroy PETSc DMPlex objects and associated Petsc local/global
            Vectors and Matrices.

            Safely quit eSCAPE model.
            """

            _UnstMesh.destroy_DMPlex(self)

            if self.showlog:
                self.log.view()

            if MPIrank == 0:
                print(
                    "\n+++\n+++ Total run time (%0.02f seconds)\n+++"
                    % (clock() - self.modelRunTime)
                )

            return

    return LandscapeEvolutionModelClass(filename, *args, **kwargs)
