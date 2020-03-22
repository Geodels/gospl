import os
import gc
import sys
import glob
import h5py
import shutil
import numpy as np

from time import clock
from mpi4py import MPI

from .flow import SPMesh as _SPMesh
from .tools import ReadYaml as _ReadYaml
from .mesher import UnstMesh as _UnstMesh
from .tools import WriteMesh as _WriteMesh

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

    class LandscapeEvolutionModelClass(_ReadYaml, _WriteMesh, _UnstMesh, _SPMesh):
        def __init__(self, filename, verbose=True, showlog=False, *args, **kwargs):

            self.showlog = showlog

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
            _UnstMesh.initExForce(self)

            # Surface processes initialisation
            _SPMesh.__init__(self, *args, **kwargs)

            if MPIrank == 0:
                print(
                    "--- Initialisation Phase (%0.02f seconds)"
                    % (clock() - self.modelRunTime),
                    flush=True,
                )

            return

        def runProcesses(self):
            """
            Run simulation.

            This function contains methods for the following operations:
             - calculating flow accumulation
             - apply surface evolution
            """

            self.newForcing = True
            self.steppaleo = 0

            while self.tNow <= self.tEnd:
                tstep = clock()

                if not self.fast:
                    # Compute Flow Accumulation
                    _SPMesh.FlowAccumulation(self)

                # Output time step for first step
                if self.tNow == self.tStart:
                    _WriteMesh.visModel(self)

                if not self.fast:
                    _SPMesh.sedChange(self)

                # Output stratal evolution
                if self.tNow >= self.saveStrat:
                    _WriteMesh.outputStrat(self)
                    self.saveStrat += self.strat

                # Update Tectonic, Sea-level & Climatic conditions
                if self.backward and self.tNow < self.tEnd:
                    _UnstMesh.applyTectonics(self)

                # Output time step
                _WriteMesh.visModel(self)

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
                        % (clock() - tstep),
                        flush=True,
                    )

            return

        def reInitialise(self):
            """
            Reinitialise model for paleo-fitting experiments.
            """

            _UnstMesh.reInitialiseModel(self)

            return

        def matchPaleo(self):

            _WriteMesh.forceFit(self)

            return

        def destroy(self):
            """
            Destroy PETSc DMPlex objects and associated Petsc local/global
            Vectors and Matrices.

            Safely quit model.
            """

            _UnstMesh.destroy_DMPlex(self)

            return

    return LandscapeEvolutionModelClass(filename, *args, **kwargs)
