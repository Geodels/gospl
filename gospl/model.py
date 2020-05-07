import os

from mpi4py import MPI
from time import process_time

if "READTHEDOCS" not in os.environ:
    from .flow import FAMesh as _FAMesh
    from .sed import SEDMesh as _SEDMesh
    from .tools import ReadYaml as _ReadYaml
    from .mesher import UnstMesh as _UnstMesh
    from .tools import WriteMesh as _WriteMesh

    class parentModel(_ReadYaml, _WriteMesh, _UnstMesh, _FAMesh, _SEDMesh):
        def __init__(self):
            pass


else:

    class parentModel(object):
        def __init__(self):
            pass


MPIrank = MPI.COMM_WORLD.Get_rank()


class Model(parentModel):
    """
    Instantiates model object and performs surface processes evolution.

    This object contains methods for the following operations:

     - initialisation of gospl mesh based on input file options
     - computation of surface processes
     - cleaning/destruction of PETSC objects

    :arg filename: YAML input file
    :arg verbose: output option for model main functions
    :arg showlog: Output option for PETSC logging file
    """

    def __init__(self, filename, verbose=True, showlog=False, *args, **kwargs):

        self.showlog = showlog

        self.modelRunTime = process_time()
        self.verbose = verbose

        # Read input dataset
        _ReadYaml.__init__(self, filename)

        # Initialise output mesh
        _WriteMesh.__init__(self)

        # Define unstructured mesh
        _UnstMesh.__init__(self)

        # River flow initialisation
        _FAMesh.__init__(self, *args, **kwargs)

        # Sediment initialisation
        _SEDMesh.__init__(self, *args, **kwargs)

        # Check if simulations just restarted
        if self.rStep > 0:
            _WriteMesh.readData(self)

        # Get external forces
        _UnstMesh.initExtForce(self)

        if MPIrank == 0:
            print(
                "--- Initialisation Phase (%0.02f seconds)"
                % (process_time() - self.modelRunTime),
                flush=True,
            )

        return

    def runProcesses(self):
        """
        Run simulation.

        This function contains methods for the following operations:

         - get flow accumulation
         - apply surface evolution
        """

        self.newForcing = True
        self.steppaleo = 0

        while self.tNow <= self.tEnd:
            tstep = process_time()

            if not self.fast:
                # Compute Flow Accumulation
                _FAMesh.flowAccumulation(self, filled=False)

            # Output time step for first step
            if self.tNow == self.tStart:
                _WriteMesh.visModel(self)

            if not self.fast:
                # Perform River Incision
                _FAMesh.riverIncision(self)
                # Find Continental Sediment Fluxes
                _SEDMesh.getSedFlux(self)
                # Compute Filled Elevation Flow Accumulation
                _FAMesh.flowAccumulation(self, filled=True)
                # Continental and Marine Deposition and Sedimentation
                _SEDMesh.sedChange(self)

            # Update Tectonic, Sea-level & Climatic conditions
            if self.backward and self.tNow < self.tEnd:
                _UnstMesh.applyTectonics(self)

            # Output stratal evolution
            if self.tNow >= self.saveStrat:
                self.stratStep += 1
                self.saveStrat += self.strat

            # Output time step
            _WriteMesh.visModel(self)
            if self.newForcing and self.paleodata is not None:
                _UnstMesh.updatePaleomap(self)

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
                    % (process_time() - tstep),
                    flush=True,
                )

        return

    def reInitialise(self):
        """
        Reinitialise model for paleo-fitting experiments.
        """

        _UnstMesh.reInitialiseModel(self)

        return

    def destroy(self):
        """
        Destroy PETSc DMPlex objects and associated Petsc local/global
        Vectors and Matrices.

        Safely quit model.
        """

        _UnstMesh.destroy_DMPlex(self)

        return
