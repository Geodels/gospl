import os

from mpi4py import MPI
from time import process_time

if "READTHEDOCS" not in os.environ:
    from .flow import FAMesh as _FAMesh
    from .sed import SEDMesh as _SEDMesh
    from .sed import STRAMesh as _STRAMesh
    from .tools import ReadYaml as _ReadYaml
    from .mesher import UnstMesh as _UnstMesh
    from .tools import WriteMesh as _WriteMesh

else:

    class _ReadYaml(object):
        def __init__(self, filename):
            print("Fake print statement for readthedocs", filename)

    class _UnstMesh(object):
        def __init__(self):
            pass

    class _WriteMesh(object):
        def __init__(self):
            pass

    class _FAMesh(object):
        def __init__(self):
            pass

    class _SEDMesh(object):
        def __init__(self):
            pass

    class _STRAMesh(object):
        def __init__(self):
            pass


MPIrank = MPI.COMM_WORLD.Get_rank()


class Model(_ReadYaml, _WriteMesh, _UnstMesh, _FAMesh, _SEDMesh, _STRAMesh):
    """
    Instantiates model object and performs surface processes evolution.

    This object contains methods for the following operations:

     - initialisation of gospl mesh based on input file options
     - computation of surface processes over time
     - cleaning/destruction of PETSC objects

    :arg filename: YAML input file
    :arg verbose: output flag for model main functions
    :arg showlog: output flag for PETSC logging file
    :arg carbctrl: carbonate control option

    """

    def __init__(
        self, filename, verbose=True, showlog=False, carbctrl=None, *args, **kwargs
    ):

        self.showlog = showlog

        self.modelRunTime = process_time()
        self.verbose = verbose

        self.carbOn = False
        if carbctrl is not None:
            self.carbOn = True
            self.carbCtrl = carbctrl

        # Read input dataset
        _ReadYaml.__init__(self, filename)

        # Stratigraphy initialisation
        _STRAMesh.__init__(self)

        # Define unstructured mesh
        _UnstMesh.__init__(self)

        # Initialise output mesh
        _WriteMesh.__init__(self)

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
        Runs simulation over time.

        This function contains methods for the following operations:

         - computes flow accumulation based on imposed precipitation field
         - performs land surface evolution from river erosion, transport and deposition
         - executes creep processes and marine depostion (linear hillslope diffusion)
         - records stratigraphic layers evolution and associated porosity variations
         - applies user-defined tectonics forcing (horizontal and vertical displacements)

        """

        self.newForcing = True
        self.steppaleo = 0

        while self.tNow <= self.tEnd:
            tstep = process_time()

            if not self.fast:
                # Compute Flow Accumulation
                _FAMesh.flowAccumulation(self)

            # Output time step for first step
            if self.tNow == self.tStart:
                _WriteMesh.visModel(self)

            if not self.fast:
                # Perform River Incision
                _FAMesh.riverIncision(self)
                # Find Continental Sediment Fluxes
                _SEDMesh.getSedFlux(self)
                # Downstream sediment deposition
                _SEDMesh.sedChange(self)
                # Hillslope diffusion
                _SEDMesh.getHillslope(self)

            # Update Tectonic, Sea-level & Climatic conditions
            if self.backward and self.tNow < self.tEnd:
                _UnstMesh.applyTectonics(self)

            # Create new stratal layer
            if self.tNow >= self.saveStrat:
                # Stratigraphic Layer Porosity and Thicknesses under Compaction
                _STRAMesh.getCompaction(self)
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

    def reInitialiseZ(self):
        """
        Reinitialise model elevation.

        This function clears PETSc vectors and forcing conditions without having to reset the mesh structure.
        """

        _UnstMesh.reInitialiseElev(self)

        return

    def destroy(self):
        """
        Destroy PETSc DMPlex objects and associated Petsc local/global Vectors and Matrices.

        Safely quit model.
        """

        _UnstMesh.destroy_DMPlex(self)

        return
