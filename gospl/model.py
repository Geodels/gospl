import os

from mpi4py import MPI
from time import process_time

if "READTHEDOCS" not in os.environ:
    from .flow import FAMesh as _FAMesh
    from .flow import IceMesh as _IceMesh
    from .flow import PITFill as _PITFill
    from .eroder import SPL as _SPL
    from .eroder import nlSPL as _nlSPL
    from .eroder import soilSPL as _soilSPL
    from .sed import SEDMesh as _SEDMesh
    from .sed import hillSLP as _hillSLP
    from .sed import SEAMesh as _SEAMesh
    from .sed import STRAMesh as _STRAMesh
    from .tools import ReadYaml as _ReadYaml
    from .mesher import UnstMesh as _UnstMesh
    from .mesher import VoroBuild as _VoroBuild
    from .tools import GridProcess as _GridProcess
    from .mesher import Tectonics as _Tectonics
    from .tools import WriteMesh as _WriteMesh

else:

    class _ReadYaml(object):
        def __init__(self, filename):
            print("Fake print statement for readthedocs", filename)

    class _VoroBuild(object):
        def __init__(self):
            pass

    class _UnstMesh(object):
        def __init__(self):
            pass

    class _Tectonics(object):
        def __init__(self):
            pass

    class _WriteMesh(object):
        def __init__(self):
            pass

    class _FAMesh(object):
        def __init__(self):
            pass

    class _IceMesh(object):
        def __init__(self):
            pass

    class _SPL(object):
        def __init__(self):
            pass

    class _nlSPL(object):
        def __init__(self):
            pass

    class _soilSPL(object):
        def __init__(self):
            pass

    class _PITFill(object):
        def __init__(self):
            pass

    class _SEDMesh(object):
        def __init__(self):
            pass

    class _hillSLP(object):
        def __init__(self):
            pass

    class _SEAMesh(object):
        def __init__(self):
            pass

    class _STRAMesh(object):
        def __init__(self):
            pass

    class _GridProcess(object):
        def __init__(self):
            pass

MPIrank = MPI.COMM_WORLD.Get_rank()


class Model(
    _ReadYaml,
    _WriteMesh,
    _UnstMesh,
    _VoroBuild,
    _GridProcess,
    _Tectonics,
    _FAMesh,
    _IceMesh,
    _SPL,
    _nlSPL,
    _soilSPL,
    _PITFill,
    _SEDMesh,
    _hillSLP,
    _SEAMesh,
    _STRAMesh,
):
    """
    Instantiates model's objects and initialise classes.

    This object contains methods for the following operations:

     - initialisation of goSPL mesh based on input file options
     - computation of surface processes over time
     - cleaning/destruction of PETSC objects

    :arg filename: YAML input file
    :arg verbose: output flag for model main functions
    :arg showlog: output flag for PETSC logging file

    """

    def __init__(
        self, filename, verbose=True, showlog=False, *args, **kwargs
    ):

        self.showlog = showlog

        self.modelRunTime = process_time()
        self.verbose = verbose

        # Read input dataset
        _ReadYaml.__init__(self, filename)

        # Stratigraphy initialisation
        _STRAMesh.__init__(self)

        # Define voronoi mesh
        _VoroBuild.__init__(self)

        # Define unstructured mesh
        _UnstMesh.__init__(self)

        # Initialise output mesh
        _WriteMesh.__init__(self)

        # River flow initialisation
        _FAMesh.__init__(self, *args, **kwargs)

        # Ice flow initialisation
        _IceMesh.__init__(self, *args, **kwargs)

        # SPL initialisation
        _SPL.__init__(self, *args, **kwargs)

        # Non-linear SPL initialisation
        _nlSPL.__init__(self, *args, **kwargs)

        # Non-linear SPL with soil generation initialisation
        _soilSPL.__init__(self, *args, **kwargs)

        # Pit filling initialisation
        _PITFill.__init__(self, *args, **kwargs)

        # Continental sediment transfer initialisation
        _SEDMesh.__init__(self, *args, **kwargs)

        # Hillslope processes (linear and non-linear) initialisation
        _hillSLP.__init__(self, *args, **kwargs)

        # Marine sediment transport and deposition initialisation
        _SEAMesh.__init__(self, *args, **kwargs)

        # Define additional grid processes (flexure, orographic rain)
        _GridProcess.__init__(self)

        # Get external forces driving landscape dynamics over time (tectonic, rainfall...)
        _UnstMesh.applyForces(self)

        # Initialise horizontal tectonics forcings (multiple advection techniques)
        _Tectonics.__init__(self)

        # Check if simulation just restarted
        if self.rStep > 0:
            _WriteMesh.readData(self)

        if not self.fast:
            # Compute flow accumulation
            _FAMesh.flowAccumulation(self)

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
         - executes creep processes and marine depostion (linear and non-linear diffusion)
         - records stratigraphic layers evolution and associated porosity variations
         - applies user-defined tectonics forcing (horizontal and vertical displacements)

        """

        while self.tNow <= self.tEnd:
            tstep = process_time()

            # Output time step
            _Tectonics.updatePaleoZ(self)
            _WriteMesh.visModel(self)
            if self.tNow == self.tEnd:
                return

            # Create new stratal layer
            if self.tNow >= self.saveStrat:
                self.stratStep += 1
                self.saveStrat += self.strat

            # Perform advection and tectonics
            _Tectonics.getTectonics(self)

            if not self.fast:
                if self.iceOn:
                    # Compute ice accumulation
                    _IceMesh.iceAccumulation(self)
                # Compute flow accumulation
                _FAMesh.flowAccumulation(self)

                # Perform River Incision/Deposition based on Stream Power Law (different flavors)
                if self.cptSoil:
                    # Non-linear coupled to soil production
                    _soilSPL.erodepSPLsoil(self)
                elif self.spl_n == 1.0:
                    # Linear slope dependencies
                    _SPL.erodepSPL(self)
                else:
                    # Non-linear slope dependencies
                    _nlSPL.erodepSPLnl(self)

                if not self.nodep:
                    # Downstream sediment deposition over the continents
                    _FAMesh.flowAccumulation(self)
                    _SEDMesh.sedChange(self)
                    if self.seaDepo:
                        # Downstream sediment deposition in marine environments
                        _SEAMesh.seaChange(self)

                # Hillslope diffusion (linear and non-linear)
                _hillSLP.getHillslope(self)

            if self.tNow >= self.saveStrat:
                # Stratigraphic layer porosity and thicknesses under compaction
                _STRAMesh.getCompaction(self)

            # Apply flexural isostasy (local and global)
            if self.flexOn:
                _GridProcess.applyFlexure(self)

            # Update tectonic, sea-level & climatic conditions
            if self.tNow < self.tEnd:
                _UnstMesh.applyForces(self)

            # Advance time
            self.tNow += self.dt

            if MPIrank == 0:
                print(
                    "--- Computational Step (%0.02f seconds) | Time Step: %d years"
                    % (process_time() - tstep, self.tNow),
                    flush=True,
                )

        return

    def destroy(self):
        """
        Destroy PETSc DMPlex objects and associated Petsc local/global Vectors and Matrices.

        Safely quit model.
        """

        _UnstMesh.destroy_DMPlex(self)

        return
