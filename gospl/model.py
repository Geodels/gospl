import os
import sys

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
    from .tools import Profiler as _Profiler

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

    class _Profiler(object):
        def __init__(self, *args, **kwargs):
            self.enabled = False

        def phase(self, name):
            from contextlib import nullcontext

            return nullcontext()

        def report(self, *args, **kwargs):
            return None

MPIrank = MPI.COMM_WORLD.Get_rank()
MPIsize = MPI.COMM_WORLD.Get_size()

_mpi_excepthook_installed = False


def _install_mpi_abort_excepthook():
    """
    Install a ``sys.excepthook`` that ``MPI_Abort``s **every** rank on any
    uncaught exception (parallel runs only; a no-op at ``MPIsize == 1``).

    Why: without it, an exception on ONE rank — especially one raised inside a
    collective (a failed KSP solve, a PETSc collective error such as the
    PETSc-3.21 ``garbage_cleanup`` ``MPI_ERR_BUFFER`` bug, …) — leaves the OTHER
    ranks blocked in MPI forever, so the whole job **hangs** and must be killed
    by hand. ``MPI.COMM_WORLD.Abort(1)`` turns that deadlock into an immediate,
    clean termination of all ranks that still surfaces the traceback. Idempotent
    and chains to the previous hook so the traceback prints as usual first.
    """
    global _mpi_excepthook_installed
    if _mpi_excepthook_installed or MPIsize == 1:
        return
    _prev_hook = sys.excepthook

    def _hook(exc_type, exc, tb):
        try:
            _prev_hook(exc_type, exc, tb)
            sys.stderr.flush()
            sys.stdout.flush()
        finally:
            # A lone rank's exception otherwise deadlocks the others waiting in
            # a collective — terminate the entire job.
            MPI.COMM_WORLD.Abort(1)

    sys.excepthook = _hook
    _mpi_excepthook_installed = True


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
        self, filename, verbose=True, showlog=False, profile=False, *args, **kwargs
    ):

        self.showlog = showlog

        # Make any uncaught exception abort ALL ranks instead of deadlocking the
        # job (parallel runs only). Installed first so even an init-time failure
        # terminates cleanly. See _install_mpi_abort_excepthook.
        _install_mpi_abort_excepthook()

        self.modelRunTime = process_time()
        self.verbose = verbose

        # Read input dataset
        _ReadYaml.__init__(self, filename)

        # Wall-clock phase profiler. Enabled by the `profile` kwarg OR the
        # YAML `output: profile: true` flag (parsed in _readOut). When off it
        # is a zero-overhead no-op. See gospl/tools/profiler.py.
        prof_on = profile or bool(getattr(self, "profileFlag", False))
        self.profiler = _Profiler(enabled=prof_on)

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

        runStart = MPI.Wtime()

        while self.tNow <= self.tEnd:
            tstep = process_time()

            # Flexure interval counter (see GridProcess + flex_interval): the
            # load reference is snapshotted at interval starts (in the eroder)
            # and flexure applied at interval ends below.
            if self.flexOn:
                self.flexCount += 1

            # Output time step
            with self.profiler.phase("tectonics"):
                _Tectonics.updatePaleoZ(self)
            with self.profiler.phase("output"):
                _WriteMesh.visModel(self)
            if self.tNow == self.tEnd:
                break

            # Create new stratal layer
            newLayer = self.tNow >= self.saveStrat
            if newLayer:
                self.stratStep += 1
                self.saveStrat += self.strat

            # Perform advection and tectonics
            with self.profiler.phase("tectonics"):
                _Tectonics.getTectonics(self)

            if not self.fast:
                if self.iceOn:
                    # Compute ice accumulation
                    with self.profiler.phase("ice"):
                        _IceMesh.iceAccumulation(self)
                # Compute flow accumulation
                with self.profiler.phase("flow"):
                    _FAMesh.flowAccumulation(self)

                # Perform River Incision/Deposition based on Stream Power Law (different flavors)
                with self.profiler.phase("erosion"):
                    if self.cptSoil:
                        # Non-linear coupled to soil production
                        _soilSPL.erodepSPLsoil(self)
                    elif self.spl_n == 1.0:
                        # Linear slope dependencies
                        _SPL.erodepSPL(self)
                    else:
                        # Non-linear slope dependencies
                        _nlSPL.erodepSPLnl(self)

                # Glacial till: abrasion-produced till transported by ice and
                # deposited (melt-out) as moraine in the ablation zone (glacial
                # till only; no-op otherwise). Done before fluvial deposition
                # so the moraine is reworked by meltwater/rivers this step.
                if self.iceOn:
                    with self.profiler.phase("till"):
                        _IceMesh.glacialTill(self)

                if not self.nodep:
                    # Downstream sediment deposition over the continents
                    with self.profiler.phase("flow"):
                        _FAMesh.flowAccumulation(self)
                    with self.profiler.phase("sed"):
                        _SEDMesh.sedChange(self)
                    if self.seaDepo:
                        # Downstream sediment deposition in marine environments
                        with self.profiler.phase("sea"):
                            _SEAMesh.seaChange(self)

                # Hillslope diffusion (linear and non-linear)
                with self.profiler.phase("hillslope"):
                    _hillSLP.getHillslope(self)

            if newLayer:
                # Stratigraphic layer porosity and thicknesses under compaction
                with self.profiler.phase("strat"):
                    _STRAMesh.getCompaction(self)

            # Apply flexural isostasy (local and global). With flex_interval > 1
            # the load accumulates and flexure is solved only at interval ends
            # (interval = 1 → every step, the default/unchanged behaviour).
            if self.flexOn and (
                self.flexCount % self.flex_interval == self.flex_interval - 1
            ):
                with self.profiler.phase("flexure"):
                    _GridProcess.applyFlexure(self)

            # Update tectonic, sea-level & climatic conditions
            if self.tNow < self.tEnd:
                with self.profiler.phase("forcing"):
                    _UnstMesh.applyForces(self)

            # Advance time
            self.tNow += self.dt

            if MPIrank == 0:
                print(
                    "--- Computational Step (%0.02f seconds) | Time Step: %d years"
                    % (process_time() - tstep, self.tNow),
                    flush=True,
                )

        # Cross-rank wall-clock profile + machine-readable profile.json
        # (collective; a no-op when profiling is disabled).
        self.profiler.report(self.outputDir, total_wall=MPI.Wtime() - runStart)

        return

    def destroy(self):
        """
        Destroy PETSc DMPlex objects and associated Petsc local/global Vectors and Matrices.

        Safely quit model.
        """

        # When `showlog` is set, dump the PETSc Log summary (KSP/SNES/TS solver
        # timings, flop counts, MPI reductions) so solver-level performance is
        # captured alongside the wall-clock phase profile. Best-effort: a
        # missing viewer/log API must not break a clean shutdown.
        if getattr(self, "showlog", False) and getattr(self, "log", None) is not None:
            try:
                import petsc4py

                viewer = petsc4py.PETSc.Viewer().createASCII(
                    os.path.join(self.outputDir, "petsc_log.txt"),
                    comm=petsc4py.PETSc.COMM_WORLD,
                )
                petsc4py.PETSc.Log.view(viewer)
                viewer.destroy()
            except Exception:
                pass

        _UnstMesh.destroy_DMPlex(self)

        return
