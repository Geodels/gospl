import petsc4py,sys
petsc4py.init(sys.argv)
from time import clock
from mpi4py import MPI

from .mesher import UnstMesh as _UnstMesh
from .tools import ReadYaml as _ReadYaml
from .tools import WriteMesh as _WriteMesh
from .flow import SPMesh as _SPMesh

from petsc4py import PETSc as _PETSc

MPIrank = MPI.COMM_WORLD.Get_rank()

# refinement = 8 dist = 32606. npoints   664578
# refinement = 9 dist = 15794. npoints  2639874
# refinement = 10 dist = 7625. npoints 10522626


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
            if self.showlog:
                self.log = _PETSc.Log()
                self.log.begin()

            self.modelRunTime = clock()
            self.verbose = True #verbose

            # Read input dataset
            _ReadYaml.__init__(self, filename)

            # Initialise output mesh
            _WriteMesh.__init__(self)

            # Define unstructured mesh
            _UnstMesh.__init__(self)

            # Get external forces
            _UnstMesh.applyForces(self)

            # Surface processes initialisation
            _SPMesh.__init__(self,*args, **kwargs)

            if MPIrank == 0:
                print('--- Initialisation Phase (%0.02f seconds)'% (clock() - self.modelRunTime))

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

            while(self.tNow<=self.tEnd):
                tstep = clock()

                # Check paleomap forcing
                _UnstMesh.updatePaleomap(self)

                # Compute Flow Accumulation
                _SPMesh.FlowAccumulation(self)

                # Output time step for first step
                if self.tNow == self.tStart:
                    _WriteMesh.outputMesh(self)
                    self.saveTime += self.tout

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

                # Output time step
                if self.tNow >= self.saveTime:
                    _WriteMesh.outputMesh(self)
                    self.saveTime += self.tout

                # Update Tectonic, Sea-level & Climatic conditions
                if self.tNow < self.tEnd:
                    _UnstMesh.applyForces(self)

                # Advance time
                self.tNow += self.dt

                if MPIrank == 0:
                    print('--- Computational Step (%0.02f seconds)'% (clock() - tstep))

            return

        def destroy(self):
            """
            Destroy PETSc DMPlex objects and associated Petsc local/global Vectors and Matrices.
            Safely quit eSCAPE model.
            """

            _UnstMesh.destroy_DMPlex(self)

            if self.showlog:
                self.log.view()

            if MPIrank == 0:
                print('\n+++\n+++ Total run time (%0.02f seconds)\n+++'% (clock() - self.modelRunTime))

            return

    return LandscapeEvolutionModelClass(filename, *args, **kwargs)
