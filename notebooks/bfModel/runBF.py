from mpi4py import MPI
from gospl.model import Model as sim

from scripts import mergeBack as merger

MPIrank = MPI.COMM_WORLD.Get_rank()
MPIcomm = MPI.COMM_WORLD

forin = "forward.yml"
backin = ["backward15Ma.yml", "backward10Ma.yml"]
backout = "output-backward"

# Running the backwards models by periods
for k in range(len(backin)):
    mod = sim(backin[k], False, False)
    mod.runProcesses()
    mod.destroy()
    if MPIrank == 0:
        print("", flush=True)

# Merging all backward models into a single outputs
if MPIrank == 0:
    merger.mergeBackModels(backin, backout)
    print("", flush=True)
MPIcomm.Barrier()

# Running the forward model forced with backward simulations
mod = sim(forin, False, False)
mod.runProcesses()
mod.destroy()
