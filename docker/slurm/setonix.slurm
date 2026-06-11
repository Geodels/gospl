#!/bin/bash --login
###############################################################################
# Pawsey Setonix (Slurm) job script for goSPL inside a Singularity container.
# Hybrid MPI: Setonix's Cray MPI (MPICH ABI) is bind-mounted over the
# container's MPICH. Ubuntu 24.04 base is mandatory (CPE 25.03+).
# See AGENTS.md > "Docker / HPC Container".
#
# Edit --account, the resource requests, and CONTAINER/INPUT before submitting:
#   sbatch docker/slurm/setonix.slurm
###############################################################################
#SBATCH --account=<project>
#SBATCH --partition=work
#SBATCH --nodes=1
#SBATCH --ntasks=128
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --job-name=gospl
#SBATCH --exclusive

set -euo pipefail

# The -mpi flavour sets SINGULARITY_BINDPATH + SINGULARITYENV_LD_LIBRARY_PATH
# automatically — no manual --bind for MPI paths is needed (AGENTS.md).
module load singularity/4.1.0-mpi

# One thread per MPI rank — AGENTS.md container invariant 2.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# Path to the .sif (scp'd from your local build) and the goSPL input YAML.
CONTAINER="${CONTAINER:-${MYSCRATCH}/gospl-hpc-v2026.06.11.sif}"
INPUT="${INPUT:-input.yml}"

# srun (not mpirun) is the Pawsey-recommended launcher with the -mpi module.
# model.destroy() is required for clean PETSc finalization (AGENTS.md
# container invariant 4) — replace the -c driver with your own script if it
# already calls runProcesses() + destroy().
srun --nodes="${SLURM_NNODES}" --ntasks="${SLURM_NTASKS}" \
     --cpus-per-task="${SLURM_CPUS_PER_TASK}" \
    singularity exec "${CONTAINER}" \
        python -c "from gospl.model import Model; m = Model('${INPUT}'); m.runProcesses(); m.destroy()"

# NOTE: Pawsey has an active known issue with MPI-parallel HDF5 collective I/O
# inside Singularity on Setonix. If goSPL fails during HDF5 output, check the
# Pawsey Known Issues page; do NOT switch to a nompi h5py (AGENTS.md).
