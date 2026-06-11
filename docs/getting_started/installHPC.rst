.. _installHPC:


=========================
Installation on HPC
=========================

The following guidelines describe how goSPL can be deployed on High Performance
Computers. Two approaches are covered:

1. **Native virtual-environment install** — load system MPI/PETSc modules and
   install goSPL Python dependencies into a per-user ``venv``.
2. **Singularity/Apptainer container** — a portable container built from the
   goSPL HPC Docker image, run in hybrid-MPI mode so the container uses the
   cluster's native high-performance fabric at runtime.

The examples below target two Australian HPC systems:

- **NCI Gadi** (*"to search for"* in the language of the Ngunnawal people), an
  `Intel-MPI / PBS <https://nci.org.au>`__ CPU-based research supercomputer.
- **Pawsey Setonix**, a Cray/Slurm AMD EPYC system.

The commands will need to be adapted for other HPC systems.


.. _installHPC.native:

Native virtual-environment install (Gadi)
------------------------------------------

.. code-block:: bash

    # Load system modules
    module load intel-mpi/2021.13.1 hdf5/1.12.1p petsc/3.21.3 \
                netcdf/4.9.2p intel-mkl/2023.2.0 python3/3.11.7

    # Expose the system mpi4py linked against Intel MPI
    export PYTHONPATH=$PETSC_BASE/lib/mpich/Intel:$PYTHONPATH

    # Create a virtual environment (system-site-packages inherits numpy, scipy etc.)
    python3 -m venv --system-site-packages ~/envi_gospl
    source ~/envi_gospl/bin/activate

    # Install goSPL Python dependencies
    python3 -m pip cache purge
    python3 -m pip install mpi4py        # 4.0+
    python3 -m pip install netCDF4       # 1.7+
    python3 -m pip install h5py          # 3.11+
    python3 -m pip install vtk           # 9.3+
    python3 -m pip install xarray
    python3 -m pip install ruamel.yaml
    python3 -m pip install numpy-indexed
    python3 -m pip install pyshtools
    python3 -m pip install gflex         # 1.2+
    python3 -m pip install meson_python  # 0.16+

    # Install goSPL from the goSPL repository (cd into your clone)
    python3 -m pip install --no-deps --no-build-isolation .

.. note::

    ``--no-deps`` is correct here because the virtual environment already
    supplies all runtime dependencies via the system modules and the ``pip
    install`` steps above. ``--no-build-isolation`` ensures the build uses the
    MPI and PETSc libraries from the loaded modules.


Testing the installation (Gadi)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Download the archived HPC test example::

    wget https://raw.githubusercontent.com/Geodels/goSPL-examples/main/hpc-setup/gospl_NCI.tar

Set up a PBS job script (replace ``PRJ`` with your NCI project code and
``input-file-name.yml`` with your goSPL input file):

.. code-block:: bash

    #!/bin/bash
    #PBS -l ncpus=48
    #PBS -l mem=60GB
    #PBS -l jobfs=10GB
    #PBS -q normal
    #PBS -P PRJ
    #PBS -l walltime=00:10:00
    #PBS -l storage=gdata/PRJ+scratch/PRJ
    #PBS -l wd

    module load intel-mpi/2021.13.1 hdf5/1.12.1p petsc/3.21.3 \
                netcdf/4.9.2p intel-mkl/2023.2.0 python3/3.11.7
    export PYTHONPATH=$PETSC_BASE/lib/mpich/Intel:$PYTHONPATH

    source ~/envi_gospl/bin/activate

    mpirun -np $PBS_NCPUS python3 runModel.py -i input-file-name.yml

    deactivate

Submit with ``qsub <script>.pbs``.

.. note::

    For an interactive session on Gadi::

        qsub -I -P PRJ -q express -l ncpus=4,mem=16GB,walltime=1:00:00,storage=scratch/PRJ


.. _installHPC.container:

Singularity / Apptainer container
-----------------------------------

Containers provide a self-contained, reproducible goSPL environment that can be
moved between systems without rebuilding. The HPC container uses **hybrid MPI**:
MPICH is compiled into the container as a placeholder, and at runtime the cluster's
native high-performance MPI (Cray on Setonix, Intel MPI on Gadi) is bind-mounted
over it by the ``singularity/...-mpi`` module. This gives full fabric performance
(Slingshot on Setonix, InfiniBand on Gadi) with a single portable image.

.. important::

    **You cannot build Singularity images on Gadi or Setonix** — building
    requires root. Build the image on your local Linux machine, then upload the
    ``.sif`` file to the cluster.

.. important::

    **Ubuntu 24.04 base is mandatory for Setonix.** After Setonix's CPE 25.03
    upgrade, containers based on older Ubuntu versions no longer run correctly.
    The ``geodels/gospl-hpc`` Docker image uses ``ubuntu:24.04``.


Building the container (local Linux machine)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The goSPL HPC container source lives in ``docker/`` in the goSPL repository.
A helper script handles the full build and conversion to ``.sif``:

.. code-block:: bash

    git clone https://github.com/Geodels/gospl.git
    cd gospl/docker/

    # Build from the latest release tag and convert to .sif
    ./build.sh v2026.06.11

    # Build + push to Docker Hub + convert (needed for singularity pull on cluster)
    ./build.sh v2026.06.11 --push

The script runs ``docker build --platform linux/amd64``, smoke-tests the image,
and converts it to ``.sif`` via ``apptainer build`` (or ``singularity build``).
It requires Docker and either ``apptainer`` or ``singularity`` on the build machine.

After the build completes, upload the ``.sif`` to the cluster::

    scp gospl-hpc-v2026.06.11.sif <user>@gadi.nci.org.au:/scratch/<PRJ>/<user>/containers/
    scp gospl-hpc-v2026.06.11.sif <user>@setonix.pawsey.org.au:/scratch/<PRJ>/<user>/containers/

Alternatively, since the image is published on Docker Hub
(`geodels/gospl-hpc <https://hub.docker.com/r/geodels/gospl-hpc>`_), pull it
directly on the cluster — this is a conversion, not a root build, so it works on
a login node::

    # On Setonix
    module load singularity/4.1.0-mpi
    singularity pull gospl-hpc-v2026.06.11.sif docker://geodels/gospl-hpc:v2026.06.11


Running on Gadi with Singularity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Example PBS job script:

.. code-block:: bash

    #!/bin/bash
    #PBS -N gospl_run
    #PBS -P PRJ
    #PBS -q normal
    #PBS -l ncpus=48
    #PBS -l mem=192GB
    #PBS -l walltime=10:00:00
    #PBS -l storage=scratch/PRJ+gdata/PRJ
    #PBS -l wd
    #PBS -j oe

    CONTAINER=/scratch/${PBS_PROJECT}/${USER}/containers/gospl-hpc-v2026.06.11.sif
    INPUT_YML=/scratch/${PBS_PROJECT}/${USER}/runs/my_gospl_input.yml

    module purge
    module load singularity
    module load intel-mpi/2021.13.1

    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    export MKL_NUM_THREADS=1

    mpirun -np ${PBS_NCPUS} \
        singularity exec \
            --bind /scratch/${PBS_PROJECT},/g/data/${PBS_PROJECT} \
            "${CONTAINER}" \
            python3 -c "from gospl.model import Model; m = Model('${INPUT_YML}'); m.runProcesses(); m.destroy()"


Running on Setonix with Singularity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Setonix uses Slurm and the ``singularity/4.1.0-mpi`` module, which automatically
bind-mounts Cray MPI paths (``SINGULARITY_BINDPATH``) so no manual ``--bind``
for MPI libraries is needed.

Example Slurm job script:

.. code-block:: bash

    #!/bin/bash
    #SBATCH --job-name=gospl_run
    #SBATCH --account=PRJ
    #SBATCH --partition=work
    #SBATCH --nodes=2
    #SBATCH --ntasks-per-node=128
    #SBATCH --cpus-per-task=1
    #SBATCH --mem=220G
    #SBATCH --time=10:00:00
    #SBATCH --output=gospl_%j.log

    CONTAINER=/scratch/${PAWSEY_PROJECT}/${USER}/containers/gospl-hpc-v2026.06.11.sif
    INPUT_YML=/scratch/${PAWSEY_PROJECT}/${USER}/runs/my_gospl_input.yml

    module purge
    module load singularity/4.1.0-mpi   # sets SINGULARITY_BINDPATH for Cray MPI

    export OMP_NUM_THREADS=1
    export OPENBLAS_NUM_THREADS=1
    export MKL_NUM_THREADS=1

    srun --export=ALL \
        singularity exec \
            --bind /scratch/${PAWSEY_PROJECT} \
            "${CONTAINER}" \
            python3 -c "from gospl.model import Model; m = Model('${INPUT_YML}'); m.runProcesses(); m.destroy()"

.. note::

    ``srun`` is preferred over ``mpirun`` on Setonix because it integrates with
    Cray PMI for process management and network fabric setup.

.. warning::

    **Known Setonix issue: parallel HDF5 I/O inside Singularity.** Pawsey has a
    known issue with MPI-parallel collective HDF5 writes inside Singularity
    containers (under investigation). If you encounter HDF5 errors on Setonix,
    check `Pawsey Known Issues <https://support.pawsey.org.au>`__ for the latest
    status and available workarounds. Do **not** switch to the ``nompi`` build
    of ``h5py`` — goSPL requires MPI-linked HDF5 for collective writes.

    Do **not** kill the job with ``SIGKILL`` before goSPL finishes — always allow
    ``model.destroy()`` to run so that PETSc finalises cleanly.
