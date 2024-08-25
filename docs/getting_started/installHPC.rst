.. _installHPC:


=========================
Installation on HPC
=========================

The following guidelines show how goSPL can be installed on **High Performance Computers**. The approach presented here might need to be adapted depending on the specs of the HPC system you are able to access.

Here, we use **Gadi** (*meaning ‘to search for’ in the language of the Ngunnawal people*), one of `NCI <https://nci.org.au/news-events/news/australias-gadi-a-recognised-powerhouse-global-supercomputing-ranking>`__ CPU-based research supercomputer. 


Installing dependencies
-----------------------------

.. code-block:: bash

    # Install goSPL environment libraries (doing this in my home directory)
    module load intel-mpi/2021.13.0 hdf5/1.12.1p petsc/3.21.3 netcdf/4.9.2p intel-mkl/2023.2.0 python3/3.11.7
    export PYTHONPATH=$PETSC_BASE/lib/mpich/Intel:$PYTHONPATH # access to mpi4py lib

    # Virtual environment creation
    python3 -m venv --system-site-packages /home/XXX/USERXXX/envi_gospl
    source /home/XXX/USERXXX/envi_gospl/bin/activate

    # Required packages
    python3 -m pip cache purge
    python3 -m pip install pyarrow        # (17.0.0)
    python3 -m pip install  mpi4py        # (4.0.0)
    python3 -m pip install netCDF4        # (1.7.1)
    python3 -m pip install  h5py          # (3.11.0)
    python3 -m pip install vtk            # (9.3.1)
    python3 -m pip install xarray         # (2024.7.0)
    python3 -m pip install meshio         # (5.3.5)
    python3 -m pip install ruamel.yaml    # (0.18.6)
    python3 -m pip install numpy-indexed  # (0.3.7)
    python3 -m pip install rioxarray      # (0.17.0)
    python3 -m pip install gflex          # (1.2.0)
    python3 -m pip install meson_python   # (0.16.0)

    # Installing isoFlex (cd to isoFlex repo)
    python3 -m pip install --no-deps .

    # Installing goSPL (cd to goSPL repo)
    python3 -m pip install --no-deps .


Testing the installation
-----------------------------

Get the archived codes and testing examples for HPC:

.. code-block:: bash    

    wget https://raw.githubusercontent.com/Geodels/goSPL-examples/main/hpc-setup/gospl_NCI.tar


Set-up your PBS script (where `PRJ` is your project name and `input-file-name.yml` is your goSPL input file name):

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

    # Load required modules
    module load intel-mpi/2021.13.0 hdf5/1.12.1p petsc/3.21.3 netcdf/4.9.2p intel-mkl/2023.2.0 python3/3.11.7
    export PYTHONPATH=$PETSC_BASE/lib/mpich/Intel:$PYTHONPATH 

    # Activate the virtual environment
    source /home/XXX/USERXXX/envi_gospl/bin/activate

    # Simulation run
    mpirun -np $PBS_NCPUS python3 runModel.py -i input-file-name.yml

    deactivate
