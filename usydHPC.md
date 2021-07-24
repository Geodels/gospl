# goSPL on USYD HPC

## Get to Artemis

When not on Sydney campus use the VPN

```bash
term@local$ ssh usydname@hpc.sydney.edu.au
```

## Load the required module

```bash
usyd@hpc$ module load gcc/4.9.3 python/3.6.5 petsc-gcc-mpich/3.11.1
```

## Cloning goSPL

### Easy installation

Create a repository named `codes` in your home folder then

```bash
usyd@hpc$ mkdir ~/codes
usyd@hpc$ cd codes
usyd@hpc$ git clone https://github.com/Geodels/gospl.git
usyd@hpc$ cd gospl
usyd@hpc$ ./install_HPC.sh
```

If it fails you might have to change the permissions of the install_HPC.sh script to execute it:

```bash
usyd@hpc$ chmod 777 install_HPC.sh
usyd@hpc$ ./install_HPC.sh
```

If it fails again then follow the steps below:

### Harder way

In Artemis, the `petsc` and `petsc4py` (version 3.11.0) have been installed with `Python3.6.5` and we will need to install `goSPL` with this version. You can check that you have the right `Python` version:

```bash
usyd@hpc$ which python
/usr/local/python/3.6.5/bin/python
```

Then you might want to install manually all the dependencies and see the ones that fails. The libraries versions that work are listed below:

| Package | version nb |  |  |  |
|:-------:|:----------:|:-:| |  |
|  numpy  |  1.19.5    | | h5py |  2.9.0 |
|  numpy-indexed  |  0.3.5 | | pandas | 1.0.3 |
|  scipy  |  1.5.4    | | ruamel.yaml | 0.16.5 |
|  mpi4py  |  3.0.0    | | meshio | 3.3.1 |
|  petsc4py  |  3.11.0    | | meshplex | 0.11.6 |
|  scikit-fuzzy  |  0.4.2    |  | vtk | 9.0.1 |


To install them successively run:

```bash
pip install Cython==0.29.22 --user
pip install numpy==1.19.5 --user
pip install meshplex==0.11.6 --user
pip install meshio==1.0.3 --user
pip install vtk==9.0.1 --user
pip install h5py==2.9.0 --user
pip install numpy-indexed==0.3.5 --user
pip install pandas==1.0.3 --user
pip install scipy==1.5.4 --user
pip install ruamel.yaml==0.16.5 --user
pip install mpi4py==3.0.0 --user
pip install petsc4py==3.11.0 --user
pip install scikit-fuzzy==0.4.2 --user
```

You might want to report any issue to goSPL developper.

If it works then you can run the local installation

```bash
usyd@hpc$ export F90=gfortran
usyd@hpc$ export F77FLAGS=-fPIC
usyd@hpc$ export FCFLAGS=-fPIC
usyd@hpc$ python setup install --user
```

## Running simulations

Example of simple `PBS` script:

```
#!/bin/bash

# Project
#PBS -P BGH

# 16 CPUs
#PBS -l select=2:ncpus=8:mpiprocs=8:mem=60GB

#PBS -l walltime=01:00:00
#PBS -M xxxx.xxxx@sydney.edu.au
#PBS -m abe
#PBS -q alloc-dm

# Set up environment
module load gcc/4.9.3 python/3.6.5 petsc-gcc-mpich/3.11.1

cd $PBS_O_WORKDIR
cd NorthAmerica

# Launching the job!
mpirun -np 16 python runModel.py -i gospl-20Ma-dt20k.yml
```

then use the `qsub` command to launch and `qstat` to check your place in the queue.
