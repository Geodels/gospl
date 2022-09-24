#!/bin/bash
module load gcc/4.9.3  python/3.6.5 petsc-gcc-mpich/3.11.1

export F90=gfortran
export F77FLAGS=-fPIC
export FCFLAGS=-fPIC

pip install Cython==0.29.22 --user
pip install numpy==1.19.5 --user
pip install meshplex==0.11.6 --user
pip install meshio==3.3.1 --user
pip install vtk==9.0.1 --user
pip install h5py==2.9.0 --user
pip install numpy-indexed==0.3.5 --user
pip install pandas==1.0.3 --user
pip install scipy==1.5.4 --user
pip install ruamel.yaml==0.16.5 --user
pip install mpi4py==3.0.0 --user
pip install petsc4py==3.11.0 --user
pip install scikit-fuzzy==0.4.2 --user
pip install xarray==0.16.2 --user

mv hpcrequirements.txt requirements.txt

python3 setup.py install --user
