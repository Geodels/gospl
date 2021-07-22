#!/bin/bash
module load gcc/4.9.3  python/3.6.5 petsc-gcc-mpich/3.11.1  vtk/8.2.0-qt5.8.0

export F90=gfortran
export F77FLAGS=-fPIC
export FCFLAGS=-fPIC

pip3 install vtk --user

python3 setup.py install --user
