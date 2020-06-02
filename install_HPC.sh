#!/bin/bash
module load gcc/4.9.3  python/3.6.5 petsc-gcc-mpich/3.11.1

export F90=gfortran
export F77FLAGS=-fPIC
export FCFLAGS=-fPIC

pip3 install -U fastfunc
pip3 install -U meshplex

python3 setup.py install --user
