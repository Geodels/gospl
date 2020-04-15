##############################
# BUILDING PYTHON PACKAGE PiPY
##############################
# python3 -m pip install --user --upgrade setuptools wheel
# python3 setup.py sdist
# python3 -m pip install --user --upgrade twine
#  ~/.local/bin/twine check dist/*
#  ~/.local/bin/twine upload dist/*
##############################

import os
import io
from setuptools import setup

try:
    from numpy.distutils.fcompiler import FCompiler

    def runtime_library_dir_option(self, dir):
        return self.c_compiler.runtime_library_dir_option(dir)

    FCompiler.runtime_library_dir_option = runtime_library_dir_option
except Exception:
    pass

this_directory = os.path.abspath(os.path.dirname(__file__))
with io.open(os.path.join(this_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

# get F90 environment variable
# ----------------------------------------
F90 = os.getenv("F90")


# raise error if F90 not defined
# !! comment out this if statement for manual install !!
# ------------------------------------------------------------
if F90 is None or F90 == "":
    l1 = "gospl requires environment variable F90 to be set. \n "
    l2 = 'Please set to one of {"ifort", "gfortran"}'
    raise RuntimeError(l1 + l2)


# specialize for different compilers
# ------------------------------------------------------------
if F90 == "ifort":
    f90_flags = [
        "-fPIC",
        "-xHost",
        "-O3",
        "-ipo",
        "-funroll-loops",
        "-heap-arrays",
        "-mcmodel=medium",
    ]

elif F90 == "gfortran":
    f90_flags = [
        "-fPIC",
        "-O3",
        "-fbounds-check",
        "-mtune=native",
    ]

elif F90 in ["pgfortran", "pgf90", "pgf95"]:
    f90_flags = ["-mp"]

else:
    l1 = "F90 = " + F90 + ". \n"
    l2 = "Environment variable F90 not recognized.  \n"
    raise RuntimeError(l1 + l2)


def configuration(parent_package="", top_path=None):
    INCLUDE_DIRS = []
    LIBRARY_DIRS = []
    LIBRARIES = []

    # PETSc
    PETSC_DIR = os.environ["PETSC_DIR"]
    PETSC_ARCH = os.environ.get("PETSC_ARCH", "")

    if PETSC_ARCH and os.path.isdir(os.path.join(PETSC_DIR, PETSC_ARCH)):
        INCLUDE_DIRS += [
            os.path.join(PETSC_DIR, PETSC_ARCH, "include"),
            os.path.join(PETSC_DIR, "include"),
        ]
        LIBRARY_DIRS += [os.path.join(PETSC_DIR, PETSC_ARCH, "lib")]
    else:
        if PETSC_ARCH:
            pass
        INCLUDE_DIRS += [os.path.join(PETSC_DIR, "include")]
        LIBRARY_DIRS += [os.path.join(PETSC_DIR, "lib")]
    LIBRARIES += ["petsc"]

    import petsc4py

    INCLUDE_DIRS += [petsc4py.get_include()]

    # Configuration
    from numpy.distutils.misc_util import Configuration

    config = Configuration("", parent_package, top_path)

    config.add_extension(
        "gospl._fortran",
        sources=["fortran/functions.pyf", "fortran/functions.F90"],
        depends=["fortran/functionsmodule.h"],
        define_macros=[],  # [('F2PY_REPORT_ON_ARRAY_COPY',0)],
        include_dirs=INCLUDE_DIRS + [os.curdir],
        libraries=LIBRARIES,
        library_dirs=LIBRARY_DIRS,
        extra_f90_compile_args=f90_flags,
        # extra_f77_compile_args=["-fPIC", "-O3", "-Wunused-variable"],
        # extra_f90_compile_args=[
        #     "-fPIC",
        #     "-O3",
        #     "-Wunused-variable",
        #     "-Wincompatible-pointer-types",
        #     "-Wcpp",
        #     "-Wunused-function",
        # ],
        # extra_f90_compile_args = ['-fPIC', '-O0', '-g', '-fbacktrace','-fcheck=all'],
        extra_link_args=["-shared"],
        runtime_library_dirs=LIBRARY_DIRS,
    )

    return config


if __name__ == "__main__":

    import numpy.distutils.core

    numpy.distutils.core.setup(
        name="gospl",
        author="Tristan Salles  ",
        author_email="tristan.salles@sydney.edu.au",
        url="https://github.com/Geodels/gospl",
        version="0.1.2",
        description="A Python interface to perform Global Landscape Evolution Model",
        long_description=long_description,
        long_description_content_type="text/markdown",
        configuration=configuration,
        packages=["gospl", "gospl.tools", "gospl.mesher", "gospl.flow"],
        install_requires=[
            "pytest",
            "Cython>=0.28.5",
            "numpy>=1.15.1",
            "scipy>=1.1.0",
            "h5py>=2.8.0",
            "pandas>=0.17.1",
            "ruamel.yaml>=0.15.64",
            "meshplex==0.11.6",
            "mpi4py>=3.0.0",
            "petsc4py>=3.8.1",
            "pre-commit",
            "fastfunc==0.2.2",
        ],
        python_requires=">=3.3",
        classifiers=[
            "Programming Language :: Python :: 3.5",
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Python :: 3.7",
            "Programming Language :: Python :: 3.8",
        ],
    )
