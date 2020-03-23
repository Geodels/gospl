import io
import numpy
import setuptools
from os import path

try:
    from numpy.distutils.fcompiler import FCompiler

    def runtime_library_dir_option(self, dir):
        return self.c_compiler.runtime_library_dir_option(dir)

    FCompiler.runtime_library_dir_option = runtime_library_dir_option
except Exception:
    pass

this_directory = path.abspath(path.dirname(__file__))
with io.open(path.join(this_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()


def configuration(parent_package="", top_path=None):
    INCLUDE_DIRS = []
    LIBRARY_DIRS = []
    LIBRARIES = []

    # PETSc
    import os

    PETSC_DIR = os.environ["PETSC_DIR"]
    PETSC_ARCH = os.environ.get("PETSC_ARCH", "")
    from os.path import join, isdir

    if PETSC_ARCH and isdir(join(PETSC_DIR, PETSC_ARCH)):
        INCLUDE_DIRS += [
            join(PETSC_DIR, PETSC_ARCH, "include"),
            join(PETSC_DIR, "include"),
        ]
        LIBRARY_DIRS += [join(PETSC_DIR, PETSC_ARCH, "lib")]
    else:
        if PETSC_ARCH:
            pass
        INCLUDE_DIRS += [join(PETSC_DIR, "include")]
        LIBRARY_DIRS += [join(PETSC_DIR, "lib")]
    LIBRARIES += ["petsc"]

    import os
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
        extra_f90_compile_args=["-fPIC", "-O3"],
        # extra_f90_compile_args = ['-fPIC', '-O0', '-g', '-fbacktrace','-fcheck=all'],
        extra_link_args=["-shared"],
        runtime_library_dirs=LIBRARY_DIRS,
    )

    return config


if __name__ == "__main__":

    from numpy.distutils.core import setup

    setup(
        name="gospl",
        author="Tristan Salles  ",
        author_email="tristan.salles@sydney.edu.au",
        url="https://github.com/Geodels/gospl",
        version="0.1",
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
        ],
        python_requires=">=3.3",
        classifiers=[
            "Programming Language :: Python :: 3.5",
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Python :: 3.7",
            "Programming Language :: Python :: 3.8",
        ],
    )
