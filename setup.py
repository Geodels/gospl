from setuptools import find_packages
from numpy.distutils.core import setup, Extension

try:
    from distutils.command import bdist_conda
except ImportError:
    pass


import os
import io
import subprocess
import platform

# in development set version to none and ...
PYPI_VERSION = "0.1.14"


def git_version():
    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ["SYSTEMROOT", "PATH"]:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env["LANGUAGE"] = "C"
        env["LANG"] = "C"
        env["LC_ALL"] = "C"
        out = subprocess.Popen(cmd, stdout=subprocess.PIPE, env=env).communicate()[0]
        return out

    try:
        out = _minimal_ext_cmd(["git", "rev-parse", "--short", "HEAD"])
        GIT_REVISION = out.strip().decode("ascii")
    except OSError:
        GIT_REVISION = "Unknown"

    return GIT_REVISION


if PYPI_VERSION is None:
    PYPI_VERSION = git_version()


ext = Extension(
    name="gospl._fortran", sources=["fortran/functions.pyf", "fortran/functions.f90"]
)


this_directory = os.path.abspath(os.path.dirname(__file__))
with io.open(os.path.join(this_directory, "README.md"), encoding="utf-8") as f:
    long_description = f.read()


if __name__ == "__main__":
    setup(
        name="gospl",
        author="Tristan Salles",
        author_email="tristan.salles@sydney.edu.au",
        url="https://github.com/Geodels/gospl",
        version=PYPI_VERSION,
        description="A Python interface to perform Global Landscape Evolution Model",
        long_description=long_description,
        long_description_content_type="text/markdown",
        ext_modules=[ext],
        packages=["gospl", "gospl.tools", "gospl.flow", "gospl.mesher", "gospl.sed"],
        install_requires=[
            "pytest",
            "Cython==0.29.21",
            "numpy==1.19.2",
            "scipy==1.5.2",
            "h5py==2.10.0",
            "pandas==1.1.2",
            "ruamel.yaml==0.16.12",
            "meshplex==0.13.3",
            "mpi4py>=3.0.3",
            "petsc4py>=3.13.0",
            "pre-commit==2.7.1",
            "fastfunc==0.2.3",
            "vtk>=9.0.0",
            "numpy-indexed==0.3.5",
            "meshio==4.2.1",
            "scikit-fuzzy>=0.4.2",
        ],
        classifiers=[
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Python :: 3.7",
            "Programming Language :: Python :: 3.8",
        ],
    )
