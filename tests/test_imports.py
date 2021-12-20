import pytest
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)


def test_scipy_import():
    import scipy

    print("\t\t You have scipy version {}".format(scipy.__version__))


def test_mpi4py_import():
    import mpi4py
    from mpi4py import MPI

    print("\t\t You have mpi4py version {}".format(mpi4py.__version__))


def test_petsc4py_import():
    import sys
    import petsc4py
    from petsc4py import PETSc

    petsc4py.init(sys.argv)
    print("\t\t You have petsc4py version {}".format(petsc4py.__version__))


def test_h5py_import():
    import h5py

    print("\t\t You have h5py version {}".format(h5py.__version__))


def test_pandas_import():
    import pandas

    print("\t\t You have pandas version {}".format(pandas.__version__))


def test_yaml_import():
    import ruamel.yaml

    print("\t\t You have yaml version {}".format(ruamel.yaml.__version__))


def test_meshio_import():
    import meshio

    print("\t\t You have meshio version {}".format(meshio.__version__))


def test_meshplex_import():
    import meshplex

    print("\t\t You have meshio version {}".format(meshplex.__version__))


# def test_skfuzzy_import():
#     import skfuzzy
#
#     print("\t\t You have skfuzzy version {}".format(skfuzzy.__version__))


def test_vtk_import():
    import vtk

    # print("\t\t You have vtk version {}".format(vtk.__version__))


def test_gospl_modules():
    import gospl
    from gospl.model import Model


def test_jupyter_available():
    from subprocess import check_output

    try:
        result = str(check_output(["which", "jupyter"]))[2:-3]

        print("\t\t You have jupyter version {}".format(result))

    except Exception:
        print("Jupyter not installed")
        print("Jupyter is needed to run the example documentation")
