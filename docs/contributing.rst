.. _contributing:

=========================
Contributing to goSPL
=========================

Contributions of any kind to goSPL are more than welcome. That does not mean new code only, but also improvements of documentation and user guide, additional tests (ideally filling the gaps in existing suite) or bug report or idea what could be added or done better.

All contributions should go through our GitHub repository. Bug reports, ideas or even questions should be raised by opening an issue on the GitHub tracker. Suggestions for changes in code or documentation should be submitted as a pull request. However, if you are not sure what to do, feel free to open an issue.

All discussion will then take place on GitHub to keep the development of
goSPL transparent.

If you decide to contribute to the codebase, ensure that you are using an
up-to-date `master` branch. The latest development version will always be there, including the documentation (powered by `sphinx`_).


Seven Steps for Contributing
----------------------------

There are seven basic steps to contributing to goSPL:

1. Fork the goSPL git repository
2. Create a development environment
3. Install goSPL dependencies
4. Make a development build of goSPL
5. Update the documentation
6. Format code
7. Submit a Pull Request

Each of the steps is detailed below.

1. Fork the goSPL git repository
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Git can be complicated for new users, but you no longer need to use command line to work with ``git``. If you are not familiar with ``git``, we recommend using tools on GitHub.org, GitHub Desktop or tools with included git like Atom. However, if you want to use command line, you can fork goSPL repository using following::

    git clone git@github.com:your-user-name/gospl.git gospl-yourname
    cd gospl-yourname
    git remote add upstream git://github.com/Geodels/gospl

This creates the directory gospl-yourname and connects your repository to
the upstream (main project) gospl repository.

Then simply create a new branch of master branch.


2. Create a development environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A development environment is a virtual space where you can keep an independent installation of goSPL. This makes it easy to keep both a stable version of Python in one place you use for work, and a development version (which you may break while playing with the code) in another.

.. note:: 
    
    For Windows user, this is recommended to first install `WSL <https://learn.microsoft.com/en-us/windows/wsl/install>`_ as some of the required libraries are not available as windows packages on Anaconda. Then you should be able to proceed as detailed below.


An easy way to create a goSPL development environment is as follows:

- Install either `Anaconda <http://docs.continuum.io/anaconda/>`_ or
  `miniconda <http://conda.pydata.org/miniconda.html>`_
- Make sure that you have cloned the repository
- ``cd`` to the goSPL source directory

Tell conda to create a new environment, named ``gospl_dev``, or any other name you would like for this environment, by running::

      conda create -n gospl_dev python=3.11

This will create the new environment, and not touch any of your existing environments, nor any existing Python installation.

To work in this environment, macOS and Linux users should use::

      conda activate gospl_dev

You will then see a confirmation message to indicate you are in the new development environment.

To view your environments::

      conda info -e

To return to you home root environment::

      conda deactivate

See the full conda docs `here <http://conda.pydata.org/docs>`__.

At this point you can easily do a *development* install, as detailed in the next sections.

3. Installing Dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^

To run goSPL in an development environment, you must first install
goSPL's dependencies. We suggest doing so using the following commands
(executed after your development environment has been activated)
to ensure compatibility of all dependencies::

    conda config --env --add channels conda-forge
    conda config --env --add channels defaults
    conda config --env --set channel_priority strict
    conda install numpy pip scipy numpy-indexed cython compilers
    conda install pandas h5py meshio ruamel.yaml
    conda install vtk mpi4py petsc4py gflex

.. note::

    This should install all necessary dependencies. Alternatively you could directly use the `environment.yml <https://raw.githubusercontent.com/Geodels/gospl/master/environment.yml>`_ file provided in the goSPL GitHub repository::

        conda env create -f environment.yml


4. Making a development build
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once dependencies are in place, make an in-place build by navigating to the git clone of the goSPL repository and running::

    pip install --no-build-isolation -e .

This will install goSPL into your environment but allows any further changes without the need of reinstalling new version.

5. Updating the Documentation and User Guide
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

goSPL documentation resides in the `docs` folder. Changes to the docs are
make by modifying the appropriate file within `docs`.

goSPL docs us reStructuredText syntax, `which is explained here <http://www.sphinx-doc.org/en/stable/rest.html#rst-primer>`_ and the docstrings follow the `Numpy Docstring standard <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_.

Once you have made your changes, you may try if they render correctly by building the docs using sphinx.

To do so, you can navigate to the doc folder and type::

    READTHEDOCS=True make clean html

The resulting html pages will be located in ``docs/build/html``. In case of any errors,
you can try to use make html within a new environment based on requirements.txt specification in the doc folder.

For minor updates, you can skip whole make html part as reStructuredText syntax is
usually quite straightforward.


6. Formatting the code
^^^^^^^^^^^^^^^^^^^^^^

Python (PEP8 / black)
~~~~~~~~~~~~~~~~~~~~~

goSPL follows the `PEP8 <http://www.python.org/dev/peps/pep-0008/>`_ standard
and uses `Black`_ to ensure a consistent code format throughout the project.

Travis CI will run ``black --check`` and fails if there are files which would be
auto-formatted by ``black``. Therefore, it is helpful before submitting code to
auto-format your code::

    black gospl

Additionally, many editors have plugins that will apply ``black`` as you edit files.
If you don't have black, you can install it using pip::

    pip install black

7. Submitting a Pull Request
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once you've made changes and pushed them to your forked repository, you then
submit a pull request to have them integrated into the goSPL code base.

You can find a pull request (or PR) tutorial in the `GitHub's Help Docs <https://help.github.com/articles/using-pull-requests/>`_.

.. _sphinx: https://www.sphinx-doc.org/

.. _Black: https://black.readthedocs.io/en/stable/
