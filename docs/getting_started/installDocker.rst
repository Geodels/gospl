.. _installDocker:


=========================
Installation via Docker
=========================


Getting Docker
--------------

Docker is an open platform for developing, shipping, and running applications.
Install it for your operating system from:

+  `https://docs.docker.com/get-docker/ <https://docs.docker.com/get-docker/>`_

Once Docker Desktop is installed, the Docker CLI is available in your terminal.


goSPL on Docker Hub
-------------------

The `geodels/gospl-examples <https://hub.docker.com/r/geodels/gospl-examples>`_
image packages the full ``gospl-smoke`` conda environment — goSPL, PETSc, MPI,
GMT, VTK and JupyterLab — and is the recommended way to run the examples without
a local conda install.

.. important::

  The ``geodels/gospl-examples`` image ships the runtime environment only
  (no notebooks). You run it against a local clone of the
  `goSPL-examples <https://github.com/Geodels/goSPL-examples>`_ repository
  mounted at ``/work`` inside the container.

The image is built automatically from ``environment.yml`` (via ``mamba``) by the
``Build and push Docker image`` GitHub Actions workflow on every release tag. It
targets goSPL ``v2026.06.11`` and is published as a **multi-arch manifest**
(``linux/amd64`` + ``linux/arm64``).

.. note::

  **Apple Silicon (M1/M2/M3/M4) users.** A plain ``docker pull`` automatically
  selects the native ``linux/arm64`` image on Apple Silicon — no ``--platform``
  flag and no Rosetta emulation needed. Running the native arm64 image is
  **significantly faster** than forcing amd64 emulation. Confirm you got the
  native build with::

    docker run --rm geodels/gospl-examples:latest \
      python -c "import platform; print(platform.machine())"
    # expect: aarch64 on Apple Silicon, x86_64 on Intel


Included packages
^^^^^^^^^^^^^^^^^

On top of the goSPL runtime dependencies, the image includes:

- **MPI / HPC**: ``mpi4py``, ``petsc4py``, ``h5py`` (MPI-linked), ``netCDF4``
- **Mesh / geometry**: ``stripy``, ``meshplex``, ``jigsawpy``, ``triangle``, ``uxarray``
- **Geoscience**: ``pygmt``, ``cartopy``, ``pyproj``, ``xesmf``, ``rasterio``
- **Visualisation**: ``vtk``, ``pyvista``, ``pyevtk``, ``xarray``, ``seaborn``
- **Notebooks**: JupyterLab

.. note::

  Depending on your operating system, configure Docker Desktop's resource
  allocation (CPUs, memory, swap, disk image size) to match the size of the
  simulations you want to run.


Main command lines
------------------

Pulling the image
^^^^^^^^^^^^^^^^^

Clone the examples repository, then pull the image::

  git clone https://github.com/Geodels/goSPL-examples.git
  cd goSPL-examples
  docker pull geodels/gospl-examples:latest

You can list all images on your system::

  docker images

An image can be removed::

  docker rmi geodels/gospl-examples:latest


Starting the container from a terminal
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Mount the cloned repository at ``/work`` and start JupyterLab::

  docker run -it --rm -p 8888:8888 -v "$PWD":/work geodels/gospl-examples:latest

JupyterLab will start and print a ``http://127.0.0.1:8888/lab?token=…`` URL —
open it in your browser. The full repository is visible under ``/work`` and the
``gospl-smoke`` environment is already active (no ``conda activate`` needed).

.. note::

  **Apple Silicon users**: no ``--platform`` flag is needed — Docker selects the
  native ``linux/arm64`` image automatically.

.. important::

  In the Docker Desktop *Settings → Resources* panel, increase the CPU count,
  memory and swap limits to ensure simulations do not run out of resources.

You can list all containers on your machine::

  docker ps -a

Stop a running container (releases CPU/RAM)::

  docker stop <container-id>

Remove a stopped container::

  docker rm <container-id>

.. important::

  Keep track of how many containers have been created — they accumulate and can
  consume significant disk space. Use ``docker ps -a`` and ``docker rm`` to clean
  up containers you no longer need.


Running the examples
--------------------

Each example in the repository follows a three-stage workflow:

1. **Build inputs** — run all cells of ``build_inputs.ipynb`` /
   ``model_setup.ipynb`` to generate the mesh and forcing files.
2. **Run the model** — open a terminal in JupyterLab (*File ▸ New ▸ Terminal*)
   and launch goSPL under MPI.
3. **Analyse outputs** — run all cells of ``sims-analysis.ipynb`` /
   ``extract_strata.ipynb`` to remap and visualise the results.

Example — running the ``stratigraphic_record`` example inside the container::

  cd /work/Local-examples/stratigraphic_record
  mpirun -np 4 python runModel.py -i input-strati.yml

Replace ``-np 4`` with the number of MPI ranks you want to use.

.. note::

  **Performance — threading.** The image sets ``OMP_NUM_THREADS``,
  ``OPENBLAS_NUM_THREADS``, ``MKL_NUM_THREADS`` and ``NUMEXPR_NUM_THREADS``
  to ``1`` by default. goSPL gets its parallelism from MPI ranks (``mpirun
  -np N``), so single-threaded BLAS per rank avoids CPU oversubscription.
  Set ``-np`` to the number of physical cores you want to use. For a
  non-MPI workload you can re-enable threading with
  ``docker run -e OMP_NUM_THREADS=8 …``.


Quick end-to-end test (headless)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``stratigraphic_record`` example ships pre-built inputs so you can verify
the full stack without opening a notebook::

  docker run -it --rm -v "$PWD":/work geodels/gospl-examples:latest \
    bash -lc "cd /work/Local-examples/stratigraphic_record && \
              mpirun -np 4 python runModel.py -i input-strati.yml"

A successful run prints goSPL's per-step progress and writes HDF5 outputs into
the example folder, ready for post-processing with ``extract_strata.ipynb``.


Using Docker Desktop Dashboard
--------------------------------

As an alternative to the terminal, you can use the Docker Desktop Dashboard to
pull and start the ``geodels/gospl-examples`` image. When starting the container
from the Dashboard:

- **Ports**: map host port ``8888`` to container port ``8888``.
- **Volumes** *Host Path*: set to your local clone of ``goSPL-examples``.
- **Volumes** *Container Path*: set to ``/work``.

Then open ``http://localhost:8888`` in your browser to access JupyterLab.

.. note::

  If you need additional Python packages, install them from the JupyterLab
  terminal — the ``gospl-smoke`` environment is already active, so
  ``pip install <package>`` or ``conda install <package>`` will work directly.
