#################
Hands-on examples
#################

To get started on how to use *gospl*, a series of examples and tutorials is provided in the following `Examples`_ folder.

.. warning::
  Be aware that the `Examples`_ folder is quite large as it comes with several global input files...


Launching Docker
------------------------

The tutorials proposed in the `Examples`_ repository are there to help you build the input files (topography, climate and precipitation meshes). The easiest way to get started is by using the docker container (see the installation page to start):

.. code:: bash

  docker pull geodels/gospl:v0.1.3

And to attach a volume to the container.

.. code:: bash

  docker run -d \
     --name my_container \
     -p 8888:8888 \
     -v my_vol:/live/share \
     geodels/gospl

Then you will access the container via your browser at the following address: :code:`http://localhost:8888`.

The landing page is composed of 2 folders: :code:`lib` & :code:`share`. You will be able to access the `Examples`_ from your local folder via the :code:`share` folder.


Available Jupyter Notebooks
---------------------------

In the `Examples`_ folder, you will find 3 notebooks:

a. :code:`A-create-base-mesh.ipynb` that is used to build the initial spherical meshes used in **gospl** based on netCDF files. A refinement level can be set to increase the mesh resolution. For simulation on a local computer one will use a refinement level of *8* to *9*. For higher models a refinement of *10* or *11* is necessary but will need to be run on a super computer.

b. :code:`B-map-paleodata.ipynb` that builds the different compressed Numpy arrays (**.npz**) from paleo-topography and precipitation maps.

c. :code:`C-create-disp.ipynb` that builds the displacements mesh over time using 3D cartesian velocities again as compressed Numpy arrays.

.. important::
  Paleotopography maps can be obtained from `Scotese`_ Project as well as precipitation grids. Another good source is from the `gPlates`_ portal where velocities and dynamic topography maps can also be exported.


Here are some *indicative* wall times for HR models ran for a 10 million years simulation:

.. list-table:: :code:`ref_lvl` set to 10 (min:7.6 km, max:10.3 km, mean: 9.1 km)
   :widths: 40 40 40 40
   :header-rows: 1

   * - pts number
     - CPU number
     - Wall time
     - CPU nodes
   * - 10612062
     - 8
     - 17:08:58
     - 1
   * - 10612062
     - 16
     - 08:58:32
     - 2
   * - 10612062
     - 32
     - 04:57:58
     - 4
   * - 10612062
     - 64
     - 03:32:15
     - 6
   * - 10612062
     - 128
     - 02:40:53
     - 7

.. list-table:: :code:`ref_lvl` set to 11 (min:4.8 km, max:7.6 km, mean: 6.0 km)
  :widths: 40 40 40 40
  :header-rows: 1

  * - pts number
    - CPU number
    - Wall time
    - CPU nodes
  * - 17004184
    - 64
    - 07:28:41
    - 4
  * - 17004184
    - 96
    - 06:38:11
    - 4
  * - 17004184
    - 128
    - 05:29:51
    - 6
  * - 17004184
    - 144
    - 04:59:49
    - 6
  * - 17004184
    - 168
    - 03:31:14
    - 7


Additional files
---------------------------

YML files
^^^^^^^^^^^^

A series of 4 input files are provided and can be used to run backward and forward models. The 3 backward models are divided in steps:

1. :code:`backward20Ma.yml` runs from 20 Ma to 15 Ma in *model time* but corresponds in *real time* in a backward model starting at 0 Ma and running back to 5 Ma ago;
2. :code:`backward15Ma.yml` goes from 15 Ma to 5 Ma in *model time*, *i.e.* 5 Ma ago to 15 Ma ago in *real time*;
3. :code:`backward5Ma.yml` runs for the remaining time from 5 Ma to 0 Ma in *model time* and 15 to 20 Ma ago in *real time*.

.. warning::
  When looking at these input files pay attention at the order of the tectonic meshes that force the model. It is also worth mentioning that for these backward models we do not compute the surface processes as a result we set :code:`fast` and :code:`backward` to **True**.


The last file (:code:`forward.yml`) is the forward model that needs to be ran at the end at it requires the outputs from the backward models.


Python3 scripts
^^^^^^^^^^^^^^^^

Two python scripts are provided to run these 4 files.

The first one can be used to run each of them individually (:code:`runModel.py`):

.. code:: python

  import argparse
  from gospl.model import Model as sim


  # Parsing command line arguments
  parser = argparse.ArgumentParser(
      description="This is a simple entry to run eSCAPE model.", add_help=True
  )
  parser.add_argument("-i", "--input", help="Input file name (YAML file)", required=True)
  parser.add_argument(
      "-v",
      "--verbose",
      help="True/false option for verbose",
      required=False,
      action="store_true",
      default=False,
  )
  parser.add_argument(
      "-l",
      "--log",
      help="True/false option for PETSC log",
      required=False,
      action="store_true",
      default=False,
  )

  args = parser.parse_args()

  # Reading input file
  model = sim(args.input, args.verbose, args.log)

  # Running forward model
  model.runProcesses()

  # Cleaning model
  model.destroy()

And can be run from the following command line in your terminal or the docker one...

.. code:: bash

  mpirun -np 4 python3 runModel.py -i backward5Ma.yml


The second (:code:`runBackwardForward.py`) automatise all the previous operations and will run all these models at once:

.. code:: python

  from mpi4py import MPI
  from gospl.model import Model as sim

  from scripts import mergeBack as merger

  MPIrank = MPI.COMM_WORLD.Get_rank()
  MPIcomm = MPI.COMM_WORLD

  forin = "forward.yml"
  backin = ["backward20Ma.yml", "backward15Ma.yml", "backward5Ma.yml"]
  backout = "output-backward"

  # Running the backwards models by periods
  for k in range(len(backin)):
    mod = sim(backin[k], False, False)
    mod.runProcesses()
    mod.destroy()
    if MPIrank == 0:
        print("", flush=True)

  # Merging all backward models into a single outputs
  if MPIrank == 0:
    merger.mergeBackModels(backin, backout)
    print("", flush=True)
  MPIcomm.Barrier()

  # Running the forward model forced with backward simulations
  mod = sim(forin, False, False)
  mod.runProcesses()
  mod.destroy()

And can be ran in a terminal using:

.. code:: bash

  mpirun -np 4 python3 runBackwardForward.py



Outputs & Paraview visualisation
--------------------------------

The model outputs are located in the output folder (:code:`dir` key as shown in the inputfile documentation) and consist of a time series file named :code:`gospl.xdmf` and 2 other folders (`h5` and `xmf`). The **XDMF** file is the main entry point for visualising the output and should be sufficient for most users.

The file can be opened with the `Paraview`_ software.


.. important::
  A series of Paraview states are provided and can be loaded in Paraview. When doing so ensure that you are setting the correct path to your model output folder!


.. _`gPlates`: http://portal.gplates.org
.. _`Scotese`: http://www.scotese.com
.. _`Paraview`: https://www.paraview.org/download/
.. _`YAML`: https://circleci.com/blog/what-is-yaml-a-beginner-s-guide/
.. _`Examples`: https://unisyd-my.sharepoint.com/:f:/g/personal/tristan_salles_sydney_edu_au/En8Wf56W_j9Jmqovx__PicgBczIcUogo6WuR-TVzZMHIMg?e=2pFtqT
