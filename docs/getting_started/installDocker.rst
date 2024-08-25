.. _installDocker:


=========================
Installation via Docker
=========================


Getting Docker application
--------------------------

Docker is an open platform for developing, shipping, and running applications. You can install it from:

+  `https://docs.docker.com/get-docker/ <https://docs.docker.com/get-docker/>`_

Once you have installed Docker Desktop for your operating system then enable the docker comand line (`Docker CLI`).


`Docker containers <https://hub.docker.com/r/geodels/gospl>`_ provide an easy-way to set up and distribute applications. They are also safe and consistent environments which
facilitate debugging and reproducibility of models / simulations.

goSPL on DockerHub
^^^^^^^^^^^^^^^^^^

The goSPL image contains all the dependencies and configuration files required to run models. Users can start running model as soon as they have downloaded the image, independently of the underlying operating system available on their machine.

.. important::
  
  The goSPL docker image is an easy-to-use approach for running our code especially if you have limited experience with Anaconda environments. 

Different flavors could be pulled using a tag and the recommended one is the ``gospl:latest`` that points to goSPL version ``v2024.09.01`` and uses the latest release of the code.

.. note::
  
  Depending on your operating system, you will be able to configure the docker application to set your resources: CPUs, memory, swap, or Disk image size.

On top of the libraries required for goSPL, additional ones have been included for pre- and post-processing.
Amongst then, you will find:

- xarray, rioxarray, xarray-spatial
- pygmt, cartopy, pyproj
- xesmf, maps, uxarray
- stripy, triangle,jigsawpy
- vtk, hdf5, netCDF4


Docker Dashboard
-------------------

To download and run goSPL from Docker, you might want to use the Docker **Dashboard** and define the Docker settings from there.

.. important::
  
  More specifically, you will need in the *Docker Settings Resources* to setup your CPU, Memory and Swap limits, and this is going to be important to ensure that your simulations do not run out of memory.

When starting your image for the first time, you will need to specify the container settings:

- **Ports**: Enter `0` to assign randomly generated host ports.
- **Volumes** *Host Path*: set a local directory that will be your workspace for running the simulations (it could be the folder containing the examples provided in this `repository <https://github.com/Geodels/goSPL-examples>`_). 
- **Volumes** *Container Path*: set it to `/notebooks`, this will be the folder in the Jupyter container that will be connected to the local directory specified above. 


Main command lines
-------------------

Alternatively, you can use the Terminal to pull and run the goSPL Docker image based on the command lines below.

Pulling the image
^^^^^^^^^^^^^^^^^

Once you have installed Docker on your system, you can ``pull`` the
`goSPL official image <https://hub.docker.com/u/geodels>`_ as follow::

  docker pull geodels/gospl:latest


You can list all the images available on your system as follow::

  docker images


An image can be deleted as follow::

  docker rmi geodels/gospl:latest


Starting the container from a terminal
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can then start a docker container (an instance of an image)::

  docker run -it -p 8888:8888 -d -v localDIR:/notebooks

where ``localDIR`` is the directory that contains the Jupyter Notebooks for your simulation or the examples provided in this `repository <https://github.com/Geodels/goSPL-examples>`_.

Once Docker is running, you could open the Jupyter notebooks on a web browser at the following address: `http://localhost:8888 <http://localhost:8888>`_. Going into the `/notebooks` folder you will access your ``localDIR`` directory.

You can list the containers currently existing on your machine by running::

  docker ps -a


The ``-a`` means "all container". The ``docker ps`` command only list
running containers.


Docker containers can be stop (so that they do not use CPU or RAM resource)::

  docker stop geodels/gospl:latest


They can also be deleted::

  docker rm geodels/gospl:latest


.. important::

  It's a good idea to keep track of how many containers have been created as
  they can rapidly take a lot of space on your machine.


To run goSPL from Docker, it is also recommended to use the terminal from the Jupyter interface. To activate the goSPL environment where all the libraries are installed you will have to run the following command::

  conda activate gospl


.. note::

  If you need additional libraries you could install them from the image Jupyter terminal by using either the `conda install` command or `pip install` command. You will need to first activate the conda environment in the terminal `conda activate gospl`.