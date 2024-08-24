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

Different flavors could be pulled using a tag and recommended ones are:

1. The ``gospl:2024.09.01`` tag points to the ``2024.09.01`` goSPL branch,
2. The ``gosplcontainer`` tag points to the github ``v2023`` goSPL branch,
3. The ``gosplmaster`` tag points to the github master branch and uses the latest
   release.

.. note::
  
  Depending on your operating system, you will be able to configure the docker application to set your resources: CPUs, memory, swap, or Disk image size.

On top of the libraries required for goSPL, the following main ones are available:

- xarray, rioxarray, xarray-spatial
- pygmt, cartopy, pygplates
- xesmf, metpy
- stripy, triangle
- vtk, hdf5, netCDF4


Main command lines
-------------------

Pulling the image
^^^^^^^^^^^^^^^^^

Once you have installed Docker on your system, you can ``pull`` the
`goSPL official image <https://hub.docker.com/u/geodels>`_ as follow::

  docker pull geodels/gospl:2024.09.01


You can list all the images available on your system as follow::

  docker images


An image can be deleted as follow::

  docker rmi geodels/gospl:2024.09.01


Starting the container from a terminal
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can then start a docker container (an instance of
an image)::

  docker run -it -p 8888:8888 -d -v localDIR:/notebooks

where ``localDIR`` is the directory that contains the Jupyter Notebooks for your simulation.

Once Docker is running, you could open the Jupyter notebooks on a web browser at the following address: `http://localhost:8888 <http://localhost:8888>`_. Going into the `/notebooks` folder you will access your ``localDIR`` directory.

You can list the containers currently existing on your machine by running::

  docker ps -a


The ``-a`` means "all container". The ``docker ps`` command only list
running containers.


Docker containers can be stop (so that they do not use CPU or RAM resource)::

  docker stop geodels/gospl:2024.09.01


They can also be deleted::

  docker rm geodels/gospl:2024.09.01


.. important::

  It's a good idea to keep track of how many containers have been created as
  they can rapidly take a lot of space on your machine.


.. note::

  If you need additional libraries you could install them from the Jupyter terminal directly (when opening the terminal change the shell to `bash` by using running `bash`) within the container by using either the `conda install` command or `pip install` command.
