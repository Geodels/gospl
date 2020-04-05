Installation
============

Docker_ installation
--------------------

Docker containers provide an easy-way to set up and distribute
applications. They also provide a safe and consistent environment which
facilitate debugging and reproducibility of models.

The **gospl** image contains all the dependencies and configuration files
required to run models. Users can start running model as soon as they
have downloaded the image, independently of the underlying operating system
available on on their machine.

It is strongly encouraged to run UWGeodynamics using the docker images
we provide on `Docker Hub`_

Different version of the `geodels/gospl` image can be pulled using a tag:

1. The *latest* tag points to the github master branch and uses the latest
   *gospl* release.
2. release tags such as *0.9.8* points to the specified version.

Command line
~~~~~~~~~~~~

Once you have installed docker on your system you can *pull* the
*gospl* official image as follow:

.. code:: bash

  docker pull geodels/gospl

You can list all the images available on your system as follow:

.. code:: bash

  docker images

An image can be deleted as follow:

.. code:: bash

  docker rmi geodels/gospl

You can then start a docker container. (An instance of
an image).

.. code:: bash

  docker run -d \
     --name my_container \
     -p 8888:8888 \
     -v my_vol:/live/share \
     geodels/gospl


You can access the container via your browser at the following
address: http://localhost:8888

It is also possible to ssh into the container as follow:

.. code:: bash

  # Bash
  docker run -it -v my_vol:/live/share \
     --entrypoint /bin/bash  \
     geodels/gospl


You can list the containers currently existing on your machine by running:

.. code:: bash

  docker ps -a

The "a" means "all container". The :code:`docker ps` command only list
running containers.

Docker containers can be stop (so that they do not use CPU or RAM resource):

.. code:: bash

  docker stop my_container

They can also be deleted:

.. code:: bash

  docker rm my_container

.. warning::

  It's a good idea to keep track of how many containers have been created as
  they can rapidly take a lot of space on your machine.

Kitematic_
~~~~~~~~~~

Kitematic_ is a program that provides a graphical user interface to
the *docker* daemon and to Docker Hub.

The software is available for Windows, MacOsx and Linux. Be aware that on
linux the installation may differ depending on the distribution you
are running.

1. Download and Install Kitematic_
2. Open Kitematic and search for the **gospl** image.
3. Create a container by clicking on the create button.

You should now have a container appearing on the left side of your
kitematic window. The first thing to do now is to create a link between
a local directory (A directory on your physical hard drive) and a volume
directory inside the docker container. A volume is a special directory
that can be accessed from outside the container. It is the location you
will use to save your results.

Local Installation
------------------

This is not recommended and involves installing *gospl* and all
its dependencies. Docker is highly recommended!!!

Requirements
~~~~~~~~~~~~

-  Python **>= 3.6**
-  Numpy **>= 1.18.2**
-  Scipy **>= 1.4.1**
-  Cython **>= 0.29.15**
-  mpi4py **== 3.0.3**
-  petsc4py **== 3.12.0**
-  h5py **== 2.10.0**
-  pandas **>= 0.24.2**
-  ruamel.yaml **== 0.16.10**
-  fastfunc **== 0.2.2**
-  meshio **== 4.0.10**
-  meshplex **== 0.12.3**
-  pre-commit **>= 1.21.0**


Install
~~~~~~~

**from Pip**

The **gospl** module can be installed directly from the Python
Package Index:

.. code:: bash

  pip3 install gospl

**from sources**

The module source files are available through github_

.. code:: bash

  git clone https://github.com/Geodels/gospl

It can then be installed locally on your system using

.. code:: bash

  export F90=gfortran
  python3 setup.py install --user

If you wish to uninstall **gospl** you can do:

.. code:: bash

  python3 setup.py install --record files.txt


.. role:: bash(code)
   :language: bash

To record a list of installed files in :bash:`files.txt`.
Once you want to uninstall you can use :bash:`xargs` to do the removal:


.. code:: bash

  xargs rm -rf < files.txt



HPC Run
----------------


After installation on HPC, you can submit a job to the queue system using for example:

.. code:: bash

   qsub job.pbs


Here is a minimal PBS script:

.. code:: bash

  #!/bin/bash

  # Project
  #PBS -P PROJECTNAME

  # 96 CPUs
  #PBS -l select=12:ncpus=8:mpiprocs=8:mem=60GB


  #PBS -l walltime=4:00:00
  #PBS -M email@address
  #PBS -m abe
  #PBS -q ALLOCNAME

  # set up environment
  module load gcc/4.9.3 python/3.6.5 petsc-gcc-mpich/3.11.1

  cd $PBS_O_WORKDIR

  # Launching the job!
  mpirun -np 96 python3 runMinimal.py

.. role:: python(code)
   :language: python

where the :python:`runMinimal.py` file will be of the form:


.. code:: python

  from gospl.model import Model as sim

  input = "imput.yml"

  # Minimal model
  model = sim(input, True, False)
  model.runProcesses()
  model.destroy()


That's it!!!

.. _Jupyter: http://jupyter.org/
.. _Docker: https://www.docker.com
.. _Docker Hub: https://hub.docker.com/repository/docker/geodels/gospl
.. _Kitematic: https://kitematic.com/
.. _github: https://github.com/Geodels/gospl
