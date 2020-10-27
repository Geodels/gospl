.. _installDocker:

=========================
Installation via Docker
=========================

`Docker containers <https://hub.docker.com/r/geodels/gospl>`_ provide an easy-way to set up and distribute
applications. They are also a safe and consistent environment which
facilitate debugging and reproducibility of models.

The :mod:`gospl` image contains all the dependencies and configuration files
required to run models. Users can start running model as soon as they
have downloaded the image, independently of the underlying operating system
available on on their machine.

It is strongly encouraged to run :mod:`gospl` using the docker image
we provide on `Docker containers <https://hub.docker.com/r/geodels/gospl>`_.

Different version of the image can be pulled using a tag:

1. The ``latest`` tag points to the github master branch and uses the latest
   release.
2. release tags such as ``v0.2.0`` points to the specified version.


Depending on your operating system, you will be able to configure the docker
application to set your resources: CPUs, memory, swap, or Disk image size.


Main command lines
-------------------

Pulling the image
^^^^^^^^^^^^^^^^^

Once you have installed Docker on your system, you can ``pull`` the
`gospl official image <https://hub.docker.com/r/geodels/gospl>`_ as follow::

  docker pull geodels/gospl:v0.2.0


You can list all the images available on your system as follow::

  docker images


An image can be deleted as follow::

  docker rmi geodels/gospl


Starting the container from a terminal
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can then start a docker container (an instance of
an image)::

  docker run -d \
     --name my_container \
     -p 8888:8888 \
     -v my_vol:/live/share \
     geodels/gospl


You can access the container via your browser at the following
address: `http://localhost:8888 <http://localhost:8888>`_

It is also possible to ``ssh`` into the container as follow::

  # Bash
  docker run -it -v my_vol:/live/share \
     --entrypoint /bin/bash  \
     geodels/gospl:v0.2.0


You can list the containers currently existing on your machine by running::

  docker ps -a


The ``-a`` means "all container". The ``docker ps`` command only list
running containers.


Docker containers can be stop (so that they do not use CPU or RAM resource)::

  docker stop my_container


They can also be deleted::

  docker rm my_container


.. note::

  It's a good idea to keep track of how many containers have been created as
  they can rapidly take a lot of space on your machine.


Kitematic
------------

`Kitematic <https://kitematic.com/>`_ is a program that provides a graphical user interface to
the ``docker`` daemon and to Docker Hub.

Downloading a tagged image
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The software is available for Windows, Mac OS X and Linux. Be aware that on
linux the installation may differ depending on the distribution you
are running.

1. Download and Install Kitematic_
2. Open Kitematic and search for the tag corresponding to a specific the ``gospl`` version.
3. Download the container by clicking on the create button.

Linking a local volume
^^^^^^^^^^^^^^^^^^^^^^^

You should now have a container appearing on the left side of your
kitematic window.

The first thing to do now is to create a link between
a local directory (A directory on your physical hard drive) and a volume
directory inside the docker container. A ``volume`` is a special directory
that can be accessed from outside the container. It is the location you
will use to save your own results.
