=====================
HPC script example
=====================

Installation quick notes
--------------------------

Follows the :ref:`installation guide <installPypi>` for running :mod:`gospl` on HPC.

In most cases, ``MPI``, ``PETSc`` and ``Hdf5`` libraries for ``Python3`` will be installed.

Using ``pip`` you will be able to install the remaining libraries locally (``pip install XXXX --user``).


PBS script example
--------------------------

After installation on HPC, you might be able to submit a job to the queue system using for example:

.. code:: ipython

   qsub job.pbs

Here is a minimal ``PBS`` script provided as an example, it is worth noting that the environment will likely
depend on the HPC service you are using:

.. code:: ipython

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



A minimal python file for ``gospl``
------------------------------------

where the ``runMinimal.py`` file will be of the form:


.. code:: ipython

  from gospl.model import Model as sim

  input = "input.yml"

  # Minimal model
  model = sim(input, True, False)
  model.runProcesses()
  model.destroy()


Runtime
-------------------

Here are some *indicative* wall times for 2 specific resolution ran for a 10 million years simulation:

Case 1
^^^^^^^

Resolution min:7.6 km, max:10.3 km, mean: 9.1 km

.. csv-table::
    :header: "Pts number", "CPU number", "Wall time", "CPU nodes"
    :widths: 25, 20, 25, 15

    10612062, 8, 17:08:58, 1
    10612062, 16, 08:58:32, 2
    10612062, 32, 04:57:58, 4
    10612062, 64, 03:32:15, 6
    10612062, 128, 02:40:53, 7


Case 2
^^^^^^^

Resolution min:4.8 km, max:7.6 km, mean: 6.0 km

.. csv-table::
    :header: "Pts number", "CPU number", "Wall time", "CPU nodes"
    :widths: 25, 20, 25, 15

    17004184, 64, 07:28:41, 4
    17004184, 96, 06:38:11, 4
    17004184, 128, 05:29:51, 6
    17004184, 144, 04:59:49, 6
    17004184, 168, 03:31:14, 7
