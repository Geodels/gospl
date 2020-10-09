=====================
HPC script example
=====================


Follows the installation guide for running :mod:`gospl` on HPC.

In most cases, ``MPI``, ``PETSc`` and ``Hdf5`` libraries for ``Python3`` will be installed.

Using ``pip`` you will be able to install the remaining libraries locally (``pip install XXXX --user``).

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


where the ``runMinimal.py`` file will be of the form:


.. code:: ipython

  from gospl.model import Model as sim

  input = "input.yml"

  # Minimal model
  model = sim(input, True, False)
  model.runProcesses()
  model.destroy()
