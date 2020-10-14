.. _advance:

.. note::

  The documentation provided below assumes that you already have a data folder containing relevant files for your simulation:

  - paleo topography maps
  - paleo precipitations grids
  - horizontal displacements files


==============================================
Remeshing workflows
==============================================


The workflows described below are available from the `Github repository <https://github.com/Geodels/gospl/tree/master/notebooks>`_ but you can also directly download them as a **tar** file from `here <https://drive.google.com/file/d/1SvRj27NBF4aA2E8svyniysQtDHuiVNIf/view?usp=sharing>`_ or using the following command in your terminal::

    gdown https://drive.google.com/uc?id=1SvRj27NBF4aA2E8svyniysQtDHuiVNIf
    tar xvf notebooks.tar


1. Initial unstructured mesh generation
---------------------------------------

We start by creating the unstructured mesh using ``JIGSAW`` library. It is ran from the terminal using ``buildMesh.py`` script::

    python3 script/buildMesh.py -t=20 -d=data -s=100,30,15


This script takes 3 arguments:

- ``t`` the time interval in Ma of the starting simulation time (here 20 Ma),
- ``d`` the folder containing the input conditions. The code assumes that the paleo-elevations are located in the folder **data** under **paleomap** and are `netCDF` files of the form ``XXMa.nc`` with ``XX`` the specified time. It also assumes that the displacement maps are located in **velocity** and are off the form ``velocity_ XXMa.xy``. Lastly for the paleo-precipitation, the assumed file is under **precipitation** and are ``netCDF`` files of the form ``XXMa.nc`` as well.
- ``s`` is the space conditions for the jigsaw algorithm and consists of 3 values: the spacing in km for the mesh in the deep ocean (<=-1000 m), the spacing in km across shelf margin (>=-1000 m and < 0m) and the spacing in km in the continental domain.

Before going further, you can check the mesh validity by loading the created ``VTK`` file (``inputXX/meshXX.vtk``) in ``Paraview`` or by using ``pyvista`` library as in the Jupyter notebooks examples.


2. Creating backward grids
---------------------------------------

To constrain the surface evolution over time and impose *uplift/subsidence maps*, we use a series of backward grids that are compared with obtained elevations from the surface process model. The strategy consists first in taking the next available paleomap and then applying to it the horizontal displacements backward.

Considering two consecutive maps, for example a paleo-elevations at 20 Ma and one at 15 Ma, we use the following Jupyter notebook (``backElev.ipynb``) to build backward elevations from 15 Ma to 20 Ma. These regular maps of 0.1 degree resolution are stored in a Numpy compressed file under the data folder.


3. Running gospl over 1 Ma
---------------------------------------

The script above will generate the input conditions required to run the surface process model over a 1 Ma time scale. You will need to specify in the ``YAML`` input file:

- ``npdata: inputXX/XXMa`` where ``XX`` is the time specified above (set to 20 above). This file contains the mesh information (coordinates, neighbours, cells, elevations)
- ``map: ['inputXX/rainXXMa','r']`` which is the paleo-precipitation mapped on the irregular mesh.

We do not use the displacement map yet. The model is ran with the following command::

    mpirun -np NB python3 script/runModel.py -i model20Ma.yml


*NB* is the number of processors to use and the script required the input file for :mod:`gospl` as an argument.

4. Define vertical displacements
---------------------------------------

We will now take our surface processes result after 1 Ma simulation and compare it with the corresponding backward elevation calculated in **2**.

.. important::

  The difference in elevations will be used to define a vertical displacements file for the corresponding period. Assuming the validity of paleo-elevation maps, this vertical displacement should already account for the effect of dynamic topography, flexural isostasy responses and other tectonic forcing.

However we do not apply the total computed displacements in one go as the difference reflects the evolution over 5 Ma. Therefore we scale the calculated differences over time with a factor <=1. For example, we could use a linear approach, that implies a linear change over time between 2 intervals and use a factor of 0.2 for 20 to 19Ma, 0.4 for 19 to 18Ma, 0.6 for 18 to 17Ma, 0.8 for 17 to 16Ma and 1.0 for the last interval. This vertical displacement factor can easily be changed manually both temporally and spatially to reflect non-linear tectonic histories.

To create the vertical tectonic file, we use the ``vertDisp.ipynb`` Jupyter notebook. The displacements are in *m/yr* and stored in a Numpy compressed file under the `inputXX` folder and are named `vtecXXMa.npz`.

5. Running tectonically constrained gospl over 1 Ma
-----------------------------------------------------------

We now rerun the :mod:`gospl` model but this time including the vertical displacements imposed based on the backward model differences. To do so we add the following in the ``YAML`` input file:


.. code:: ipython

  tectonic:
    - start: -20000000.
      end: -19000000.
      mapV: 'input20/vtec20Ma'

Obviously this will need to be modified according to the simulated time interval. A good thing to do is also to modify the output folder file to keep previous simulation on disk. Once modified we run the surface process model using the previous command::

    mpirun -np NB python3 script/runModel.py -i model20Ma.yml


6. Perform horizontal displacements and remeshing
-----------------------------------------------------------

We now apply the horizontal displacements on the final surface process output from the last run. We will also extract the required input file for the following run from 19 Ma to 18 Ma in our case. This is done by using the following command line::

    python3 script/npzMesh.py -t=19 -d=data -s=100,30,15 -i=model20Ma.yml -n=100 -a=1 -r=20


Where the arguments ``t``, ``d`` and ``s`` are the same as in step 1. In addition, the following arguments are required:

- ``i`` the ``YAML`` input file from the previous simulation
- ``n`` the final time step number from :mod:`gospl` model output
- ``a`` the applied displacement time interval in Ma (here set to 1 Ma for example)
- ``r`` the paleo-precipitation file time step to use (see ``d`` in **1** for some explanations), paleo-precipitation is supposed uniformed between 2 increments (we have values at 20 & 15 Ma in our case)

This command will create 3 compressed Numpy files that are stored in the ``inputXX`` folder where ``XX`` is the value provided with the ``t`` argument. The elevation is given by ``inputXX/XXMa.npz``, the erosion deposition values are in the file ``inputXX/erodepXXMa.npz``, and the rainfall in ``inputXX/rainXXMa.npz``. These 3 files are then specified in the next ``YAML`` input file:


.. code:: ipython

    domain:
        npdata: 'input19/19Ma'
        flowdir: 5
        fast: False
        backward: False
        interp: 1
        npvalue: 'input19/erodep19Ma'

    climate:
      - start: -19000000.
        map: ['input19/rain19Ma','r']

With next input file created, steps **3** to **6** are iteratively repeated to simulate surface evolution over time.


7. Visualisation
-----------------------------------------------------------

To visualise the output over time in ``Paraview`` one need to merge all successive :mod:`gospl` outputs together. This is done by using the Jupyter notebook ``combXDMF.ipynb``.
