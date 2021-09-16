.. _inputfile:

==============================
Input file parameters
==============================


.. note::

  Input files for  :mod:`gospl` are based on `YAML`_ syntax.
  The YAML structure is shown through indentation (one or more spaces) and sequence items are denoted by a dash. At the moment the following component are available:

.. role:: yaml(code)
   :language: yaml



:yaml:`domain`
---------------

.. raw:: html

   <div class="container">
   <div id="accordion" class="shadow tutorial-accordion">

       <div class="card tutorial-card">
           <div class="card-header collapsed card-link" data-toggle="collapse" data-target="#collapseOne">
               <div class="d-flex flex-row tutorial-card-header-1">
                   <div class="d-flex flex-row tutorial-card-header-2">
                       <button class="btn btn-dark btn-sm"></button>
                       Initial mesh definition and simulation declaration
                   </div>
               </div>
           </div>
           <div id="collapseOne" class="collapse" data-parent="#accordion">
               <div class="card-body">

**Declaration example**:

.. code:: ipython

 name: Global model from 20 Ma to present

 domain:
     npdata: 'input8/elev20Ma'
     flowdir: 5
     flowexp: 0.5
     fast: False
     backward: False
     interp: 1
     overlap: 1
     rstep: 25
     nperodep: 'strat8/erodep20Ma'
     npstrata: 'strat8/sed20Ma'

The following parameters are **required**:

a. the initial spherical surface mesh :yaml:`npdata` as well as
b. the flow direction method to be used :yaml:`flowdir` that takes an integer value between 1 (for SFD) and 6 (for MFD)
c. the exponent used in the flow direction approach. Default value is set to 0.42.

In addition the following optional parameters can be set:

d. the :yaml:`fast` key allows you to run a model without applying any surface processes on top. This is used to run backward model in a quick way, but can also potential be set to *True* if you want to check your input files prior to running a forward model with all options.
e. when running a backward model the :yaml:`backward` key has to be set to *True* as well!
f. the :yaml:`interp` key is set when running model with 3D cartesian displacements and allows you to choose the number of points that will be used when interpolating the spherical mesh after displacements. The key has 2 possible values: **1** or **3**. A value of **3** will take the 3 closest nodes to perform the interpolation and will tend to smooth the topography over time. A value of **1** will pick the closest point when performing the interpolation thus limiting the smoothing but potentially increasing the distorsion.
g. the :yaml:`overlap` key is set when running model with 3D cartesian displacements and specifies the number of ghost nodes used when defining the PETSc partition. It needs to be set so that all the points belonging to a single processors will not move further than the distances between the maximum horizontal displacement distance. The value will change depending of the resolution of your mesh.
h. to restart a simulation use the :yaml:`rstep` key and specify the time step number.
i. to start a simulation using a previous erosion/deposition map use the :yaml:`nperodep` key and specify a file containing for each vertex of the mesh the cumulative erosion deposition values in metres.
j. to start a simulation using an initial stratigraphic layer use the :yaml:`npstrata` key and specify a file containing for each vertex of the mesh the stratigraphic layer thickness, the percentage of fine lithology inside each layer and the porosities of the coarse and fine sediments (the multi-lithology option is only available for model without horizontal displacement and when the `backward` key is set to `False`).


.. raw:: html

               </div>
           </div>
       </div>
    </div>
    </div>


.. warning::

  It is worth noting that all the input files require to run a *gospl* simulation must be defined as numpy zip array (**.npz**). This allows to directly and efficiently load the dataset during initialisation. This is specially efficient when running large models.


:yaml:`time`
--------------------

.. raw:: html

   <div class="container">
   <div id="accordion" class="shadow tutorial-accordion">

       <div class="card tutorial-card">
           <div class="card-header collapsed card-link" data-toggle="collapse" data-target="#collapseTwo">
               <div class="d-flex flex-row tutorial-card-header-1">
                   <div class="d-flex flex-row tutorial-card-header-2">
                       <button class="btn btn-dark btn-sm"></button>
                       Setting model temporal evolution
                   </div>
               </div>
           </div>
           <div id="collapseTwo" class="collapse" data-parent="#accordion">
               <div class="card-body">

**Declaration example**:

.. code:: ipython

  time:
      start: -20000000.
      end: 0.
      tout: 1000000.
      dt: 250000.
      tec: 1000000.
      strat: 500000.

:yaml:`time` is also a required component of every input file. The following parameters are needed:

a. :yaml:`start` is the model start time in years,
b. :yaml:`end` is the model end time in years,
c. :yaml:`tout` is the output interval used to create model outputs,
d. :yaml:`dt` is the model internal time step (the approach in *gospl* uses an implicit time step.
e. :yaml:`tec` is the tectonic timestep interval used to update the tectonic meshes and perform the required horizontal displacements (vertical displacements are done every :yaml:`dt`).
f. :yaml:`strat` is the stratigraphic timestep interval used to update the stratigraphic record.

.. raw:: html

               </div>
           </div>
       </div>
    </div>
    </div>


.. important::

  In cases where the specify :yaml:`dt`, :yaml:`strat` and :yaml:`tec` parameters are greater than :yaml:`tout`, they will automatically be rescaled to match with the output interval. The :yaml:`tec` parameter should be set to similar to the temporal time step used in your reconstruction (usually around 1Ma). This time step is used to perform the horizontal displacements. The vertical displacements are updated for each time step. When turn-on the stratal records will be output at the same time as the output ones, but the file will potentially contain multiple stratigraphic layers per output if :yaml:`strat` is lower than :yaml:`tout`.


:yaml:`spl`
--------------------

.. raw:: html

   <div class="container">
   <div id="accordion" class="shadow tutorial-accordion">

       <div class="card tutorial-card">
           <div class="card-header collapsed card-link" data-toggle="collapse" data-target="#collapseThree">
               <div class="d-flex flex-row tutorial-card-header-1">
                   <div class="d-flex flex-row tutorial-card-header-2">
                       <button class="btn btn-dark btn-sm"></button>
                       Stream Power Law parameters
                   </div>
               </div>
           </div>
           <div id="collapseThree" class="collapse" data-parent="#accordion">
               <div class="card-body">

**Declaration example**:

.. code:: ipython

  spl:
      K: 3.e-8
      d: 0.42

This part of the input file define the parameters for the fluvial surface processes based on the *Stream Power Law* (SPL) and is composed of:

a. :yaml:`K` representing the erodibility coefficient which is scale-dependent and its value depend on lithology and mean precipitation rate, channel width, flood frequency, channel hydraulics. It is used in the SPL law: :math:`E = K (\bar{P}A)^m S^n`

.. warning::
  It is worth noting that the coefficient *m* and *n* are fixed in this version of *gospl* and take the value of *0.5* & *1* respectively.

b. Studies have shown that the physical strength of bedrock which varies with the degree of chemical weathering, increases systematically with local rainfall rate. Following `Murphy et al. (2016) <https://doi.org/10.1038/nature17449>`_, the stream power equation is adapted to explicitly incorporate the effect of local mean annual precipitation rate, P, on erodibility: :math:`E = (K_i P^d) (\bar{P}A)^m S^n`. :yaml:`d` (:math:`d` in the equation) is a positive exponent that has been estimated from field-based relationships to 0.42. Its default value is set to 0.


.. raw:: html

               </div>
           </div>
       </div>
    </div>
    </div>


:yaml:`diffusion`
----------------------


.. raw:: html

   <div class="container">
   <div id="accordion" class="shadow tutorial-accordion">

       <div class="card tutorial-card">
           <div class="card-header collapsed card-link" data-toggle="collapse" data-target="#collapseFour">
               <div class="d-flex flex-row tutorial-card-header-1">
                   <div class="d-flex flex-row tutorial-card-header-2">
                       <button class="btn btn-dark btn-sm"></button>
                       Hillslope and marine deposition parameters
                   </div>
               </div>
           </div>
           <div id="collapseFour" class="collapse" data-parent="#accordion">
               <div class="card-body">

**Declaration example**:

.. code:: ipython

  diffusion:
      hillslopeKa: 0.02
      hillslopeKm: 0.2
      clinSlp: 5.e-5
      smthS: 2.e5
      smthD: 1.e5
      offset: 500.
      nldep: False
      nlf: 1.e-3
      nlK: 3.e5
      nlKf: 5.e5
      nlKw: 7.e5


Hillslope processes in *gospl* is defined using a classical *diffusion law* in which sediment deposition and erosion depend on slopes (*simple creep*). The following parameters can be tuned based on your model resolution:

a. :yaml:`hillslopeKa` is the diffusion coefficient for the aerial domain,
b. :yaml:`hillslopeKm` is the diffusion coefficient for the marine domain,
c. :yaml:`clinSlp` is the maximum slope of clinoforms (needs to be positive), this slope is then used to estimate the top of the marine deposition based on distance to shore,
d. :yaml:`smthS` is the initial surface smoothing used to define the downstream transport of the marine sediments coming from rivers,
e. :yaml:`smthD` is the smoothing of the surface added to the freshly deposited sediments thicknesses used to define the downstream transport of the marine sediments coming from rivers
f. :yaml:`offset` is the offset in meters used to evaluate from the smoothed surface the maximum marine deposition thicknesses as sediments move on the continal slope and deep offshore basins.

.. warning::
  The following parameters are used to specify non-linear diffusion of rivers' sediments entering the ocean. This option is quite slow when not used on multi-processors and you might want to first look at the results of the simulation without this option turned on.

g. :yaml:`nldep` boolean set to *True* to account for non linear marine deposition,
h. :yaml:`nlf` nonlinear marine diffusion exponential factor for the freshly river deposited thicknesses (only accounted for if :yaml:`nldep` is True),
i. :yaml:`nlK` is the non linear diffusion coefficient for sediment deposited by rivers entering the marine environment (only accounted for if :yaml:`nldep` is True),
j. :yaml:`nlKf` is the diffusion coefficient for fine sediment deposited by rivers entering the marine environment. This parameter is only used when the multi-lithology and :yaml:`nlf` options are turned on,
k. :yaml:`nlKw` is the diffusion coefficient for weathered sediment deposited by hillslope processes and transported by rivers into the marine environment. This parameter is only used when the multi-lithology and :yaml:`nlf` options are turned on.

.. raw:: html

               </div>
           </div>
       </div>
    </div>
    </div>

:yaml:`sea`
--------------------

.. raw:: html

   <div class="container">
   <div id="accordion" class="shadow tutorial-accordion">

       <div class="card tutorial-card">
           <div class="card-header collapsed card-link" data-toggle="collapse" data-target="#collapseFive">
               <div class="d-flex flex-row tutorial-card-header-1">
                   <div class="d-flex flex-row tutorial-card-header-2">
                       <button class="btn btn-dark btn-sm"></button>
                       Sea-level (eustatic) forcing
                   </div>
               </div>
           </div>
           <div id="collapseFive" class="collapse" data-parent="#accordion">
               <div class="card-body">

**Declaration example**:

.. code:: ipython

  sea:
      position: 0.
      curve: 'data/sealevel.csv'


The sea-level declaration is defined with 2 optional parameters:

a. the relative sea-level :yaml:`position` in meters (optional),
b. a sea-level :yaml:`curve` *e.g.* a file containing 2 columns (time and sea-level position).

.. raw:: html

               </div>
           </div>
       </div>
    </div>
    </div>

:yaml:`tectonic`
----------------------

.. raw:: html

   <div class="container">
   <div id="accordion" class="shadow tutorial-accordion">

       <div class="card tutorial-card">
           <div class="card-header collapsed card-link" data-toggle="collapse" data-target="#collapseSix">
               <div class="d-flex flex-row tutorial-card-header-1">
                   <div class="d-flex flex-row tutorial-card-header-2">
                       <button class="btn btn-dark btn-sm"></button>
                       Tectonic forcing parameters
                   </div>
               </div>
           </div>
           <div id="collapseSix" class="collapse" data-parent="#accordion">
               <div class="card-body">

**Declaration example**:

.. code:: ipython

  tectonic:
    - start: -20000000.
      end: -19000000.
      mapH: 'input8/disp20Ma'
    - start: -19000000.
      end: -18000000.
      mapH: 'input8/disp19Ma'
    - start: -18000000.
      end: -17000000.
      mapH: 'input8/disp18Ma'
    - start: -17000000.
      end: -16000000.
      mapH: 'input8/disp17Ma'
      mapV: 'input8/dispv17Ma'
    - start: -16000000.
      end: -15000000.
      mapV: 'input8/dispv16Ma'

Follows the tectonic forcing conditions with a sequence of events defined by a starting time (:yaml:`start`) and either a vertical only forcing (*e.g.* uplift and/or subsidence defined with :yaml:`mapV`) or a fully 3D displacement mesh :yaml:`mapH`. These displacements are set in metres per year.

.. raw:: html

               </div>
           </div>
       </div>
    </div>
    </div>

.. important::

  As mentioned above and for the next key parameter as well, these forcing files are defined as numpy zip array (**.npz**).


:yaml:`compaction`
--------------------


.. raw:: html

   <div class="container">
   <div id="accordion" class="shadow tutorial-accordion">

       <div class="card tutorial-card">
           <div class="card-header collapsed card-link" data-toggle="collapse" data-target="#collapseSeven">
               <div class="d-flex flex-row tutorial-card-header-1">
                   <div class="d-flex flex-row tutorial-card-header-2">
                       <button class="btn btn-dark btn-sm"></button>
                       Compaction & porosity variables defintion
                   </div>
               </div>
           </div>
           <div id="collapseSeven" class="collapse" data-parent="#accordion">
               <div class="card-body">

**Declaration example**:

.. code:: ipython

  compaction:
      phis: 0.49
      phif: 0.63
      phiw: 0.65
      z0s: 3700.0
      z0f: 1960.0
      z0w: 1580.0

The compaction module is turned-on when a multi-lithology model is ran (_i.e._ the :yaml:`npstrata` key is defined). We assume  different depth-porosity relationships for the 3 considered lithology types, the following parameters are required:

a. lithology one (coarser lithology) porosity at the surface :yaml:`phis`,
b. lithology two (finer lithology) porosity at the surface :yaml:`phif`,
c. lithology three (weathered lithology) porosity at the surface :yaml:`phiw`,
d. e-folding depth :yaml:`z0s` of lithology one (in metres)
e. e-folding depth :yaml:`z0f` of lithology two (in metres)
f. e-folding depth :yaml:`z0w` of lithology three (in metres)

.. raw:: html

               </div>
           </div>
       </div>
    </div>
    </div>


:yaml:`climate`
--------------------

.. raw:: html

   <div class="container">
   <div id="accordion" class="shadow tutorial-accordion">

       <div class="card tutorial-card">
           <div class="card-header collapsed card-link" data-toggle="collapse" data-target="#collapseEight">
               <div class="d-flex flex-row tutorial-card-header-1">
                   <div class="d-flex flex-row tutorial-card-header-2">
                       <button class="btn btn-dark btn-sm"></button>
                       Climatic (rainfall) forcing conditions
                   </div>
               </div>
           </div>
           <div id="collapseEight" class="collapse" data-parent="#accordion">
               <div class="card-body">

**Declaration example**:

.. code:: ipython

  climate:
    - start: -20000000.
      map: ['input8/rain20Ma','r']
    - start: -15000000.
      uniform: 1.


The climatic forcing is defined in a similar fashion as the tectonic one with again a sequence of events by a starting time (:yaml:`start`) and either an uniform rainfall over the entire mesh (:yaml:`uniform`) or with a precipitation mesh :yaml:`map`. The rainfall values have to be in metres per year.


.. raw:: html

               </div>
           </div>
       </div>
    </div>
    </div>


:yaml:`forcepaleo`
-----------------------


.. raw:: html

   <div class="container">
   <div id="accordion" class="shadow tutorial-accordion">

       <div class="card tutorial-card">
           <div class="card-header collapsed card-link" data-toggle="collapse" data-target="#collapseNine">
               <div class="d-flex flex-row tutorial-card-header-1">
                   <div class="d-flex flex-row tutorial-card-header-2">
                       <button class="btn btn-dark btn-sm"></button>
                       Forcing paleo-topography definition
                   </div>
               </div>
           </div>
           <div id="collapseNine" class="collapse" data-parent="#accordion">
               <div class="card-body">

**Declaration example**:

.. code:: ipython

  forcepaleo:
      dir: 'output-backward'
      steps: [5,10,5]

For simulations that require to be forced with paleo-topography maps obtained from backward models, you will also have to set this key composed of 2 parameters:

a. :yaml:`dir` the directory containing the outputs of the backward model,
b. :yaml:`steps` the steps from the model outputs that will be used to force the forward model topography.

.. important::

  The :yaml:`steps` often correspond to the time where you have a paleotopography dataset that you want to match for example from a Scotese paleotopography map.

.. raw:: html

               </div>
           </div>
       </div>
    </div>
    </div>

:yaml:`output`
--------------------

.. raw:: html

   <div class="container">
   <div id="accordion" class="shadow tutorial-accordion">

       <div class="card tutorial-card">
           <div class="card-header collapsed card-link" data-toggle="collapse" data-target="#collapseTen">
               <div class="d-flex flex-row tutorial-card-header-1">
                   <div class="d-flex flex-row tutorial-card-header-2">
                       <button class="btn btn-dark btn-sm"></button>
                       Output folder definition
                   </div>
               </div>
           </div>
           <div id="collapseTen" class="collapse" data-parent="#accordion">
               <div class="card-body">

**Declaration example**:

.. code:: ipython

  output:
      dir: 'forward'
      makedir: False

Finally, you will need to specify the output folder, with 2 possible parameters:

a. :yaml:`dir` gives the output directory name and
b. the option :yaml:`makedir` gives the ability to delete any existing output folder with the same name (if set to False) or to create a new folder with the given `dir` name plus a number at the end (*e.g.* outputDir_XX if set to True with XX the run number). It allows you to avoid overwriting on top of previous runs.


.. raw:: html

               </div>
           </div>
       </div>
    </div>
    </div>


.. _`Paraview`: https://www.paraview.org/download/
.. _`YAML`: https://circleci.com/blog/what-is-yaml-a-beginner-s-guide/
