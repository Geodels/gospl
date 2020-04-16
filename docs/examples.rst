##################
Hand's on examples
##################

Input file
------------


Input files for *gospl* are based on `YAML`_ syntax.

A typical file will look like this:

```YAML
name: Global model from 20 Ma to present

domain:
    npdata: 'input8/elev20Ma'
    flowdir: 5
    fast: False
    backward: False
    interp: 1

time:
    start: -20000000.
    end: 0.
    tout: 1000000.
    dt: 1000000.
    tec: 1000000.

spl:
    K: 3.e-8
    Ff: 0.2
    wgth: 1.e-3

diffusion:
    hillslopeKa: 0.02
    hillslopeKm: 0.2
    sedimentK: 1000.

sea:
    position: 0.

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
  - start: -16000000.
    end: -15000000.
    mapH: 'input8/disp16Ma'

climate:
  - start: -20000000.
    map: ['input8/rain20Ma','r']
  - start: -15000000.
    map: ['input8/rain15Ma','r']

forcepaleo:
    dir: 'output-backward'
    steps: [5,10,5]

output:
    dir: 'forward'
    makedir: False

```

The YAML structure is shown through indentation (one or more spaces) and sequence items are denoted by a dash. At the moment the following component are available:

+ `domain`: definition of the unstructured grid containing the vtk grid `filename` and the associated field (here called `Z`) as well as the flow direction method to be used `flowdir` that takes an integer value between 1 (for SFD) and 12 (for Dinf) and the boundary conditions (`bc`: 'flat', 'fixed' or 'slope')
+ `time`: the simulation time parameters defined by `start`, `end`, `tout` (the output interval) and `dt` (the internal time-step).

Follows the optional forcing conditions:

+ `sea`: the sea-level declaration with the relative sea-level `position` (m) and the sea-level `curve` which is a file containing 2 columns (time and sea-level position).
+ `climatic` & `tectonic` have the same structure with a sequence of events defined by a starting time (`start`) and either a constant value (`uniform`) or a `map`.

Then the parameters for the surface processes to simulate:

+ `sp_br`: for the _stream power law_ with a unique parameter `Kbr` representing the The erodibility coefficient which is scale-dependent and its value depend on lithology and mean precipitation rate, channel width, flood frequency, channel hydraulics. It is worth noting that the coefficient _m_ and _n_ are fixed in this version and take the value 0.5 & 1 respectively.
+ `diffusion`: hillslope and marine diffusion coefficients. `hillslopeK` sets the _simple creep_ transport law which states that transport rate depends linearly on topographic gradient. The marine sediment are transported based on a diffusion coefficient `sedimentK`.

Finally, you will need to specify the output folder:

+ `output`: with `dir` the directory name and the option `makedir` that gives the possible to delete any existing output folder with the same name (if set to False) or to create a new folder with the give `dir` name plus a number at the end (e.g. outputDir_XX if set to True with XX the run number)

Additional informations
^^^^^^^^^^^^^^^^^^^^^^^

The tutorials proposed in the `Examples`_ repository propose several examples to create the required `vtk` input files (used in the input files for the domain, the climate, and the uplift) for running *gospl*.

In case of spatial change in forcing conditions, the user needs to specify (using the `map` element) for each vertices of the TIN, a series of fields for each forcing (tectonic displacements and precipitation) at chosen time intervals.

In cases where the forcing conditions are uniform over the entire domain one can chose to define its values in the YAML input file directly using the `uniform` element.

Sea level position can be set at a given relative position using the `position` element or using a `curve` which is a file containing 2 columns (time and sea-level position). In this later case, sea level position at any timestep is obtained by linear interpolation between closest sea level times.

Model-domain edge boundary conditions are set using the `bc` element in the YAML input file and the following options are available:

1. `flat`: edges elevations are set to the elevations of the closest non-edge vertices
2. `fixed`: edges elevations remain at the same position during the model run
3. `slope`: edges elevations are defined based on the closest non-edge vertices average slope


Outputs & Paraview visualisation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The model outputs are located in the output folder (`dir` element) and consist of a time series file named `gospl.xdmf` and 2 other folders (`h5` and `xmf`). The **XDMF** file is the main entry point for visualising the output and should be sufficient for most users.

The file can be opened with the `Paraview`_ software.

As shown in the video, after loading the file, we perform 2 operations:

1. first we use the `wrap by scalar` filter to create a 2D representation of the surface
2. then we define a `contour` line corresponding to the sea-level position


Available examples
------------------

`Examples`_


.. _`Paraview`: https://www.paraview.org/download/
.. _`YAML`: https://circleci.com/blog/what-is-yaml-a-beginner-s-guide/
.. _`Examples`: https://unisyd-my.sharepoint.com/:f:/g/personal/tristan_salles_sydney_edu_au/En8Wf56W_j9Jmqovx__PicgBczIcUogo6WuR-TVzZMHIMg?e=2pFtqT
