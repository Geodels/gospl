.. _outputs:

====================================
Output variables & visualisation
====================================

goSPL writes its results to the directory given by the ``output: dir`` key of
the input file. Each output interval (``time: tout``) produces, per MPI rank, an
**HDF5** data file plus a small **XDMF** descriptor; a top-level ``.xmf`` /
``.xdmf`` file ties the per-step, per-rank pieces together. **Open the top-level
``gospl.xdmf`` in ParaView** (or any XDMF reader) — it exposes every surface
field below as a node-centred scalar that can be coloured, warped and filtered.

This page lists the fields goSPL can write and what they mean, then shows how to
turn the flat-mesh or global-sphere output into a 3-D surface in ParaView.

Output fields
-------------

Many fields are only written when the corresponding process is enabled (noted in
the *When written* column). Unless stated otherwise, lengths are in **metres**
and rates are **per year**; discharges/fluxes are **volumetric (m³/yr)**.

Core surface
^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 18 52 30

   * - Field
     - Meaning
     - When written
   * - ``elev``
     - Surface elevation of the bed.
     - always
   * - ``erodep``
     - **Cumulative** erosion (negative) and deposition (positive) since the start of the run.
     - always
   * - ``EDrate``
     - Erosion–deposition **rate** over the current step (m/yr; negative = incision).
     - always

Water & drainage
^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 18 52 30

   * - Field
     - Meaning
     - When written
   * - ``FA``
     - Water **flow accumulation** (drainage discharge, m³/yr) on the actual surface.
     - always
   * - ``fillFA``
     - Flow accumulation on the **depression-filled** surface — drainage routed through lakes/pits.
     - always
   * - ``waterFill``
     - Filled water-surface elevation (lake / depression water level).
     - when depressions are filled
   * - ``sedLoad``
     - River **sediment load** carried downstream (m³/yr).
     - always
   * - ``sedLoadF``
     - Fine-fraction sediment **load**/flux (dual-lithology runs only); coarse flux = ``sedLoad − sedLoadF``.
     - dual lithology
   * - ``surfFineFrac``
     - Surface exposed **fine fraction** (mud share, 0–1, of the topmost stratigraphic layer) — the in-place composition, complementing the ``sedLoadF`` flux.
     - dual lithology
   * - ``rain``
     - Precipitation rate forcing the run (m/yr).
     - when rainfall is set
   * - ``evap``
     - Evaporation rate (m/yr).
     - when evaporation is set

Soil, tectonics & flexure
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 18 52 30

   * - Field
     - Meaning
     - When written
   * - ``soilH``
     - Soil / regolith thickness.
     - ``soil`` enabled
   * - ``uplift``
     - Vertical tectonic displacement rate applied this step (m/yr; positive = uplift).
     - vertical tectonics
   * - ``flexIso``
     - Flexural isostatic vertical deflection (m; positive = up).
     - ``flexure`` enabled

Glacial (ice) fields
^^^^^^^^^^^^^^^^^^^^^

Written only when an ``ice`` section is present (see :ref:`surfproc`). The
diagnostic glacial model derives these each step from the routed ELA
accumulation — there is no ice-thickness time integration.

.. list-table::
   :header-rows: 1
   :widths: 18 62 20

   * - Field
     - Meaning
     - Units
   * - ``iceH``
     - **Ice thickness** from the Bahr discharge scaling :math:`H=\mathrm{eheight}\cdot\mathrm{fwidth}\cdot Q^{0.3}`. Zero below the terminus.
     - m
   * - ``iceFA``
     - **Ice discharge** :math:`Q` — the volume of ice routed through each cell (the glacial analogue of river ``FA``).
     - m³/yr
   * - ``iceUb``
     - **Basal sliding velocity** magnitude from Glen's sliding law on the ice thickness and bed slope (:math:`u_b\propto H^{n-1}|\nabla s|^{n-1}\nabla s`). This is the driver of abrasion.
     - m/yr
   * - ``iceAbr``
     - **Glacial (vertical) abrasion rate** :math:`E_g = K_g\,|u_b|^{l}` — bed lowering by sliding ice. (Lateral valley-wall erosion, ``Kl``, contributes to bed change via the till budget but is *not* added into this field.)
     - m/yr
   * - ``iceMelt``
     - **Glacial meltwater delivered to the rivers** — the ice accumulation released as liquid water where the ice melts out (ablation zone / terminus). With ``melt_conserve: True`` it is discharge-conserving (total melt = total accumulation), so it shows where glaciers feed streamflow. *Not* the same as the internal till-deposition weight.
     - m³/yr

.. note::

   To see **glacial valley deepening**, colour by ``erodep`` / ``EDrate`` (the
   bed change) — ``iceAbr`` reports only the *vertical* abrasion rate. To see
   **valley widening** (U-shaping), enable lateral wall erosion with
   ``abrasion: Kl > 0`` (see :ref:`surfproc`) and inspect ``erodep`` across the
   valley: vertical ``Kg`` deepens the trough, lateral ``Kl`` lowers the flanking
   walls. Increasing ``Kg`` deepens *more* but does **not** widen — widening is
   controlled by ``Kl``.

Stratigraphy
^^^^^^^^^^^^

Written when stratigraphic recording is on (``strat`` interval set), in the
separate ``h5/stratal.<step>.p<rank>.h5`` files. These are **per-layer** arrays
(one column per recorded layer): ``stratZ`` (layer elevation), ``stratH`` (layer
thickness), ``phiS`` (porosity), ``stratK`` (erodibility multiplier), ``stratP``
(per-layer provenance), and, under dual lithology, ``stratHf`` (fine-fraction
thickness) and ``phiF`` (fine porosity). The per-layer fine fraction is
``stratHf / stratH``.

The stratal HDF5 files have no XDMF wrapper of their own. Several
post-processing tools turn them into figures (all documented, with the exact
commands, in :ref:`Running goSPL <running>`):

* **3-D volume for ParaView** — ``gospl-strata-volume`` stacks the layers into a
  wedge (triangular-prism) volume coloured by lithology or provenance.
* **Publication sections / wells / Wheeler diagrams** — ``gospl-section`` draws
  vector cross-sections (along x/y/a path), horizontal depth slices, synthetic
  wells (one, or several side by side via ``well_panel``) and
  chronostratigraphic (Wheeler) charts, coloured by deposition elevation,
  thickness, lithology, porosity, age, provenance or **depositional facies**
  (water depth at deposition). Distance is shown in km, time in ky; see
  :ref:`Running goSPL <running>` for the full option list.
* **Surface composition over time** — instead of the pile, colour the main
  ``gospl.xdmf`` by ``surfFineFrac`` (in-place mud share) or ``sedLoadF`` (fine
  transport flux).

Visualising in ParaView
-----------------------

The mesh is stored flat — node coordinates plus the ``elev`` (``Z``) scalar — so
by default ParaView shows a flat sheet (planar model) or a smooth sphere (global
model) coloured by whatever field you pick. Apply one of the filters below to
turn elevation into relief.

Flat / planar model — *Warp By Scalar*
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Load the ``.xdmf`` and colour by a field (e.g. ``elev``).
2. Apply **Filters → Alphabetical → Warp By Scalar**.
3. Set **Scalars = Z** (the elevation), choose a **Scale Factor** (1 for true
   scale, larger for vertical exaggeration), and **Apply**.

The sheet now stands up as topography; re-colour by ``erodep``, ``FA``,
``iceH``, … as needed.

Global / spherical model — *Calculator* radial warp
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For a global run the nodes already lie on a sphere of radius ~6 378 137 m, so a
*Warp By Scalar* would displace them incorrectly. Instead push each node
**radially** by its (exaggerated) elevation with a **Calculator** filter:

1. Apply **Filters → Alphabetical → Calculator**.
2. Tick **Coordinate Results** and **Result Normals**.
3. Set the expression to::

       (iHat*coordsX + jHat*coordsY + kHat*coordsZ) * (1 + Z*50/6378137)

   Here ``coordsX/Y/Z`` are the node coordinates, ``Z`` is the elevation scalar,
   and ``50`` is the vertical-exaggeration factor (raise/lower to taste; ``1``
   gives true scale). The term ``1 + Z·50/6378137`` scales the radius by the
   (exaggerated) elevation as a fraction of the planet radius.
4. **Apply**, then in the **Properties** set the **Interpolation** to **Flat**
   (rather than *Gouraud*) so faceted relief reads clearly.

.. note::

   ``6378137`` is **Earth's radius in metres** — goSPL's default. If you are
   modelling **another planet**, set the planet radius (and surface gravity, used
   by the flexural-isostasy solver) in the ``domain`` block of the input file and
   substitute that radius into the formula above. For example, a Mars run
   (radius ≈ 3 389 500 m, gravity ≈ 3.71 m/s²)::

       domain:
           npdata: ['input/mesh', 'v', 'c', 'z']   # spherical mesh: vertices, cells, elevation
           radius: 3389500.0                       # planet radius (m); default 6378137 (Earth)
           gravity: 3.71                           # surface gravity (m/s^2); default 9.81 (Earth)
           flowdir: 6                              # MFD flow directions

   would use ``... * (1 + Z*50/3389500)`` in the Calculator. (The ``bc`` edge
   key is for flat/planar meshes only — a spherical planet mesh has no edges.)

Re-colour the warped surface by any output field to inspect topography,
drainage, ice or sediment in 3-D.

Gridded export for PyGMT / ArcGIS (basins, χ, river profiles)
-------------------------------------------------------------

goSPL's native output is an unstructured triangular mesh, which PyGMT and ArcGIS
do not read directly. The ``gospl-grid`` tool reassembles the global mesh,
interpolates every surface field of a step onto a **regular CF-NetCDF grid**
(lon/lat for global meshes), and runs a raster D8 hydrology pass to add
``drainage_area``, ``basin`` and ``chi`` (:math:`\chi`); a companion notebook
API extracts per-basin **river longitudinal profiles**. The base level for
catchments/χ defaults to the run's sea level. See :ref:`Running goSPL <running>`
for the command, options and the river-profile API.
