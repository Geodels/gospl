.. _provenance:

==============================================
Sediment provenance & source-to-sink attribution
==============================================

goSPL ships a post-processing tool (:mod:`gospl.analyse.provenance`) that
attributes the sediment **deposited in sink basins** to user-defined
**source-rock types**, per basin and per pixel (mesh node), with transport
distance and an optional **copper-fertility** layer. It answers questions such
as *"what fraction of the sediment in this basin came from each source
province, and by what pathway?"* — the source-to-sink ingredient for, e.g.,
detrital copper prospectivity.

The full design (data model, algorithm, caveats, and the planned in-model
variant) is in ``docs/DESIGN_PROVENANCE.md``.

Two approaches
--------------

#. **Standalone tracking** (implemented, described below): post-process an
   existing goSPL run — re-derive the flow network and route the eroded material
   to the sinks. No re-run; flexible; gives distances and pathways. Approximate
   (it sees the *net* erosion/deposition per output step).
#. **In-model provenance tracers** (see :ref:`below <prov-inmodel>`):
   generalise the dual-lithology two-fraction machinery to *N* provenance
   classes for exact, recycling-aware percentages — at the cost of re-running.

Standalone method
-----------------

For each output step the **net** ``erodep`` change gives the eroded volume
(sources) and the deposited volume (sinks). The eroded material — **recycled**
from the local deposited pile plus fresh bedrock of the node's source class — is
routed downstream and laid down where the model deposits, carrying its
provenance composition and travelled distance:

#. **Erode** — at nodes losing material, draw provenance from the stored deposit
   pile first (recycling cannibalised sediment) then from the bedrock source
   class.
#. **Route** — push the per-class flux down the **ice/water-style flow network**
   re-derived from the filled surface :math:`s`. Routing is
   **multiple-flow-direction** by default (flux splits across all lower
   neighbours by slope, so one source can feed several basins), with single
   steepest descent optional.
#. **Deposit** — capture the model's deposition volume at each node with the
   passing composition; sediment reaching the basins / sea settles there.

Conservation holds within each step. The routing is a topological sweep
(O(N+E)); it is provided as a single kernel used both as the pure-Python
reference and, when ``numba`` is installed, as a compiled fast path
(``method='auto'``).

Partitioned output
~~~~~~~~~~~~~~~~~~~

goSPL writes one HDF5 per MPI partition per step. Provenance routing is global,
so :class:`~gospl.analyse.provenance.GosplOutput` reassembles each field onto the
**global input-mesh** ordering by mapping every partition's local nodes to the
global vertices with a KDTree — the same local↔global map goSPL builds at load.
The per-vertex ``source_class`` / ``basin_id`` arrays are supplied in that global
ordering.

Inputs
------

============== =========  ====================================================
input          required   meaning
============== =========  ====================================================
mesh ``v``     yes        global vertex coordinates (defines the ordering)
mesh ``c``     yes        triangular cells (routing connectivity)
``source_class`` yes      source-rock class per vertex, in ``[0, n_classes)``
                          (point-in-polygon of the source shapefiles; outside
                          vertices get a background class — not -1)
goSPL output   yes        ``elev`` + ``erodep`` series (via ``GosplOutput``)
``area``       no         Voronoi cell area; default 1
``basin_id``   no         sink-basin label per vertex (-1 = not a sink) —
                          needed **only** for the per-basin roll-up
``cu_weight``  no         copper fertility per class — for the Cu layer
============== =========  ====================================================

Use :func:`~gospl.analyse.provenance.classes_from_shapefile` to build
``source_class`` / ``basin_id`` from polygon shapefiles (needs ``geopandas``).

Outputs
-------

Results are written **per output step** so the provenance evolution is preserved:

- ``<prefix>.h5`` — HDF5 with ``/mesh/{coords,cells}`` and a ``/step_<s>/`` group
  per step holding ``fractions`` (per-class composition), ``dominant`` (the
  leading source), ``distance`` (mean transport distance) and ``cu_fraction``;
- ``<prefix>.xdmf`` — an XDMF temporal wrapper so the HDF5 opens directly in
  **ParaView** as a time series (per-step dominant source, distance, Cu fraction
  and each source's ``frac_class<c>`` contribution);
- ``<prefix>_basins.csv`` — per-basin % through time.

Running it
----------

From a **terminal**::

    python -m gospl.analyse.provenance \
        --mesh-npz input/mesh.npz:v --cells input/mesh.npz:c \
        --h5dir output/h5 --steps 0:50 \
        --source input/source.npz:rock --basins input/sinks.npz:basin \
        --cu-weights 1,0,0.3 --routing mfd --method auto \
        --out-prefix provenance

From a **Jupyter notebook**, drive the API directly::

    from gospl.analyse.provenance import GosplOutput, ProvenanceTracker
    import numpy as np
    v = np.load("input/mesh.npz")["v"]; c = np.load("input/mesh.npz")["c"]
    src = np.load("input/source.npz")["rock"]
    reader = GosplOutput("output/h5", v)
    t = ProvenanceTracker(v, src, n_classes=int(src.max() + 1), cells=c,
                          basin_id=np.load("input/sinks.npz")["basin"])
    for s in range(51):
        t.step(reader.field(s, "elev"), reader.field(s, "erodep"))
    pct = t.basin_percentages()          # {basin: % per source class}
    frac = t.pixel_fractions()           # per-node composition

Performance
-----------

The flow-graph build is vectorised; the sequential routing sweep is a Numba-JIT
kernel when ``numba`` is available (``pip install "gospl[analysis]"``), several
times faster than pure Python on large meshes for an identical result.

Caveats
-------

- It is a **reconstruction** from the *net* erosion/deposition per *output* step,
  not the compute timestep — sub-step erode-then-redeposit is invisible. Exact
  volumes need the in-model tracer.
- Routing directions are re-derived from the saved filled elevation.
- **Provenance is not prospectivity:** the Cu layer is a fertility-weighted
  detritus proxy (which Cu-fertile sediment lands where), upstream of the
  ore-forming process (hydrodynamic concentration, facies, redox). Combine it
  with grain size (dual lithology) and depositional setting for a favourability
  layer.

.. _prov-inmodel:

In-model provenance tracers
---------------------------

For **conservation-exact, recycling-aware** provenance, goSPL can carry *N*
source-rock classes through its own erosion → transport → deposition →
stratigraphy, generalising the dual-lithology two-fraction machinery (the routed
sub-flux ``vSedF``, the per-layer composition ``stratHf``). Each stratigraphic
layer stores its source composition, so provenance percentages per pixel and per
layer are read directly off the saved stratigraphy — no post-processing
reconstruction. Provenance is a **passive label** (no erodibility / diffusivity /
porosity / sorting feedback), so it simply rides the total-sediment routing.

Enabled by a ``provenance:`` block (requires stratigraphy, ``time: strat``):

.. code:: yaml

    provenance:
        classes: 3                          # number of source-rock classes
        source: ['input/source', 'rock']    # per-vertex int class in [0, classes)
        # uniform: 0                         #   ... or a single class everywhere
        cu_weight: [1.0, 0.0, 0.3]           # optional copper fertility per class

How it works:

- **State** ``stratP[node, layer, class]`` — the per-class thickness in every
  layer (Σ over classes == the layer thickness ``stratH``), seeded so each
  initial layer and the bedrock carry the node's ``source`` class.
- **Erosion** (``erodeStrat``) splits the eroded sediment by the consumed layers'
  composition; eroded bedrock is attributed to the node's source class.
- **Transport** (``_getSedFlux``) routes each class sub-flux through the same
  upstream-integration operator as the total (linear, so the sub-fluxes sum to
  the total exactly).
- **Deposition** (``deposeStrat``) lays the deposit into the layer composition
  with the arriving (routed) fractions.
- **Advection & I/O** — ``stratP`` is advected with the pile under horizontal
  tectonics and written to / restored from the stratal HDF5 (the ``stratP``
  dataset).

**Conservation** is structural: ``stratP`` partitions ``stratH`` exactly for any
number of sources — with a single source every layer stays 100 % that class. The
**marine** sink (typically the dominant one) carries the basin-delivered source
mix (``_marineProvFraction``), so the recorded composition matches the eroded
supply to ~1e-6. The one remaining approximation — the per-class composition of
**intracontinental pit** deposits, which still use the through-flux composition
(overspill-chain mixing) — is **not** a conservation gap and is deferred; see
``docs/DESIGN_PROVENANCE.md`` §6 (phase B2b).

.. note::

    The **standalone** tool (above) and the **in-model** tracers answer the same
    question from opposite ends: the standalone post-processes an existing run
    (flexible, approximate, no re-run); the in-model tracers are exact and
    recycling-aware but require running with ``provenance:`` enabled.
