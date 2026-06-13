"""
Sediment provenance / source-to-sink attribution from goSPL output.

Standalone post-processor (Approach A in ``docs/DESIGN_PROVENANCE.md``): given a
per-vertex **source class** (from user source-rock polygons) and a goSPL run, it
tracks the eroded material down the per-step flow network and accumulates, in
each sink, the deposited volume attributable to each source rock type — per
basin and per pixel (mesh node) — plus the mean transport distance and an
optional copper-fertility layer.

The core (:class:`ProvenanceTracker`) is pure NumPy and operates on arrays, so it
runs without GIS / HDF5 dependencies; the optional I/O helpers
(:func:`classes_from_shapefile`, :func:`tracker_from_gospl`) import ``geopandas``
/ ``h5py`` lazily.

Method (deterministic continuum of particle tracking, mass-conserving per step):
for each output step the net ``erodep`` change gives the eroded volume (sources)
and deposited volume (sinks); the eroded material — recycled from the local
deposited pile, plus fresh bedrock of the node's source class — is routed
downstream and laid down where the model deposits, carrying its provenance and
travelled distance. Routing is **multiple-flow-direction** by default (flux
splits across all lower neighbours by slope, so a source can feed several
basins), with single steepest descent available (``routing='sfd'``).

Caveats (see the design doc): it sees the **net** erodep per *output* step (not
the compute dt), re-derives the routing from the saved elevation, and is a
reconstruction — exact volumes need the in-model tracer (Approach B). Provenance
is not prospectivity: the Cu layer is a fertility-weighted detritus proxy,
upstream of the ore-forming process.

Inputs (all per-vertex arrays are in the **global input-mesh** ordering — the
``npdata`` mesh ``v``; build per-vertex labels with
:func:`classes_from_shapefile` and reassemble model fields with
:class:`GosplOutput`):

==================  ========  ========  ===================================================
input               shape     required  meaning / how to obtain
==================  ========  ========  ===================================================
``coords``          (N, 2|3)  yes       global mesh vertex coordinates (mesh ``v``)
``cells``           (M, 3)    yes*      triangular cells (mesh ``c``); *or* pass ``neighbours``
``source_class``    (N,)      yes       integer source-rock class per vertex, in
                                        ``[0, n_classes)`` — point-in-polygon of the source
                                        shapefiles (assign an "other"/background class to
                                        vertices outside every polygon; do **not** use -1)
``n_classes``       scalar    yes       number of source classes
``erodep``/``elev`` (N,)      yes       per step, fed to :meth:`step` (from ``GosplOutput``)
``area``            (N,)      no        Voronoi cell area (m²); defaults to 1. Provide it for
                                        correct volumes on irregular meshes (ratios/%, which
                                        are area-independent, are unaffected)
``basin_id``        (N,)      no        sink-basin label per vertex, ``-1`` = not in a sink.
                                        **Optional** — needed *only* for the per-basin %
                                        roll-up (:meth:`basin_percentages`); per-pixel
                                        provenance is produced everywhere regardless
``cu_weight``       (C,)      no        copper fertility per source class; needed only for
                                        :meth:`cu_fraction`
==================  ========  ========  ===================================================
"""

import numpy as np


def build_neighbours(cells, npoints):
    """
    Build the vertex adjacency (CSR ``indptr``/``indices``) from triangular mesh
    cells — the undirected edge graph used for steepest-descent routing.

    :arg cells: (ntri, 3) integer vertex indices.
    :arg npoints: number of vertices.
    :return: (indptr, indices) CSR adjacency arrays.
    """
    cells = np.asarray(cells)
    e = np.vstack(
        [cells[:, [0, 1]], cells[:, [1, 2]], cells[:, [2, 0]]]
    )
    e = np.vstack([e, e[:, ::-1]])                       # undirected
    e = np.unique(e, axis=0)
    counts = np.bincount(e[:, 0], minlength=npoints)
    indptr = np.zeros(npoints + 1, dtype=np.int64)
    indptr[1:] = np.cumsum(counts)
    indices = e[:, 1].astype(np.int64)                   # e sorted by col 0 via unique
    return indptr, indices


def steepest_receivers(elev, indptr, indices):
    """
    Single-flow-direction (SFD) receivers: each node points to its lowest
    neighbour strictly below it, or to itself (a sink / outlet) if none is lower.

    :return: integer receiver array (length npoints).
    """
    n = len(elev)
    rcv = np.arange(n, dtype=np.int64)
    for i in range(n):
        lo = elev[i]
        best = i
        for k in range(indptr[i], indptr[i + 1]):
            j = indices[k]
            if elev[j] < lo:
                lo = elev[j]
                best = j
        rcv[i] = best
    return rcv


def downhill_edges(elev, coords, indptr, indices, routing="mfd", exponent=1.0):
    """
    Build the weighted downhill-flow graph as CSR edges
    ``(e_indptr, e_idx, e_weight, e_seg)``: for each node, the receivers it
    routes to, the fraction of its flux each takes, and the edge length.

    - ``routing='sfd'``: a single edge to the steepest-descent neighbour
      (weight 1); empty for sinks/outlets.
    - ``routing='mfd'``: an edge to **every** lower neighbour, weight
      :math:`\\propto \\text{slope}^{exponent}` and normalised to sum to one
      (multiple-flow-direction, as goSPL routes water/sediment). ``exponent``
      is the flow-partition exponent (1 ≈ uniform-to-slope spreading).
    """
    n = len(elev)
    elev = np.asarray(elev, dtype=np.float64)
    coords = np.asarray(coords, dtype=np.float64)
    # All adjacency edges as (src -> dst), vectorised (indices is CSR-grouped by
    # source, so `src` is the run-length expansion of the node ids).
    src = np.repeat(np.arange(n), np.diff(indptr))
    dst = indices
    dz = elev[src] - elev[dst]
    down = dz > 0.0                                    # keep downhill edges only
    src, dst, dz = src[down], dst[down], dz[down]
    seg = np.linalg.norm(coords[dst] - coords[src], axis=1)
    slope = np.divide(dz, seg, out=dz.copy(), where=seg > 0.0)
    w = slope ** exponent

    if routing == "sfd":
        # One edge per source: the steepest. Edges are grouped by source; pick
        # the max-slope edge in each group.
        order = np.lexsort((-slope, src))
        s_o = src[order]
        first = np.empty(len(s_o), dtype=bool)
        first[0] = True if len(s_o) else False
        first[1:] = s_o[1:] != s_o[:-1]
        keep = order[first]
        e_src, e_idx, e_seg = src[keep], dst[keep], seg[keep]
        e_w = np.ones(len(keep))
    else:
        # MFD: normalise weights per source to sum to one.
        wsum = np.bincount(src, weights=w, minlength=n)
        e_src, e_idx, e_seg = src, dst, seg
        e_w = w / wsum[src]

    # CSR indptr from per-source edge counts (edges stay grouped by source).
    counts = np.bincount(e_src, minlength=n)
    e_indptr = np.zeros(n + 1, dtype=np.int64)
    e_indptr[1:] = np.cumsum(counts)
    return (
        e_indptr,
        np.asarray(e_idx, dtype=np.int64),
        np.asarray(e_w, dtype=np.float64),
        np.asarray(e_seg, dtype=np.float64),
    )


def _sweep_impl(order, e_indptr, e_idx, e_w, e_seg, eroded, D):
    """
    Topological routing/deposition sweep (high → low): route the per-class
    eroded flux downstream, depositing the model's volume ``D`` at each node with
    the passing composition, and accumulate travel distance. Pure-arithmetic and
    Numba-``njit``-compatible (no Python objects), so the same source serves as
    both the reference and the compiled fast path. Returns the per-step
    increments ``(dep, dep_dist, dep_vol, exported)``.
    """
    n, nc = eroded.shape
    flux = eroded.copy()
    flux_dist = np.zeros(n)
    dep = np.zeros((n, nc))
    dep_dist = np.zeros(n)
    dep_vol = np.zeros(n)
    exported = np.zeros(nc)
    for idx in range(n):
        i = order[idx]
        tot = 0.0
        for c in range(nc):
            tot += flux[i, c]
        if tot <= 0.0:
            continue
        d = D[i] if D[i] < tot else tot                 # deposit min(D, flux)
        if d > 0.0:
            mean_here = flux_dist[i] / tot
            for c in range(nc):
                inc = flux[i, c] / tot * d
                dep[i, c] += inc
                flux[i, c] -= inc
            dep_dist[i] += mean_here * d
            dep_vol[i] += d
            flux_dist[i] -= mean_here * d
            tot -= d
        a = e_indptr[i]
        b = e_indptr[i + 1]
        if a == b:                                      # outlet / sink: export
            for c in range(nc):
                exported[c] += flux[i, c]
            continue
        fd_i = flux_dist[i]
        for k in range(a, b):
            j = e_idx[k]
            w = e_w[k]
            for c in range(nc):
                flux[j, c] += w * flux[i, c]
            flux_dist[j] += w * (fd_i + tot * e_seg[k])
        for c in range(nc):
            flux[i, c] = 0.0
        flux_dist[i] = 0.0
    return dep, dep_dist, dep_vol, exported


def _get_sweep(method):
    """Resolve the sweep kernel: Numba-compiled when available/requested, else
    the pure-Python reference."""
    if method == "python":
        return _sweep_impl
    try:
        import numba
    except ImportError:
        if method == "numba":
            raise ImportError("method='numba' needs numba (`pip install numba`)")
        return _sweep_impl                              # 'auto' fallback
    if not hasattr(_get_sweep, "_njit"):
        _get_sweep._njit = numba.njit(cache=True)(_sweep_impl)
    return _get_sweep._njit


class ProvenanceTracker:
    """
    Accumulate source-to-sink sediment provenance over a goSPL run.

    Feed cumulative ``erodep`` snapshots in time order with :meth:`step`; read
    the attribution with :meth:`basin_percentages`, :meth:`pixel_fractions`,
    :meth:`mean_distance` and :meth:`cu_fraction`.
    """

    def __init__(
        self,
        coords,
        source_class,
        n_classes,
        area=None,
        cells=None,
        neighbours=None,
        basin_id=None,
        cu_weight=None,
        routing="mfd",
        flow_exp=1.0,
        method="auto",
    ):
        """
        **Required**

        :arg coords: (npoints, 2|3) global-mesh vertex coordinates (mesh ``v``);
            used for travel distance and to define the vertex ordering all the
            other per-vertex arrays must follow.
        :arg source_class: (npoints,) integer source-rock class per vertex, each
            value in ``[0, n_classes)``. Obtain by point-in-polygon of the source
            shapefiles (:func:`classes_from_shapefile`); give vertices outside
            every polygon a dedicated "other"/background class — **not** -1.
        :arg n_classes: number of source classes.

        **One of** ``cells`` / ``neighbours`` (the routing connectivity)

        :arg cells: (ntri, 3) triangular mesh cells (mesh ``c``) — adjacency is
            built from these.
        :arg neighbours: pre-built ``(indptr, indices)`` CSR adjacency (instead
            of ``cells``).

        **Optional**

        :arg area: (npoints,) Voronoi cell area (m²); defaults to 1 (volume ==
            thickness). Provide it for correct volumes on irregular meshes;
            percentages/ratios are area-independent and unaffected.
        :arg basin_id: (npoints,) sink-basin label per vertex, ``-1`` = not in a
            sink basin. **Optional** — needed *only* for :meth:`basin_percentages`
            (the per-basin roll-up). Per-pixel provenance (:meth:`pixel_fractions`)
            is produced at every depositing node regardless.
        :arg cu_weight: (n_classes,) copper fertility per source class; needed
            only for :meth:`cu_fraction`.
        :arg routing: ``'mfd'`` (default, multiple-flow-direction — flux splits
            across all lower neighbours by slope, so a source can feed several
            basins) or ``'sfd'`` (single steepest descent).
        :arg flow_exp: MFD flow-partition exponent (slope power for the weights).
        :arg method: routing-sweep backend — ``'auto'`` (default: the
            Numba-compiled sweep when ``numba`` is installed, else pure Python),
            ``'numba'`` (require Numba), or ``'python'`` (reference). All give
            **identical** results (same topological sweep, including the
            ``min(D, flux)`` deposition cap). ``numba`` removes the Python
            per-node loop — the dominant cost on large meshes — for a
            several-fold (and growing with N) per-step speedup; the vectorised
            edge build / erosion is shared, so the overall gain is modest on
            small meshes. (A sparse transport-with-loss solve was evaluated and
            rejected — slower at realistic sizes due to LU fill-in, and only
            exact when no node is sediment-undersupplied; see
            ``docs/DESIGN_PROVENANCE.md``.)
        """
        self.coords = np.asarray(coords, dtype=np.float64)
        self.npoints = self.coords.shape[0]
        self.source_class = np.asarray(source_class, dtype=np.int64)
        self.n_classes = int(n_classes)
        if self.source_class.shape != (self.npoints,):
            raise ValueError(
                "source_class must have one value per vertex (%d), got %s"
                % (self.npoints, self.source_class.shape)
            )
        if self.source_class.min() < 0 or self.source_class.max() >= self.n_classes:
            raise ValueError(
                "source_class values must lie in [0, n_classes); assign a "
                "background/'other' class to vertices outside every source "
                "polygon (do not use -1)."
            )
        self.area = (
            np.ones(self.npoints) if area is None else np.asarray(area, dtype=np.float64)
        )
        if neighbours is not None:
            self.indptr, self.indices = neighbours
        elif cells is not None:
            self.indptr, self.indices = build_neighbours(cells, self.npoints)
        else:
            raise ValueError("provide either `cells` or `neighbours`")
        self.basin_id = None if basin_id is None else np.asarray(basin_id, dtype=np.int64)
        self.cu_weight = None if cu_weight is None else np.asarray(cu_weight, dtype=np.float64)
        if routing not in ("mfd", "sfd"):
            raise ValueError("routing must be 'mfd' or 'sfd'")
        if method not in ("auto", "numba", "python"):
            raise ValueError("method must be 'auto', 'numba' or 'python'")
        self.routing = routing
        self.flow_exp = float(flow_exp)
        self.method = method
        self._sweep = _get_sweep(method)

        # State.
        self.pile = np.zeros((self.npoints, self.n_classes))     # stored deposit comp.
        self.dep = np.zeros((self.npoints, self.n_classes))      # cumulative deposited
        self._dep_dist = np.zeros(self.npoints)                  # Σ dist·vol deposited
        self._dep_vol = np.zeros(self.npoints)                   # Σ vol deposited
        self.exported = np.zeros(self.n_classes)                 # left via outlets
        self._erodep_prev = None

    def step(self, elev, erodep):
        """
        Process one output step.

        :arg elev: (npoints,) surface elevation at this step (filled elevation
            preferred — ``fillFA``'s surface — so routing has no spurious pits).
        :arg erodep: (npoints,) cumulative erosion-deposition at this step.
        """
        erodep = np.asarray(erodep, dtype=np.float64)
        if self._erodep_prev is None:
            self._erodep_prev = erodep.copy()
            return
        vol = (erodep - self._erodep_prev) * self.area          # +dep / -ero (m^3)
        self._erodep_prev = erodep.copy()

        # 1. Erosion sources: recycle the local pile first, then fresh bedrock.
        E = np.maximum(-vol, 0.0)
        pile_tot = self.pile.sum(axis=1)
        take = np.minimum(E, pile_tot)
        rfrac = np.divide(take, pile_tot, out=np.zeros_like(take), where=pile_tot > 0)
        recycled = self.pile * rfrac[:, None]
        self.pile -= recycled
        eroded = recycled
        bed = E - take
        np.add.at(eroded, (np.arange(self.npoints), self.source_class), bed)

        # 2. Route downstream and 3. deposit where the model does.
        elev = np.asarray(elev, dtype=np.float64)
        e_indptr, e_idx, e_w, e_seg = downhill_edges(
            elev, self.coords, self.indptr, self.indices, self.routing, self.flow_exp
        )
        D = np.maximum(vol, 0.0)                                 # deposition (m^3)
        order = np.argsort(elev)[::-1].astype(np.int64)          # high -> low
        dep_inc, dist_inc, vol_inc, exp_inc = self._sweep(
            order, e_indptr, e_idx, e_w, e_seg, eroded, D
        )
        self.dep += dep_inc
        self.pile += dep_inc                                     # deposited -> local pile
        self._dep_dist += dist_inc
        self._dep_vol += vol_inc
        self.exported += exp_inc

    # ----- results -----

    def pixel_fractions(self):
        """Per-node deposited provenance fractions (npoints, n_classes); rows
        that received no deposit are all-zero."""
        tot = self.dep.sum(axis=1, keepdims=True)
        return np.divide(self.dep, tot, out=np.zeros_like(self.dep), where=tot > 0)

    def mean_distance(self):
        """Per-node flux-weighted mean source->deposition transport distance."""
        return np.divide(
            self._dep_dist, self._dep_vol,
            out=np.full(self.npoints, np.nan), where=self._dep_vol > 0,
        )

    def basin_percentages(self):
        """
        Dict ``{basin_id: array(n_classes)}`` of the % contribution of each
        source class to each sink basin. Requires ``basin_id``.
        """
        if self.basin_id is None:
            raise ValueError("basin_id was not provided")
        out = {}
        for b in np.unique(self.basin_id):
            if b < 0:
                continue
            v = self.dep[self.basin_id == b].sum(axis=0)
            s = v.sum()
            out[int(b)] = (100.0 * v / s) if s > 0 else np.zeros(self.n_classes)
        return out

    def cu_fraction(self):
        """Per-node Cu-sourced fraction of the deposit (needs ``cu_weight``)."""
        if self.cu_weight is None:
            raise ValueError("cu_weight was not provided")
        tot = self.dep.sum(axis=1)
        cu = self.dep @ self.cu_weight
        return np.divide(cu, tot, out=np.zeros_like(cu), where=tot > 0)


# --------------------------------------------------------------------------- #
# Optional I/O (lazy GIS / HDF5 imports)
# --------------------------------------------------------------------------- #

def classes_from_shapefile(
    coords, shapefile, attribute, background="other", crs=None
):
    """
    Assign each vertex an integer class by point-in-polygon against a source-rock
    (or sink-basin) polygon shapefile. The result is directly usable as the
    tracker's ``source_class``: every vertex gets a class in ``[0, n_classes)``,
    with vertices outside every polygon assigned a **background class** (so there
    are no ``-1`` sentinels). Lazily imports ``geopandas``/``shapely``.

    :arg coords: (npoints, 2|3) vertex coordinates in the shapefile's CRS
        (only the first two columns are used).
    :arg shapefile: path to a polygon shapefile.
    :arg attribute: polygon attribute holding the class label.
    :arg background: label for vertices outside every polygon (appended as its
        own class). Set to ``None`` to instead return ``-1`` for those vertices
        (e.g. when building a ``basin_id``, where -1 means "not a sink").
    :arg crs: optional CRS to reproject the polygons to before testing.
    :return: ``(class_array, {label: int} mapping)``. For ``source_class`` use
        the array as-is; for ``basin_id`` pass ``background=None``.
    """
    try:
        import geopandas as gpd
        from shapely.geometry import Point
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise ImportError(
            "classes_from_shapefile needs geopandas + shapely "
            "(`pip install geopandas`)"
        ) from exc

    gdf = gpd.read_file(shapefile)
    if crs is not None:
        gdf = gdf.to_crs(crs)
    pts = gpd.GeoDataFrame(
        geometry=[Point(xy) for xy in coords[:, :2]], crs=gdf.crs
    )
    joined = gpd.sjoin(pts, gdf[[attribute, "geometry"]], how="left", predicate="within")
    joined = joined[~joined.index.duplicated(keep="first")]
    labels = joined[attribute].to_numpy()
    uniq = [u for u in np.unique(labels[~_isnull(labels)])]
    mapping = {u: i for i, u in enumerate(uniq)}
    out = np.full(len(coords), -1, dtype=np.int64)
    for i, lab in enumerate(labels):
        if not _isnull(lab):
            out[i] = mapping[lab]
    if background is not None and (out < 0).any():
        bg = len(mapping)
        mapping[background] = bg
        out[out < 0] = bg
    return out, mapping


def _isnull(v):
    try:
        return bool(np.isnan(v))
    except (TypeError, ValueError):
        return v is None


class GosplOutput:
    """
    Reader that reassembles **global** per-vertex fields from goSPL's
    *partitioned* HDF5 output.

    goSPL writes one file per MPI partition per step
    (``<base>.<step>.p<rank>.h5``) holding only that rank's local nodes, with the
    partition meshes in ``topology.p<rank>.h5``. Provenance routing is global, so
    this reader maps each partition's local nodes onto the **global input mesh**
    ordering with a KDTree — the same local↔global map goSPL builds at load time
    (``tree.query(lcoords)``) — and scatters each field into a global array.
    Pass the global mesh coordinates (the ``npdata`` mesh ``v``) so the assembled
    arrays line up with your per-vertex ``source_class`` / ``basin_id``.
    """

    def __init__(self, h5dir, global_coords, file_base="gospl"):
        import glob
        import os

        import h5py
        from scipy.spatial import cKDTree

        self.h5dir = h5dir
        self.file_base = file_base
        self.npoints = len(global_coords)
        tfiles = sorted(glob.glob(os.path.join(h5dir, "topology.p*.h5")))
        if not tfiles:
            raise FileNotFoundError("no topology.p*.h5 in %s" % h5dir)
        tree = cKDTree(np.asarray(global_coords))
        self.ranks, self.maps = [], []
        for tf in tfiles:
            rank = int(os.path.basename(tf).split(".p")[1].split(".h5")[0])
            with h5py.File(tf, "r") as f:
                lc = np.asarray(f["coords"])
            _, gid = tree.query(lc)            # local node -> global vertex id
            self.ranks.append(rank)
            self.maps.append(gid.astype(np.int64))

    def field(self, step, name):
        """Assemble a global per-vertex array for ``name`` at output ``step``."""
        import os

        import h5py

        out = np.zeros(self.npoints)
        for rank, m in zip(self.ranks, self.maps):
            df = os.path.join(
                self.h5dir, "%s.%d.p%d.h5" % (self.file_base, step, rank)
            )
            with h5py.File(df, "r") as f:
                v = np.asarray(f[name])
            out[m] = v[:, 0] if v.ndim > 1 else v   # ghost overlaps agree
        return out


def _spec(s):
    """Split a ``file.npz:key`` argument into (path, key)."""
    path, _, key = s.rpartition(":")
    if not path:
        raise ValueError("expected file.npz:key, got %r" % s)
    return path, key


def main(argv=None):
    """
    CLI: attribute deposited sediment in sink basins to source rock types over a
    goSPL run, and write a per-basin CSV + a per-pixel ``.npz``.

    The global input mesh (``--mesh-npz``, the ``npdata`` file) defines the
    vertex ordering; goSPL's *partitioned* output (``--h5dir`` with
    ``<base>.<step>.p<rank>.h5`` + ``topology.p*.h5``) is reassembled onto it by
    :class:`GosplOutput`. Source classes and basins are per-vertex ``.npz``
    arrays in that same global ordering (rasterise your shapefiles with
    :func:`classes_from_shapefile`). Routing defaults to MFD.
    """
    import argparse

    p = argparse.ArgumentParser(description="goSPL sediment provenance attribution")
    p.add_argument("--mesh-npz", required=True, type=_spec,
                   help="global input mesh coords, as file.npz:key (e.g. mesh.npz:v)")
    p.add_argument("--cells", type=_spec, default=None,
                   help="global mesh cells, as file.npz:key (e.g. mesh.npz:c)")
    p.add_argument("--h5dir", required=True,
                   help="goSPL output h5/ directory (partitioned files)")
    p.add_argument("--file-base", default="gospl", help="output file base name")
    p.add_argument("--steps", required=True,
                   help="output steps, inclusive range 'a:b' or comma list")
    p.add_argument("--source", required=True, type=_spec,
                   help="per-vertex source-class array, file.npz:key (int)")
    p.add_argument("--n-classes", type=int, default=None)
    p.add_argument("--basins", type=_spec, default=None,
                   help="per-vertex sink-basin id array file.npz:key (-1 = none)")
    p.add_argument("--cu-weights", default=None,
                   help="comma-separated Cu fertility per class")
    p.add_argument("--routing", choices=["mfd", "sfd"], default="mfd")
    p.add_argument("--flow-exp", type=float, default=1.0)
    p.add_argument("--method", choices=["auto", "numba", "python"], default="auto",
                   help="routing-sweep backend (numba is ~50-100x faster on large meshes)")
    p.add_argument("--elev-field", default="elev",
                   help="HDF5 field used as the routing surface (e.g. waterFill)")
    p.add_argument("--out-prefix", required=True)
    args = p.parse_args(argv)

    mp, mk = args.mesh_npz
    coords = np.load(mp)[mk]
    cells = None
    if args.cells is not None:
        cp, ck = args.cells
        cells = np.load(cp)[ck]
    if ":" in args.steps and args.steps.count(":") == 1 and "," not in args.steps:
        a, b = (int(x) for x in args.steps.split(":"))
        steps = range(a, b + 1)
    else:
        steps = [int(x) for x in args.steps.split(",")]

    sp, sk = args.source
    source_class = np.load(sp)[sk].astype(np.int64)
    n_classes = args.n_classes or int(source_class.max() + 1)
    basin_id = None
    if args.basins is not None:
        bp, bk = args.basins
        basin_id = np.load(bp)[bk].astype(np.int64)
    cu_weight = None
    if args.cu_weights is not None:
        cu_weight = np.array([float(x) for x in args.cu_weights.split(",")])

    import h5py

    reader = GosplOutput(args.h5dir, coords, file_base=args.file_base)
    t = ProvenanceTracker(
        coords, source_class, n_classes, cells=cells,
        basin_id=basin_id, cu_weight=cu_weight,
        routing=args.routing, flow_exp=args.flow_exp, method=args.method,
    )

    # Per-step output: one HDF5 holding the cumulative provenance state after
    # each step (so the evolution is saved and can be viewed in ParaView), and a
    # per-basin time series CSV.
    steps = list(steps)
    basin_rows = []
    with h5py.File(args.out_prefix + ".h5", "w") as h:
        g = h.create_group("mesh")
        g["coords"] = coords
        if cells is not None:
            g["cells"] = cells
        h["steps"] = np.asarray(steps)
        for s in steps:
            t.step(reader.field(s, args.elev_field), reader.field(s, "erodep"))
            grp = h.create_group("step_%d" % s)
            frac = t.pixel_fractions().astype("float32")
            grp.create_dataset("fractions", data=frac, compression="gzip")
            grp.create_dataset(
                "dominant",
                data=np.where(frac.sum(1) > 0, frac.argmax(1), -1).astype("int32"),
                compression="gzip",
            )
            grp.create_dataset(
                "distance", data=t.mean_distance().astype("float32"), compression="gzip"
            )
            if cu_weight is not None:
                grp.create_dataset(
                    "cu_fraction", data=t.cu_fraction().astype("float32"),
                    compression="gzip",
                )
            if basin_id is not None:
                for b, pct in t.basin_percentages().items():
                    basin_rows.append((s, b, pct))
    print("wrote %s.h5 (%d steps)" % (args.out_prefix, len(steps)))

    if basin_id is not None:
        with open(args.out_prefix + "_basins.csv", "w") as fh:
            fh.write("step,basin," + ",".join("class%d" % c for c in range(n_classes)) + "\n")
            for s, b, pct in basin_rows:
                fh.write("%d,%d,%s\n" % (s, b, ",".join("%.3f" % v for v in pct)))
        print("wrote %s_basins.csv (per-step basin %%)" % args.out_prefix)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
