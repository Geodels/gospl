"""
Post-processing: rasterise a goSPL surface to a **regular CF-NetCDF grid** for
PyGMT / ArcGIS, with drainage **basins**, **chi** (:math:`\\chi`) and drainage
area computed on the grid — plus per-basin **river longitudinal profiles**.

goSPL writes an unstructured triangular mesh (XDMF + per-partition HDF5) that
PyGMT and ArcGIS do not read directly. This tool reassembles the global mesh,
interpolates every surface field of a chosen output step onto a regular grid
(honouring the mesh triangulation), runs a standard raster D8 hydrology pass on
the gridded elevation, and writes one NetCDF holding every field plus
``drainage_area``, ``basin`` (integer outlet id) and ``chi``.

Raster D8 hydrology (self-contained, no heavy deps): a priority-flood + epsilon
fill guarantees drainage, then steepest-descent receivers give the drainage
area, the basin each cell drains to, the flow distance to the outlet, and
:math:`\\chi = \\int (A_0/A)^{m/n}\\,dl` integrated upstream from base level.
These are **raster/D8** quantities (the convention for chi and basins), so they
can differ slightly from goSPL's internal MFD drainage.

The companion API extracts and plots river profiles for a chosen basin:

* :func:`basin_rivers` — the channel network (cells with drainage area above a
  threshold) of one basin, split into the **main stem** and its **tributaries**,
  each with along-channel distance, elevation, chi and area.
* :func:`plot_long_profile` — longitudinal profile (distance or chi vs
  elevation), main stem highlighted.
* :func:`plot_basin_map` — the channel network drawn on the surface.

Options
-------
``--h5dir DIR``
    goSPL output ``h5`` directory (**required**).
``--mesh FILE:vkey:ckey``
    global mesh ``.npz`` and its vertex / cell keys (**required**), e.g.
    ``input/mesh.npz:v:c`` — supplies the global coordinates the partitions map
    onto and the triangulation used for interpolation.
``--out FILE``
    output NetCDF path (default ``surface.nc``).
``--step N``
    output step to grid (default: the last one found).
``--spacing DX[,DY]``
    grid spacing in mesh units (default: the median node spacing).
``--fields a,b,...``
    surface fields to grid (default: all present in the step).
``--mn FLOAT`` / ``--a0 FLOAT``
    chi concavity ``m/n`` (default 0.5) and reference area (default 1.0).
``--base-level FLOAT``
    elevation defining the coast/outlet — catchments drain to it and chi is
    measured from it (default: the run's **sea level**, read from the step's
    ``.xmf``; falls back to 0). Cells at/below it are marine and excluded.
``--latlim FLOAT``
    geographic meshes only: crop ``|latitude|`` to this limit, dropping the
    singular polar caps (default 89.9).
``--file-base BASE``
    mesh-output file base name (default ``gospl``).

A **global (spherical)** mesh is auto-detected and gridded in **lon/lat**
(degrees, CF ``lon``/``lat``). A near-full-longitude global grid is treated as
**periodic in longitude**, so the D8 flow and the rasterisation wrap across the
antimeridian — a continent/river crossing ±180° stays intact, regardless of
where the seam falls. Hydrology metrics are latitude-weighted
(``dx = R·cos(lat)·dλ``). The **poles** are a genuine singularity of any lon/lat
grid (meridians converge); they are treated as boundaries and the polar caps are
cropped (``--latlim``). Drainage that genuinely crosses a pole needs a polar
projection (out of scope) — regional/continental basins are unaffected.

Runnable as ``gospl-grid`` (installed) or ``python -m gospl.analyse.gridexport``::

    gospl-grid --h5dir myrun/h5 --mesh input/mesh.npz:v:c --out surface.nc \\
        --spacing 1000 --mn 0.5

    # then, in a notebook:
    from gospl.analyse.gridexport import grid_export, basin_rivers, \\
        plot_long_profile, plot_basin_map
    g = grid_export("myrun/h5", "input/mesh.npz", step=20, spacing=1000.)
    riv = basin_rivers(g, basin_id=g["main_basin"], area_threshold=5e6)
    plot_long_profile(riv); plot_basin_map(g, riv)
"""

import os
import heapq
import argparse

import numpy as np


# ---------------------------------------------------------------------------
# Reading / rasterising
# ---------------------------------------------------------------------------

def _list_fields(h5dir, file_base, step):
    """Surface-field dataset names present in the step's mesh HDF5 (rank 0)."""
    import glob
    import h5py

    f0 = sorted(glob.glob(os.path.join(h5dir, "%s.%d.p*.h5" % (file_base, step))))
    if not f0:
        raise FileNotFoundError(
            "no %s.%d.p*.h5 in %s" % (file_base, step, h5dir)
        )
    skip = {"coords", "cells"}
    with h5py.File(f0[0], "r") as f:
        return [k for k in f.keys() if k not in skip and f[k].ndim <= 2]


def _read_sealevel(h5dir, file_base, step):
    """
    Read the model sea level for a step from its companion ``xmf`` file (goSPL
    encodes it in the ``sea`` attribute's Function constant). Returns the float,
    or ``None`` if the file is absent / unparseable.
    """
    import re

    xmf = os.path.join(
        os.path.dirname(os.path.normpath(h5dir)), "xmf",
        "%s%d.xmf" % (file_base, step),
    )
    try:
        with open(xmf) as f:
            txt = f.read()
    except OSError:
        return None
    m = re.search(r'Function="\$0 \* [0-9.eE+-]+ \+ (-?[0-9.eE+]+)"', txt)
    return float(m.group(1)) if m else None


def _read_time(h5dir, file_base, step):
    """
    Read the model display time (years) for a step from its companion ``xmf``
    (goSPL writes ``<Time Value="..."/>``). Returns the float, or ``None`` if the
    file is absent / unparseable.
    """
    import re

    xmf = os.path.join(
        os.path.dirname(os.path.normpath(h5dir)), "xmf",
        "%s%d.xmf" % (file_base, step),
    )
    try:
        with open(xmf) as f:
            txt = f.read()
    except OSError:
        return None
    m = re.search(r'<Time\s+Value="(-?[0-9.eE+]+)"', txt)
    return float(m.group(1)) if m else None


# Seeds tried per grid point when the nearest mesh node's incident triangles do
# not contain it. Only ever applied to the residual of the first pass, so a
# generous value costs almost nothing; 6 sufficed for the randomised Delaunay
# regression, 16 leaves margin for a strongly refined/anisotropic mesh.
_LOCATE_FALLBACK_K = 16


def _seam_mask(lon_of_cells):
    """Triangles whose longitude span exceeds 180° (cross the framing seam)."""
    return (lon_of_cells.max(axis=1) - lon_of_cells.min(axis=1)) > 180.0


def _node_tri_csr(cells, npts):
    """
    Node -> incident-triangle adjacency in CSR form: the triangles touching node
    ``i`` are ``tids[ptr[i]:ptr[i + 1]]``.

    Built with an UNSTABLE sort on purpose: only the *set* of triangles per node
    matters, never their order, and dropping ``kind="stable"`` is ~4x faster on a
    multi-million-triangle mesh (12.1 s -> 3.2 s for 11.8M triangles).
    """
    flat = np.asarray(cells).ravel()
    tids = np.repeat(np.arange(len(cells), dtype=np.int64), 3)
    order = np.argsort(flat)
    deg = np.bincount(flat, minlength=npts)
    ptr = np.zeros(npts + 1, dtype=np.int64)
    np.cumsum(deg, out=ptr[1:])
    return ptr, tids[order], int(deg.max()) if deg.size else 0


def _nearest_mesh_node(tree, gx, gy, R=None):
    """
    Nearest mesh node for every grid point, used to seed the point location.

    For a geographic mesh the query is done in **3-D on the sphere**, so the
    result is independent of where the longitude seam is put — which is what
    lets both longitude framings of a periodic grid share one query.

    :return: ``(nn, refine)`` where ``nn`` is the nearest node per grid point and
        ``refine(idx, k)`` returns the ``k`` nearest nodes for the grid-point
        subset ``idx`` (used only for the residual of the first location pass, so
        the tree is never re-queried over the whole grid).
    """
    if R is not None:
        la = np.radians(gy.ravel())
        lo = np.radians(gx.ravel())
        cl = np.cos(la)
        pts = np.column_stack([R * cl * np.cos(lo), R * cl * np.sin(lo),
                               R * np.sin(la)])
    else:
        pts = np.column_stack([gx.ravel(), gy.ravel()])

    def _query(p, k=1):
        try:
            return tree.query(p, k=k, workers=-1)[1]
        except TypeError:                   # scipy < 1.6 has no `workers`
            return tree.query(p, k=k)[1]

    nn = _query(pts)

    def refine(idx, k):
        return np.atleast_2d(_query(pts[idx], k=k))

    return nn, refine


def _build_tri_interp(x, y, cells, gx, gy, seed, csr, tri_mask=None):
    """
    Precompute a reusable linear (barycentric) interpolation operator from the
    mesh TIN onto a regular grid, and return a closure ``interp(values) -> 2-D
    grid``.

    The expensive part of rasterising — locating every grid point in its
    containing triangle — depends only on the **geometry**, not on the field
    values, so it is done once here and each subsequent field becomes a cheap
    gather + weighted sum.

    Point location does NOT use ``matplotlib.tri``. Its
    ``TrapezoidMapTriFinder`` is a general-purpose structure for arbitrary
    triangulations and building it costs O(minutes) on a multi-million-triangle
    mesh: measured on a 5.9M-node / 11.8M-triangle global mesh gridded at 0.1°,
    ``get_trifinder()`` alone took **158 s**, and a periodic grid needs TWO of
    them (one per longitude framing) for a total of ~320 s per step before any
    field is touched. We do not need that generality: a goSPL mesh is a
    near-uniform Delaunay triangulation, so the triangle containing a grid point
    is incident to that point's nearest mesh NODE. So seed with one KDTree query
    (the tree is already built for the partition reassembly) and test only the
    handful of triangles incident to that node — the same location, ~10x
    cheaper: 174 s -> 16 s, of which the shared parts (adjacency + KDTree query)
    are paid once for both framings.

    Verified against ``TrapezoidMapTriFinder`` on that mesh: identical located
    set (6 468 566 of 6 480 000 grid points) and **100.000000 % identical
    containing triangle**, zero misses in either direction. A ``k``-nearest
    fallback still runs for any point the first pass leaves unresolved, so an
    anisotropic mesh degrades in cost rather than in correctness.

    Points outside the mesh / in masked triangles are ``NaN``, as before.

    :arg seed: ``(nn, refine)`` from ``_nearest_mesh_node`` — the nearest mesh
        node per grid point plus a k-nearest query for the residual
    :arg csr: ``(ptr, tids, maxdeg)`` node -> triangle adjacency over ALL
        triangles (shared across framings; ``tri_mask`` is applied per framing
        here rather than baked into the adjacency)
    """
    nn, refine = seed
    ptr, tids, maxdeg = csr
    shape = gx.shape
    pxv = np.ascontiguousarray(gx.ravel(), dtype=np.float64)
    pyv = np.ascontiguousarray(gy.ravel(), dtype=np.float64)
    npix = pxv.size

    tmask = (np.zeros(len(cells), dtype=bool) if tri_mask is None
             else np.asarray(tri_mask, dtype=bool))

    found = np.full(npix, -1, dtype=np.int64)
    w0 = np.zeros(npix)
    w1 = np.zeros(npix)
    w2 = np.zeros(npix)
    pending = np.ones(npix, dtype=bool)
    eps = -1.0e-12

    def _test(idx, tidx):
        """Barycentric test of grid points `idx` against candidate `tidx`."""
        tv = cells[tidx]
        i0, i1, i2 = tv[:, 0], tv[:, 1], tv[:, 2]
        x0, y0 = x[i0], y[i0]
        x1, y1 = x[i1], y[i1]
        x2, y2 = x[i2], y[i2]
        px, py = pxv[idx], pyv[idx]
        det = (y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2)
        good = det != 0.0                            # skip degenerate triangles
        inv = np.where(good, 1.0 / np.where(good, det, 1.0), 0.0)
        a = ((y1 - y2) * (px - x2) + (x2 - x1) * (py - y2)) * inv
        b = ((y2 - y0) * (px - x2) + (x0 - x2) * (py - y2)) * inv
        c = 1.0 - a - b
        ok = good & ~tmask[tidx] & (a >= eps) & (b >= eps) & (c >= eps)
        hit = idx[ok]
        found[hit] = tidx[ok]
        w0[hit] = a[ok]
        w1[hit] = b[ok]
        w2[hit] = c[ok]
        pending[hit] = False

    def _sweep(seed):
        """Try every triangle incident to the seed node of each pending point."""
        for j in range(maxdeg):
            idx = np.flatnonzero(pending)
            if idx.size == 0:
                return
            nd = seed[idx]
            slot = ptr[nd] + j
            keep = slot < ptr[nd + 1]
            if not keep.any():
                continue
            _test(idx[keep], tids[slot[keep]])

    _sweep(nn)

    # Fallback for anything the nearest node did not resolve, widening the seed
    # to the _LOCATE_FALLBACK_K nearest nodes. Most of the residual is genuinely
    # outside the mesh (or under a masked seam stripe) and stays NaN, so this
    # pass cannot be triggered selectively -- but it runs ONLY on the residual
    # (~0.2 % of the grid on a real global run), so it is near-free.
    #
    # It is not decorative. On a near-uniform Delaunay mesh -- what goSPL builds
    # -- the containing triangle is always incident to the nearest node. On an
    # IRREGULAR triangulation it need not be: a sliver triangle can contain a
    # point comfortably in its interior (barycentric 0.007 / 0.473 / 0.520 in the
    # randomised case that motivated this) while not touching any of the 4
    # nearest nodes. Without the widened retry that point would silently become a
    # hole in the raster.
    idx = np.flatnonzero(pending)
    if idx.size:
        knn = refine(idx, min(_LOCATE_FALLBACK_K, len(x)))
        for m in range(1, knn.shape[1]):
            still = pending[idx]
            if not still.any():
                break
            # default to nn so an unexpected gap can never index ptr[-1]
            seed_m = nn.copy()
            seed_m[idx[still]] = knn[still, m]
            _sweep(seed_m)

    return _finish(found, w0, w1, w2, cells, pxv, shape)


def _finish(found, w0, w1, w2, cells, pxv, shape):
    """Freeze a located grid into an `interp(values) -> 2-D grid` closure."""
    valid = found >= 0
    tv = cells[found[valid]]
    i0, i1, i2 = tv[:, 0], tv[:, 1], tv[:, 2]
    a, b, c = w0[valid], w1[valid], w2[valid]

    def interp(values):
        v = np.asarray(values)
        out = np.full(pxv.shape, np.nan, dtype=np.float64)
        out[valid] = a * v[i0] + b * v[i1] + c * v[i2]
        return out.reshape(shape)

    return interp


def grid_export(h5dir, mesh, step=None, vkey="v", ckey="c", spacing=None,
                nx=None, ny=None, fields=None, mn=0.5, a0=1.0,
                base_level=None, file_base="gospl", latlim=None):
    """
    Build the regular-grid surface (fields + D8 hydrology) for one step.

    :return: a dict of 2-D grids (``ny, nx``) keyed by field name plus
        ``drainage_area``, ``basin``, ``chi``, ``flowdist``, the 1-D axes
        ``x``/``y``, ``mask``, ``receiver`` (flat index, -1 at outlets),
        ``spacing`` ``(dx, dy)`` and ``main_basin`` (largest basin id).
    """
    import glob
    import h5py
    from scipy.spatial import cKDTree

    data = np.load(mesh)
    coords = np.asarray(data[vkey], dtype=np.float64)
    cells = np.asarray(data[ckey], dtype=np.int64)

    # Geographic (spherical) mesh? Grid in lon/lat; otherwise planar x/y.
    # A near-full-longitude geographic mesh is treated as PERIODIC so flow wraps
    # across the antimeridian (continents crossing the seam stay intact).
    r = np.linalg.norm(coords, axis=1)
    geographic = bool(r.mean() > 1.0e5 and (r.std() / max(r.mean(), 1.0) < 1.0e-3))
    if geographic:
        R = float(r.mean())
        x = np.degrees(np.arctan2(coords[:, 1], coords[:, 0]))            # lon
        y = np.degrees(np.arcsin(np.clip(coords[:, 2] / R, -1.0, 1.0)))   # lat
    else:
        R = None
        x, y = coords[:, 0], coords[:, 1]

    if step is None:
        step = max(
            int(os.path.basename(f).split(".")[1])
            for f in glob.glob(os.path.join(h5dir, "%s.*.p*.h5" % file_base))
        )

    # Reassemble each global field by mapping partition nodes onto the global
    # mesh with a KDTree (the same local<->global map goSPL builds at load).
    # 3-D coords -> robust for planar AND spherical meshes.
    tree = cKDTree(coords)
    npts = coords.shape[0]
    parts = sorted(glob.glob(os.path.join(h5dir, "%s.%d.p*.h5" % (file_base, step))))
    names = fields or _list_fields(h5dir, file_base, step)
    glob_fields = {n: np.full(npts, np.nan) for n in names}
    for pth in parts:
        p = os.path.basename(pth).split(".p")[-1].split(".")[0]
        tfile = os.path.join(h5dir, "topology.p%s.h5" % p)
        with h5py.File(tfile, "r") as tf:
            lc = np.asarray(tf["coords"], dtype=np.float64)
        _, idx = tree.query(lc)
        with h5py.File(pth, "r") as f:
            for n in names:
                if n in f:
                    glob_fields[n][idx] = np.asarray(f[n])[:, 0]

    # Regular grid (degrees for a geographic mesh, mesh units otherwise).
    if spacing is None and nx is None:
        spacing = _median_spacing(np.column_stack([x, y]), cells)
    if spacing is not None:
        dx = dy = float(spacing) if np.isscalar(spacing) else float(spacing[0])
        if not np.isscalar(spacing):
            dy = float(spacing[1])
        xs = np.arange(x.min(), x.max() + dx, dx)
        ys = np.arange(y.min(), y.max() + dy, dy)
    else:
        xs = np.linspace(x.min(), x.max(), int(nx))
        ys = np.linspace(y.min(), y.max(), int(ny or nx))
        dx = xs[1] - xs[0]
        dy = ys[1] - ys[0]

    # Poles: the lon/lat grid is singular there (meridians converge). Crop the
    # latitude range to drop the polar caps — avoids the singularity and the
    # redundant pole-row cells. Default keeps everything below |89°|.
    if geographic:
        lim = 89.9 if latlim is None else float(latlim)
        ys = ys[(ys >= -lim) & (ys <= lim)]

    # Periodic longitude: a geographic grid spanning ~360°. Drop the duplicated
    # wrap meridian so column nx-1 + dx wraps onto column 0.
    periodic = bool(geographic and (x.max() - x.min() > 350.0))
    if periodic and xs[-1] - xs[0] >= 360.0 - 1.0e-9:
        xs = xs[:-1]
    gx, gy = np.meshgrid(xs, ys)

    # Build the interpolation operator(s) ONCE (point location + barycentric
    # weights), then reuse for every field — the per-field cost drops to a
    # gather + weighted sum (was a full triangulation + re-location per field).
    # Node -> incident-triangle adjacency and the per-grid-point nearest mesh
    # node: both depend only on the geometry, so they are built ONCE and shared
    # by both longitude framings below. For a geographic mesh the nearest-node
    # query reuses the KDTree already built above for the partition reassembly
    # and runs in 3-D, so it is seam-independent.
    csr = _node_tri_csr(cells, npts)
    if geographic:
        seed = _nearest_mesh_node(tree, gx, gy, R=R)
    else:
        seed = _nearest_mesh_node(cKDTree(np.column_stack([x, y])), gx, gy)

    if periodic:
        # Two longitude framings (seam at ±180° and at 0°/360°), each with its
        # seam-spanning triangles masked so their no-data stripes fall on
        # DIFFERENT meridians; merging takes frame A where valid and frame B for
        # its seam stripe, so the field is gap-free across the antimeridian.
        # Each operator is precomputed once and reused for every field; the
        # adjacency and the nearest-node seed are shared between the two.
        _interp_a = _build_tri_interp(x, y, cells, gx, gy, seed, csr,
                                      _seam_mask(x[cells]))
        lonb = x % 360.0
        _interp_b = _build_tri_interp(lonb, y, cells, gx % 360.0, gy, seed, csr,
                                      _seam_mask(lonb[cells]))

        def _raster(vals):
            za = _interp_a(vals)
            zb = _interp_b(vals)
            return np.where(np.isfinite(za), za, zb)
    else:
        # Mask any antimeridian-spanning triangle (regional geographic grids).
        tri_mask = None
        if geographic:
            tlon = x[cells]
            tri_mask = (tlon.max(axis=1) - tlon.min(axis=1)) > 180.0
        _interp = _build_tri_interp(x, y, cells, gx, gy, seed, csr, tri_mask)

        def _raster(vals):
            return _interp(vals)

    grids = {}
    for n in names:
        grids[n] = np.ma.filled(_raster(glob_fields[n]), np.nan)
    elev = grids.get("elev")
    if elev is None:
        raise ValueError("the mesh output has no 'elev' field to build hydrology")
    mask = np.isfinite(elev)

    # Hydrology runs on the SUBAERIAL cells: outlets are then the shoreline
    # (cells next to sub-base-level / marine) and the domain edge — giving
    # proper river basins draining to base level, with chi measured from it.
    # Base level defaults to the run's SEA LEVEL (catchments drain to the coast,
    # chi is measured from it), read from the step's xmf; falls back to 0.
    if base_level is None:
        base_level = _read_sealevel(h5dir, file_base, step)
        if base_level is None:
            base_level = 0.0
    hydro_mask = mask & (elev > base_level)

    # Cell metrics for the hydrology: planar -> uniform dx, dy (mesh units);
    # geographic -> dy = R.dlat, per-row dx = R.cos(lat).dlon (metres), so area
    # / chi / distance are physical despite the lon/lat raster.
    if geographic:
        hdy = R * np.radians(dy)
        hdx = np.clip(R * np.cos(np.radians(ys)) * np.radians(dx), 1.0, None)
    else:
        hdx, hdy = dx, dy
    hydro = _d8_hydrology(elev, hydro_mask, hdx, hdy, mn, a0, periodic=periodic)
    grids.update(hydro["grids"])

    out = {
        "x": xs, "y": ys, "mask": mask, "spacing": (dx, dy),
        "geographic": geographic, "periodic": periodic,
        "base_level": float(base_level),
        "receiver": hydro["receiver"], "order": hydro["order"],
        "main_basin": hydro["main_basin"],
    }
    out.update(grids)
    return out


def _median_spacing(xy, cells):
    """Median triangle-edge length — a sensible default grid spacing."""
    e = np.vstack([cells[:, [0, 1]], cells[:, [1, 2]], cells[:, [2, 0]]])
    d = np.linalg.norm(xy[e[:, 0]] - xy[e[:, 1]], axis=1)
    return float(np.median(d))


# ---------------------------------------------------------------------------
# Raster D8 hydrology
# ---------------------------------------------------------------------------

# 8-neighbour offsets (drow, dcol) and their unit-cell distances factor.
_NB = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]


def _d8_hydrology(elev, mask, dx, dy, mn, a0, periodic=False):
    """
    Priority-flood (+epsilon) fill -> D8 receivers -> drainage area, basins,
    chi and flow distance, on the cells flagged in ``mask`` (pass the SUBAERIAL
    mask, ``valid & elev > base_level``, so the shoreline + domain edge are the
    outlets and marine cells are excluded — the standard basin/chi definition).

    ``periodic`` (global lon/lat grids) wraps the **longitude** neighbours
    (column 0 <-> nx-1) so flow routes continuously across the antimeridian; the
    left/right columns are then NOT outlets — only the poles (top/bottom rows)
    and coastlines are. Returns flat ``receiver`` (-1 at outlets), the
    processing ``order`` and the field grids (NaN / -1 outside ``mask``).
    """
    ny, nx = elev.shape
    n = ny * nx
    z = np.where(mask, elev, np.inf).astype(np.float64)
    # Neighbour distances. `dx` may be a per-ROW array (geographic grids, where
    # the east-west cell size R*cos(lat)*dlon shrinks with latitude); `dy` is
    # uniform. Each cell uses its own row's distances (the usual D8 diagonal
    # approximation).
    dx_row = np.full(ny, float(dx)) if np.isscalar(dx) else np.asarray(dx, float)
    diag_row = np.hypot(dx_row, dy)
    dyc = np.full(ny, dy)
    dist_row = np.stack(                                # (ny, 8) in _NB order
        [diag_row, dyc, diag_row, dx_row, dx_row, diag_row, dyc, diag_row],
        axis=1,
    )
    cellarea_row = dx_row * dy                          # (ny,)
    eps = max(1.0e-6, (np.nanmax(elev) - np.nanmin(elev)) * 1.0e-6)

    def _neighbours(r, c):
        """Yield (k, rr, cc) valid 8-neighbours; longitude wraps if periodic."""
        for k, (dr, dc) in enumerate(_NB):
            rr = r + dr
            if rr < 0 or rr >= ny:             # poles never wrap in latitude
                continue
            cc = c + dc
            if periodic:
                cc %= nx
            elif cc < 0 or cc >= nx:
                continue
            yield k, rr, cc

    filled = np.where(mask, z, np.inf).copy()
    done = ~mask
    is_seed = np.zeros((ny, nx), dtype=bool)   # domain-boundary outlets
    heap = []
    # Seeds: valid cells on the domain boundary — next to an invalid cell, or on
    # the grid edge (poles always; left/right columns too UNLESS periodic, where
    # longitude wraps). They drain off-domain and become basin outlets.
    for r in range(ny):
        for c in range(nx):
            if not mask[r, c]:
                continue
            edge = r == 0 or r == ny - 1
            if not periodic:
                edge = edge or c == 0 or c == nx - 1
            if not edge:
                edge = any(not mask[rr, cc] for _, rr, cc in _neighbours(r, c))
            if edge:
                heapq.heappush(heap, (filled[r, c], r * nx + c))
                done[r, c] = True
                is_seed[r, c] = True
    while heap:
        e, flat = heapq.heappop(heap)
        r, c = divmod(flat, nx)
        for _, rr, cc in _neighbours(r, c):
            if not done[rr, cc]:
                filled[rr, cc] = max(z[rr, cc], e + eps)
                done[rr, cc] = True
                heapq.heappush(heap, (filled[rr, cc], rr * nx + cc))

    # Steepest-descent (D8) receiver on the filled surface.
    fflat = filled.ravel()
    recv = np.full(n, -1, dtype=np.int64)
    recv_dist = np.zeros(n)                     # distance cell -> its receiver
    valid_flat = mask.ravel()
    seed_flat = is_seed.ravel()
    vidx = np.where(valid_flat)[0]
    for flat in vidx:
        if seed_flat[flat]:
            continue                           # boundary outlet: drains off-domain
        r, c = divmod(int(flat), nx)
        best, bslope, bdist = -1, 0.0, 0.0
        fz = fflat[flat]
        for k, rr, cc in _neighbours(r, c):
            if mask[rr, cc]:
                slope = (fz - filled[rr, cc]) / dist_row[r, k]
                if slope > bslope:
                    bslope, best, bdist = slope, rr * nx + cc, dist_row[r, k]
        recv[flat] = best                      # -1 stays -> outlet
        recv_dist[flat] = bdist

    # Process in increasing filled elevation (downstream-first) for basins /
    # chi / distance, and the reverse for area accumulation.
    order = vidx[np.argsort(fflat[vidx], kind="stable")]
    area = np.where(valid_flat, cellarea_row[np.arange(n) // nx], 0.0)
    for flat in order[::-1]:                   # high -> low: add to receiver
        rc = recv[flat]
        if rc >= 0:
            area[rc] += area[flat]

    basin = np.full(n, -1, dtype=np.int64)
    chi = np.zeros(n)
    fdist = np.zeros(n)
    nextid = 0
    for flat in order:                         # low -> high: inherit from recv
        rc = recv[flat]
        if rc < 0:                             # boundary / shoreline outlet
            basin[flat] = nextid
            nextid += 1
            chi[flat] = 0.0
            fdist[flat] = 0.0
        else:
            basin[flat] = basin[rc]
            dl = recv_dist[flat]
            chi[flat] = chi[rc] + (a0 / max(area[flat], 1.0e-12)) ** mn * dl
            fdist[flat] = fdist[rc] + dl

    def _grid(flatarr, fill=np.nan, dtype=float):
        g = np.full(n, fill, dtype=dtype)
        g[valid_flat] = flatarr[valid_flat]
        return g.reshape(ny, nx)

    # Largest basin (by cell count) for convenience.
    blab = basin[valid_flat]
    main_basin = int(np.bincount(blab[blab >= 0]).argmax()) if blab.size else -1

    grids = {
        "drainage_area": _grid(area),
        "basin": _grid(basin, fill=-1, dtype=np.int64),
        "chi": _grid(chi),
        "flowdist": _grid(fdist),
        # Priority-flood-filled (hydrologically-conditioned) elevation: pits /
        # lakes raised to their spill level. A long profile drawn on THIS is
        # strictly monotonic upstream (each cell drains to a lower receiver),
        # unlike the raw `elev` which keeps real depressions + interpolation
        # roughness as small peaks/dips.
        "filled": _grid(fflat),
    }
    return {"grids": grids, "receiver": recv, "order": order,
            "main_basin": main_basin}


# ---------------------------------------------------------------------------
# NetCDF export
# ---------------------------------------------------------------------------

# (units, long_name/definition) attached to each NetCDF variable so the file is
# self-describing (CF-style). Unknown fields are written without attributes.
_VAR_META = {
    "elev": ("m", "surface elevation"),
    "erodep": ("m", "cumulative erosion (negative) / deposition (positive)"),
    "EDrate": ("m/yr", "erosion (negative) / deposition (positive) rate"),
    "FA": ("m3/yr", "flow accumulation (water discharge)"),
    "fillFA": ("m3/yr", "flow accumulation over the depression-filled surface"),
    "waterFill": ("m", "filled water-surface elevation (lake / depression level)"),
    "sedLoad": ("m3/yr", "river sediment load (total)"),
    "sedLoadF": ("m3/yr", "river sediment load (fine fraction)"),
    "drainage_area": ("m2", "drainage area from D8 routing on the gridded surface"),
    "basin": ("1", "drainage-basin id (subaerial; -1 = marine / outside)"),
    "chi": ("m", "chi, integral of (A0/A)^(m/n) dl from the outlet"),
    "flowdist": ("m", "flow distance along the network to the outlet"),
    "filled": ("m", "priority-flood-filled (hydrologically-conditioned) elevation"),
    "flexIso": ("m", "cumulative isostatic (flexural) response"),
    "rain": ("m/yr", "rainfall (precipitation) rate"),
    # Groundwater + Level-B geochemistry fields.
    "recharge": ("m/yr", "groundwater recharge"),
    "wtable": ("m", "water-table elevation"),
    "wtdepth": ("m", "water-table depth below the surface"),
    "baseflow": ("m3/yr", "groundwater baseflow to rivers"),
    "duricrust": ("m", "duricrust thickness"),
    "induration": ("1", "duricrust induration degree (0-1)"),
    "Karmor": ("1", "erodibility armoring multiplier from the duricrust"),
    "solute": ("kg/m3", "dissolved solute concentration (summed over species)"),
    "soluteflux": ("m3/yr", "groundwater dissolved-solute export / baseflow (summed over species)"),
    "riverSolute": ("m3/yr", "river dissolved load routed downstream (summed over species)"),
    "marineSoluteInput": ("m3/yr", "dissolved solute entering the ocean at coast / outlet exits"),
    "crust_type": ("1", "dominant crust-forming solute species (index; -1 = none)"),
    "crust_source": ("1", "dominant source-rock class of the crust (index; -1 = none)"),
    # spoken-name aliases, in case a field is written under these names
    "flexiso": ("m", "cumulative isostatic (flexural) response"),
    "rainfall": ("m/yr", "rainfall (precipitation) rate"),
    "precipitation": ("m/yr", "rainfall (precipitation) rate"),
}


def to_netcdf(result, path, time=None):
    """
    Write the gridded surface to a CF-1.x NetCDF (1-D coordinate variables +
    2-D fields), readable by PyGMT (``grdimage`` etc.) and ArcGIS. Geographic
    grids use CF ``lon``/``lat`` (``degrees_east``/``degrees_north``) so they
    are recognised as geographic; planar grids use ``x``/``y``.
    """
    import netCDF4

    skip = {"x", "y", "mask", "spacing", "receiver", "order", "main_basin",
            "geographic", "periodic", "base_level"}
    geo = result.get("geographic", False)
    xname, yname = ("lon", "lat") if geo else ("x", "y")
    ny, nx = result["y"].size, result["x"].size
    with netCDF4.Dataset(path, "w", format="NETCDF4") as ds:
        ds.Conventions = "CF-1.7"
        ds.createDimension(xname, nx)
        ds.createDimension(yname, ny)
        xv = ds.createVariable(xname, "f8", (xname,))
        yv = ds.createVariable(yname, "f8", (yname,))
        xv[:] = result["x"]
        yv[:] = result["y"]
        if geo:
            xv.units, xv.standard_name = "degrees_east", "longitude"
            yv.units, yv.standard_name = "degrees_north", "latitude"
        else:
            xv.long_name, yv.long_name = "x", "y"
            xv.units = yv.units = "m"
        if time is not None:
            ds.time = float(time)
            ds.time_units = "yr"
        # Sea level used as the hydrology base level (drainage outlets + chi
        # datum). Recorded as a global attribute and a scalar variable so the
        # NetCDF is self-describing.
        if result.get("base_level") is not None:
            sl = float(result["base_level"])
            ds.sea_level = sl
            sv = ds.createVariable("sea_level", "f8", ())
            sv[...] = sl
            sv.long_name = "sea level used as hydrology base level"
            sv.units = "m"
        for name, grid in result.items():
            if name in skip or not isinstance(grid, np.ndarray) or grid.ndim != 2:
                continue
            dt = "i4" if grid.dtype.kind in "iu" else "f8"
            fill = -1 if dt == "i4" else np.nan
            v = ds.createVariable(name, dt, (yname, xname), fill_value=fill,
                                  zlib=True)
            v[:, :] = grid
            meta = _VAR_META.get(name)
            if meta is None:
                # Per-species geochem fields (solute_<name>, soluteflux_<name>,
                # riverSolute_<name>, crust_<name>) — describe by their prefix.
                for pre, unit, base in (
                    ("riverSolute_", "m3/yr", "river dissolved load"),
                    ("soluteflux_", "m3/yr", "groundwater dissolved-solute export"),
                    ("solute_", "kg/m3", "dissolved solute concentration"),
                    ("crust_", "m", "duricrust contribution"),
                ):
                    if name.startswith(pre):
                        meta = (unit, "%s — species '%s'" % (base, name[len(pre):]))
                        break
            if meta is not None:
                v.units, v.long_name = meta
    return path


# ---------------------------------------------------------------------------
# River profiles for a basin
# ---------------------------------------------------------------------------

def basin_rivers(result, basin_id=None, area_threshold=None):
    """
    Extract the channel network of one basin: cells whose drainage area exceeds
    ``area_threshold`` and that belong to ``basin_id`` (default: the largest
    basin). Returns a dict with the **main stem** (traced from the outlet up the
    largest-area donor at each step) and the list of **tributaries** (each
    channel head traced down to the main stem), every path carrying
    ``x``/``y``/``dist``/``elev``/``chi``/``area`` arrays ordered outlet->source.
    """
    ny, nx = result["y"].size, result["x"].size
    basin = result["basin"]
    area = result["drainage_area"]
    elev = result["elev"]
    chi = result["chi"]
    fdist = result["flowdist"]
    recv = result["receiver"]
    xs, ys = result["x"], result["y"]
    if basin_id is None:
        basin_id = result["main_basin"]
    if area_threshold is None:
        area_threshold = float(np.nanpercentile(area[basin == basin_id], 95))

    inb = (basin == basin_id)
    chan = inb & (area >= area_threshold)
    chan_flat = chan.ravel()

    # Donors per cell (who flows into it) to find main stem + channel heads.
    donors = [[] for _ in range(ny * nx)]
    cf = chan.ravel()
    for flat in np.where(cf)[0]:
        rc = recv[flat]
        if rc >= 0 and cf[rc]:
            donors[rc].append(int(flat))

    # Outlet of this basin = the channel cell with the smallest flow distance.
    chan_idx = np.where(cf)[0]
    if chan_idx.size == 0:
        return {"basin_id": basin_id, "main_stem": None, "tributaries": [],
                "area_threshold": area_threshold}
    fd = fdist.ravel()
    outlet = int(chan_idx[np.argmin(fd[chan_idx])])

    def _path_up_mainstem(start):
        path = [start]
        while donors[path[-1]]:
            nxt = max(donors[path[-1]], key=lambda d: area.ravel()[d])
            path.append(nxt)
        return path

    main = _path_up_mainstem(outlet)
    on_main = np.zeros(ny * nx, dtype=bool)
    on_main[main] = True

    filled = result.get("filled")
    fillflat = filled.ravel() if filled is not None else None

    def _pack(flat_list):
        a = np.array(flat_list)
        r, c = np.divmod(a, nx)
        out = {
            "x": xs[c], "y": ys[r],
            "dist": fd[a],
            "elev": elev.ravel()[a],
            "chi": chi.ravel()[a],
            "area": area.ravel()[a],
        }
        # Hydrologically-filled elevation (monotonic upstream) when available.
        if fillflat is not None:
            out["elev_filled"] = fillflat[a]
        return out

    # Tributaries: each channel head (no channel donor) not on the main stem,
    # traced DOWN its receivers until it meets the main stem.
    heads = [int(f) for f in chan_idx
             if not donors[f] and not on_main[f]]
    tribs = []
    for h in heads:
        seg = [h]
        cur = recv[h]
        while cur >= 0 and cf[cur] and not on_main[cur]:
            seg.append(int(cur))
            cur = recv[cur]
        if cur >= 0 and on_main[cur]:
            seg.append(int(cur))          # confluence point on the main stem
        tribs.append(_pack(seg[::-1]))    # outlet-ward order

    return {
        "basin_id": basin_id,
        "area_threshold": area_threshold,
        "main_stem": _pack(main[::-1]),   # outlet -> source... reverse to source->outlet
        "tributaries": tribs,
    }


def plot_long_profile(rivers, xaxis="dist", ax=None, which="elev", figsize=(7, 4)):
    """
    Longitudinal profile (``xaxis`` = ``'dist'`` distance-to-outlet, or
    ``'chi'``) vs elevation: main stem bold, tributaries thin.

    :arg which: ``'elev'`` (default) plots the **raw** gridded elevation — which
        keeps real depressions/lakes and TIN-interpolation roughness as small
        peaks/dips; ``'filled'`` plots the priority-flood-filled elevation, which
        is **strictly monotonic** upstream (each cell drains to a lower
        receiver) — the smooth, hydrologically-conditioned profile.
    """
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    key = "elev_filled" if which == "filled" else "elev"
    ms = rivers["main_stem"]
    for t in rivers["tributaries"]:
        ax.plot(t[xaxis], t.get(key, t["elev"]), color="0.6", lw=0.6)
    if ms is not None:
        ax.plot(ms[xaxis], ms.get(key, ms["elev"]), color="C3", lw=2.0,
                label="main stem")
    if xaxis == "chi":
        ax.set_xlabel("χ")
    else:
        from matplotlib.ticker import FuncFormatter
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: "%g" % (v / 1000.0)))
        ax.set_xlabel("distance to outlet (km)")
    ax.set_ylabel("elevation (m)" + (" (filled)" if which == "filled" else ""))
    ax.set_title("Basin %s longitudinal profile" % rivers["basin_id"])
    ax.legend(loc="best")
    return ax


def plot_basin_map(result, rivers, background="elev", ax=None, figsize=(7, 5),
                   sea_level=None):
    """
    Map the channel network on the surface: ``background`` field as an image,
    main stem (red) + tributaries (cyan) overlaid, plus the **sea-level
    coastline** (the ``elev == sea_level`` contour). Returns the axes.

    :arg sea_level: elevation of the coastline contour; defaults to the run's
        sea level used for the hydrology (``result['base_level']``). Pass a float
        to override, or ``None`` is used (no line) if neither is available.
    """
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(figsize=figsize)
    xs, ys = result["x"], result["y"]
    bg = result[background]
    ax.imshow(bg, origin="lower", extent=[xs[0], xs[-1], ys[0], ys[-1]],
              cmap="terrain", aspect="equal")
    if sea_level is None:
        sea_level = result.get("base_level")
    if sea_level is not None:
        X, Y = np.meshgrid(xs, ys)
        ax.contour(X, Y, result["elev"], levels=[float(sea_level)],
                   colors="steelblue", linewidths=1.4, linestyles="--")
        # Legend proxy (contour sets aren't directly legendable across mpl
        # versions; QuadContourSet.collections was removed in mpl 3.8+).
        ax.plot([], [], color="steelblue", lw=1.4, ls="--",
                label="sea level (%.0f m)" % float(sea_level))
    for t in rivers["tributaries"]:
        ax.plot(t["x"], t["y"], color="c", lw=0.7)
    ms = rivers["main_stem"]
    if ms is not None:
        ax.plot(ms["x"], ms["y"], color="r", lw=1.8, label="main stem")
    ax.legend(loc="best", fontsize=8)
    ax.set_title("Basin %s channel network" % rivers["basin_id"])
    # Planar (UTM, m) grids: label axes in km; geographic grids stay in degrees.
    if not result.get("geographic", False):
        from matplotlib.ticker import FuncFormatter
        kmf = FuncFormatter(lambda v, _: "%g" % (v / 1000.0))
        ax.xaxis.set_major_formatter(kmf)
        ax.yaxis.set_major_formatter(kmf)
        ax.set_xlabel("x (km)")
        ax.set_ylabel("y (km)")
    else:
        ax.set_xlabel("longitude")
        ax.set_ylabel("latitude")
    return ax


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    p = argparse.ArgumentParser(
        description="Rasterise a goSPL surface to a CF-NetCDF grid (fields + "
        "basins + chi) for PyGMT / ArcGIS."
    )
    p.add_argument("--h5dir", required=True, help="goSPL output 'h5' directory")
    p.add_argument("--mesh", required=True,
                   help="global mesh as FILE.npz[:vkey[:ckey]] (default keys v,c)")
    p.add_argument("--out", default="surface.nc", help="output NetCDF (default surface.nc)")
    p.add_argument("--step", type=int, default=None, help="output step (default: last)")
    p.add_argument("--spacing", default=None,
                   help="grid spacing DX[,DY] in mesh units (default: median edge)")
    p.add_argument("--fields", default=None, help="comma-separated fields (default: all)")
    p.add_argument("--mn", type=float, default=0.5, help="chi m/n concavity (default 0.5)")
    p.add_argument("--a0", type=float, default=1.0, help="chi reference area (default 1)")
    p.add_argument("--base-level", type=float, default=None,
                   help="elevation defining the coast/outlet; catchments + chi "
                        "are measured from it (default: the run's sea level, "
                        "read from the step xmf; else 0)")
    p.add_argument("--file-base", default="gospl", help="mesh file base (default gospl)")
    p.add_argument("--latlim", type=float, default=None,
                   help="geographic only: crop |latitude| to this limit, dropping "
                        "the singular polar caps (default 89.9)")
    p.add_argument("--tout", type=float, default=None,
                   help="output interval (yr) -> NetCDF time = step*tout")
    p.add_argument("--tstart", type=float, default=0.0)
    args = p.parse_args(argv)

    parts = args.mesh.split(":")
    mesh = parts[0]
    vkey = parts[1] if len(parts) > 1 else "v"
    ckey = parts[2] if len(parts) > 2 else "c"
    spacing = None
    if args.spacing is not None:
        s = [float(v) for v in args.spacing.split(",")]
        spacing = s[0] if len(s) == 1 else s
    fields = args.fields.split(",") if args.fields else None

    g = grid_export(args.h5dir, mesh, step=args.step, vkey=vkey, ckey=ckey,
                    spacing=spacing, fields=fields, mn=args.mn, a0=args.a0,
                    base_level=args.base_level, file_base=args.file_base,
                    latlim=args.latlim)
    time = None
    if args.tout is not None and args.step is not None:
        time = args.tstart + args.step * args.tout
    to_netcdf(g, args.out, time=time)
    print("wrote %s (%d x %d grid; %d basins; sea level %.3f m) "
          "— open in PyGMT / ArcGIS"
          % (args.out, g["x"].size, g["y"].size,
             int(g["basin"].max()) + 1, g["base_level"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
