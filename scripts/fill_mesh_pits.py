#!/usr/bin/env python3
"""Depression-fill (priority-flood + epsilon) an input goSPL mesh.

Standalone pre-processing utility: it reads the `npdata` npz an input YAML
points at (vertices / cells / elevation), removes every closed depression above
a chosen base level, and writes a new npz with a strictly drainable elevation
field.

This is NOT part of `Model.runProcesses`. goSPL already fills depressions
internally every step (`flow/pitfilling.py`) and keeps the *unfilled* surface as
the real elevation, so pre-filling the input permanently deletes lakes and
endorheic basins from the initial condition. Do it only when that is what you
want (e.g. a drainage-consistent starting DEM).

Two backends, same seeding contract:

* ``fortran`` (default) calls the model's own `globalngbhs` + `epsfill`
  kernels, so the result is bit-identical to what `seaplex._matOcean` computes
  internally. Its increment is one ULP (`nearest`), i.e. ~1e-13 m at 1000 m
  elevation: enough to order a double-precision descent in memory, but it does
  NOT survive a float32 round-trip and is invisible to any downstream tool that
  rounds. Fine when the filled surface is consumed immediately in float64.
* ``python`` takes ``--epsilon`` (metres) and applies the same Barnes (2014)
  priority-flood with a physical increment, which is what you want when the
  filled surface is written to disk and re-read later.

Seeding: a priority-flood needs outlets. Pick one of

* ``--sea-level Z``  every node strictly below Z is an outlet (the marine
  domain; this is exactly what the model does, with Z = sealevel + oFill).
* ``--borders``      every node on the mesh hull is an outlet (planar meshes).
* ``--outlets IDS``  explicit node ids, from a .npy file or a comma list.

Outlets are never raised, whichever mode is used, so a closed basin that dips
below ``--sea-level`` stays a sink by construction.

``--max-depth`` / ``--max-volume`` / ``--max-cells`` keep the big basins: a
depression exceeding any limit given is left at its original elevation, so a
run removes the interpolation noise and the NODATA flats without deleting the
real endorheic basins. Every depression is measured and the largest are
printed even with no limit set, so one run tells you where to put one.

``--outlet-mode`` then decides what an outlet means to its upslope neighbours:

* ``elevation`` (default) seeds the flood at each outlet's own elevation, the
  classic DEM convention: an interior cell has to be raised over the rim to
  reach an outlet, and the filled surface descends all the way to the outlet.
* ``absorbing`` treats outlets as bottomless, which is what the model does to a
  planar mesh in `seaplex._matOcean` (it forces `outletIDs` to a -1e6
  sentinel) because flow leaves an open/fixed edge regardless of its
  elevation. Less filling, but the filled surface does not descend into the
  outlet node itself.

Examples
--------
    python scripts/fill_mesh_pits.py mesh.npz -o mesh_filled.npz --sea-level 0.
    python scripts/fill_mesh_pits.py flat.npz -o flat_filled.npz --borders \
        --backend python --epsilon 1.e-4

User documentation: `docs/user_guide/inputfile.rst`, section "Conditioning the
initial topography".
"""

import argparse
import heapq
import sys

import numpy as np


# `globalngbhs` fills the module-level FVgnID(nt, 12) neighbour table that
# `epsfill` walks; a vertex with more neighbours overflows it (no bound check in
# the Fortran), so the degree is verified before the call.
MAX_FORTRAN_NGBH = 12


def read_mesh(path, coords_key, cells_key, elev_key):
    """Return (all arrays, coords, cells, elevation) from a goSPL input npz."""
    data = np.load(path)
    for key in (coords_key, cells_key, elev_key):
        if key not in data.files:
            raise KeyError(
                "'%s' is not in %s (found: %s)" % (key, path, ", ".join(data.files))
            )
    arrays = {k: data[k] for k in data.files}
    coords = np.asarray(arrays[coords_key], dtype=np.float64)
    cells = np.asarray(arrays[cells_key], dtype=np.int32)
    elev = np.asarray(arrays[elev_key], dtype=np.float64).ravel()
    if cells.ndim != 2 or cells.shape[1] != 3:
        raise ValueError("'%s' must be an (ncells, 3) triangle array" % cells_key)
    if len(elev) != len(coords):
        raise ValueError(
            "elevation ('%s', %d) and coordinates ('%s', %d) disagree"
            % (elev_key, len(elev), coords_key, len(coords))
        )
    return arrays, coords, cells, elev


def edge_list(cells):
    """Directed edge array (3*ncells, 2) of every triangle side."""
    return np.vstack((cells[:, [0, 1]], cells[:, [1, 2]], cells[:, [2, 0]]))


def neighbour_csr(npoints, cells):
    """Symmetric vertex adjacency as a CSR (indptr, indices) pair."""
    edges = edge_list(cells)
    both = np.vstack((edges, edges[:, ::-1])).astype(np.int64)
    both = np.unique(both, axis=0)
    counts = np.bincount(both[:, 0], minlength=npoints)
    indptr = np.zeros(npoints + 1, dtype=np.int64)
    np.cumsum(counts, out=indptr[1:])
    return indptr, np.ascontiguousarray(both[:, 1])


def hull_nodes(cells):
    """Nodes on a boundary edge (an edge used by a single triangle)."""
    key = np.sort(edge_list(cells).astype(np.int64), axis=1)
    uniq, counts = np.unique(key, axis=0, return_counts=True)
    return np.unique(uniq[counts == 1])


def vertex_areas(coords, cells):
    """Per-vertex cell area: a third of each incident triangle.

    This is the barycentric dual, not the Voronoi area goSPL itself builds in
    `definetin`, which needs the circumcentres. On the near-uniform Delaunay
    meshes goSPL runs the two agree to a few percent, which is far inside the
    precision a volume threshold needs.
    """
    p = coords[cells[:, 0]], coords[cells[:, 1]], coords[cells[:, 2]]
    tri = 0.5 * np.linalg.norm(np.cross(p[1] - p[0], p[2] - p[0]), axis=1) / 3.0
    areas = np.zeros(len(coords))
    for k in range(3):
        areas += np.bincount(cells[:, k], weights=tri, minlength=len(coords))
    return areas


def depression_stats(filled, elev, indptr, indices, areas):
    """Label each filled depression and measure it.

    A depression is a connected component of the raised set. Returns
    ``(labels, depth, volume, ncells)``: `labels` is -1 on unraised nodes and a
    component index elsewhere; the three arrays are per component.

    Caveat worth knowing before thresholding on the result: two depressions
    that end up contiguous after filling (a small pit spilling into a large
    basin, say) merge into ONE component and are measured together, so the
    pair is judged by its combined size. That errs toward keeping, which is
    the safe direction for a "preserve the real basins" switch.
    """
    try:
        from scipy.sparse import csr_matrix
        from scipy.sparse.csgraph import connected_components
    except ImportError:  # pragma: no cover - scipy is a goSPL dependency
        raise SystemExit(
            "the depression-size limits need scipy (a goSPL dependency): "
            "pip install scipy"
        )

    raised = filled > elev
    npoints = len(elev)
    labels = np.full(npoints, -1, dtype=np.int64)
    if not raised.any():
        return labels, np.empty(0), np.empty(0), np.empty(0, dtype=np.int64)

    # Keep only edges with BOTH ends raised, so every unraised node is an
    # isolated singleton and cannot bridge two distinct depressions.
    rows = np.repeat(raised, np.diff(indptr))
    data = (raised[indices] & rows).astype(np.int8)
    graph = csr_matrix((data, indices, indptr), shape=(npoints, npoints))
    graph.eliminate_zeros()
    _, comp = connected_components(graph, directed=False)

    # Re-index so only components that actually contain raised nodes are kept.
    lab = comp[raised]
    order = np.argsort(lab, kind="stable")
    sorted_lab = lab[order]
    starts = np.flatnonzero(np.r_[True, sorted_lab[1:] != sorted_lab[:-1]])
    dz = (filled - elev)[raised][order]
    dv = ((filled - elev) * areas)[raised][order]

    depth = np.maximum.reduceat(dz, starts)
    volume = np.add.reduceat(dv, starts)
    ncells = np.diff(np.r_[starts, len(sorted_lab)])

    remap = np.full(comp.max() + 1, -1, dtype=np.int64)
    remap[sorted_lab[starts]] = np.arange(len(starts))
    labels[raised] = remap[lab]
    return labels, depth, volume, ncells


def fill_fortran(elev, cells, outlets, base_level, absorbing):
    """Fill via the model's own `epsfill` (increment = 1 ULP).

    `epsfill` seeds on `elev < base_level` alone, so an arbitrary outlet mask is
    expressed the way `seaplex._matOcean` does it for a planar mesh: push the
    outlets to a sentinel below every real elevation, fill, restore. That is the
    absorbing reading by construction, which is why the elevation reading is
    only available here when the mask already IS the sub-`base_level` set.
    """
    from gospl._fortran import epsfill, globalngbhs

    npoints = len(elev)
    native = np.array_equal(outlets, elev < base_level)
    if native and not absorbing:
        work, cut_off = elev, base_level
    else:
        sentinel = float(elev.min()) - 1000.0
        work = elev.copy()
        work[outlets] = sentinel
        cut_off = sentinel + 500.0

    nover = globalngbhs(npoints, cells)
    if nover:
        raise SystemExit(
            "%d vertices exceed the %d-neighbour Fortran table and lost "
            "connections; use --backend python" % (nover, MAX_FORTRAN_NGBH)
        )
    filled = epsfill(cut_off, work)
    filled[outlets] = elev[outlets]
    return filled


def fill_python(elev, indptr, indices, outlets, epsilon, absorbing):
    """Barnes (2014) priority-flood + a physical epsilon, pure Python heap."""
    filled = elev.copy()
    done = np.zeros(len(elev), dtype=bool)
    done[outlets] = True

    seed_height = -np.inf if absorbing else None
    queue = []
    for node in np.flatnonzero(outlets):
        node = int(node)
        queue.append((filled[node] if seed_height is None else seed_height, node))
    heapq.heapify(queue)

    while queue:
        height, node = heapq.heappop(queue)
        floor = height + epsilon
        for nbr in indices[indptr[node]:indptr[node + 1]]:
            nbr = int(nbr)
            if done[nbr]:
                continue
            done[nbr] = True
            if filled[nbr] <= floor:
                filled[nbr] = floor
            heapq.heappush(queue, (filled[nbr], nbr))

    unreached = ~done
    if unreached.any():
        print(
            "  warning: %d node(s) never reached by the flood (disconnected "
            "from every outlet); left unchanged" % unreached.sum()
        )
    return filled


def count_sinks(elev, indptr, indices, outlets, absorbing):
    """Nodes with nowhere to send water, excluding the outlets themselves.

    Under `absorbing` an outlet swallows whatever reaches it, so touching one is
    enough to drain; otherwise a node needs a strictly lower neighbour.
    """
    # Vectorised over the CSR: a per-node reduction with `reduceat`, so this
    # stays usable on a multi-million-node mesh (a Python loop over the nodes
    # took minutes there, which defeats the point of a quick diagnostic).
    has_nbrs = np.diff(indptr) > 0
    starts = indptr[:-1][has_nbrs]
    lowest = np.full(len(elev), np.inf)
    lowest[has_nbrs] = np.minimum.reduceat(elev[indices], starts)
    drains = lowest < elev
    if absorbing:
        touches = np.zeros(len(elev), dtype=bool)
        touches[has_nbrs] = np.maximum.reduceat(
            outlets[indices].astype(np.int8), starts
        ).astype(bool)
        drains |= touches
    return int(np.count_nonzero(~drains & ~outlets & has_nbrs))


def parse_outlets(spec, npoints):
    if spec.endswith(".npy"):
        ids = np.load(spec).astype(np.int64).ravel()
    else:
        ids = np.array([int(v) for v in spec.split(",") if v.strip()], dtype=np.int64)
    if ids.size and (ids.min() < 0 or ids.max() >= npoints):
        raise ValueError("outlet ids out of range [0, %d)" % npoints)
    mask = np.zeros(npoints, dtype=bool)
    mask[ids] = True
    return mask


def normalise_numeric_argv(argv):
    """Make ``--sea-level -50.`` parse.

    argparse only accepts a leading-``-`` token as a value when it matches its
    negative-number pattern (``-50`` or ``-50.0``); the trailing-dot spelling
    goSPL uses everywhere else (``-50.``) misses it and is read as an unknown
    option, which fails with a confusing "expected one argument". Rewrite that
    one shape into an equivalent argparse does accept. The ``--sea-level=-50.``
    form is a single token and was never affected.
    """
    numeric = ("--sea-level", "--epsilon")
    out = []
    expects_value = False
    for token in argv:
        if expects_value and token.startswith("-"):
            try:
                token = repr(float(token))
            except ValueError:
                pass
        expects_value = token in numeric
        out.append(token)
    return out


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("npz", help="input goSPL mesh (the `npdata` file)")
    parser.add_argument("-o", "--out", required=True, help="output npz")
    parser.add_argument("--coords-key", default="v")
    parser.add_argument("--cells-key", default="c")
    parser.add_argument("--elev-key", default="z")
    parser.add_argument(
        "--out-key",
        default=None,
        help="store the filled field under this name instead of overwriting "
        "the elevation key (keeps the raw elevation in the file)",
    )
    seed = parser.add_mutually_exclusive_group(required=True)
    seed.add_argument(
        "--sea-level",
        type=float,
        help="nodes strictly below this elevation are outlets (and are never "
        "raised) — the model's own seeding",
    )
    seed.add_argument(
        "--borders", action="store_true", help="mesh hull nodes are the outlets"
    )
    seed.add_argument(
        "--outlets", help="explicit outlet node ids (.npy file or comma list)"
    )
    parser.add_argument(
        "--outlet-mode",
        choices=("elevation", "absorbing"),
        default="elevation",
        help="'elevation': outlets seed the flood at their own height (classic "
        "DEM fill). 'absorbing': outlets are bottomless, the model's planar "
        "convention. Default: %(default)s",
    )
    parser.add_argument("--backend", choices=("fortran", "python"), default="fortran")
    parser.add_argument(
        "--epsilon",
        type=float,
        default=1.0e-4,
        help="python backend: elevation increment per cell, in metres "
        "(default: %(default)s)",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="count the interior sinks (cells with no strictly lower "
        "neighbour) before and after filling — the quick way to tell whether a "
        "goSPL un-drained-region abort comes from the input topography",
    )
    limits = parser.add_argument_group(
        "depression size limits",
        "Keep the big basins. A depression exceeding ANY limit given is left "
        "at its original elevation, so only the small ones (interpolation "
        "noise, NODATA flats) are filled. Every depression is measured and "
        "the largest are listed, so a first run with no limit tells you where "
        "to put one.",
    )
    limits.add_argument(
        "--max-depth",
        type=float,
        default=None,
        help="fill only depressions whose deepest fill is at most this, in metres",
    )
    limits.add_argument(
        "--max-volume",
        type=float,
        default=None,
        help="fill only depressions holding at most this fill volume, in m^3",
    )
    limits.add_argument(
        "--max-cells",
        type=int,
        default=None,
        help="fill only depressions spanning at most this many cells",
    )
    if argv is None:
        argv = sys.argv[1:]
    args = parser.parse_args(normalise_numeric_argv(argv))

    if args.epsilon < 0.0:
        raise SystemExit("--epsilon must be >= 0 (got %g)" % args.epsilon)

    arrays, _coords, cells, elev = read_mesh(
        args.npz, args.coords_key, args.cells_key, args.elev_key
    )
    npoints = len(elev)
    print("mesh %s: %d nodes, %d cells, elevation %.3f .. %.3f m"
          % (args.npz, npoints, len(cells), elev.min(), elev.max()))

    if args.sea_level is not None:
        outlets = elev < args.sea_level
        base_level = args.sea_level
        seeding = "below %.3f m" % args.sea_level
    elif args.borders:
        outlets = np.zeros(npoints, dtype=bool)
        outlets[hull_nodes(cells)] = True
        base_level = float(elev.min()) - 1.0
        seeding = "mesh hull"
    else:
        outlets = parse_outlets(args.outlets, npoints)
        base_level = float(elev.min()) - 1.0
        seeding = "explicit ids"
    if not outlets.any():
        raise SystemExit(
            "no outlet node selected (%s) — the flood would have no seed and "
            "the mesh would come back unchanged" % seeding
        )
    if outlets.all():
        raise SystemExit(
            "every node is an outlet (%s) — outlets are never raised, so there "
            "is nothing left to fill" % seeding
        )
    absorbing = args.outlet_mode == "absorbing"
    print("outlets: %d node(s) (%s, %s)" % (outlets.sum(), seeding, args.outlet_mode))

    indptr, indices = neighbour_csr(npoints, cells)
    degree = np.diff(indptr)
    print("vertex degree: max %d, mean %.2f" % (degree.max(), degree.mean()))

    if args.check:
        print("interior sinks before: %d"
              % count_sinks(elev, indptr, indices, outlets, absorbing))

    if args.backend == "fortran":
        if not absorbing and args.sea_level is None:
            raise SystemExit(
                "the fortran backend seeds on 'elevation < cut-off' alone, so "
                "an explicit outlet set can only be expressed as --outlet-mode "
                "absorbing. Use --backend python for elevation-seeded outlets."
            )
        if degree.max() > MAX_FORTRAN_NGBH:
            raise SystemExit(
                "%d node(s) have more than %d neighbours (max %d). The Fortran "
                "neighbour table is fixed at %d slots and `globalngbhs` does "
                "not bound-check it, so this mesh would corrupt memory. Use "
                "--backend python."
                % ((degree > MAX_FORTRAN_NGBH).sum(), MAX_FORTRAN_NGBH,
                   degree.max(), MAX_FORTRAN_NGBH)
            )
        filled = fill_fortran(elev, cells, outlets, base_level, absorbing)
        print("backend: fortran epsfill (increment = 1 ULP)")
    else:
        filled = fill_python(elev, indptr, indices, outlets, args.epsilon, absorbing)
        print("backend: python priority-flood (increment = %g m)" % args.epsilon)

    if (filled < elev).any():
        raise SystemExit("internal error: the fill lowered %d node(s)"
                         % (filled < elev).sum())

    labels, depth, volume, ncells = depression_stats(
        filled, elev, indptr, indices, vertex_areas(_coords, cells)
    )
    if len(depth):
        print("depressions filled: %d" % len(depth))
        rank = np.argsort(volume)[::-1][:3]
        for pos, comp in enumerate(rank):
            print("  largest #%d: depth %.3f m, volume %.4g m3, %d cells"
                  % (pos + 1, depth[comp], volume[comp], ncells[comp]))

    keep = np.zeros(len(depth), dtype=bool)
    if args.max_depth is not None:
        keep |= depth > args.max_depth
    if args.max_volume is not None:
        keep |= volume > args.max_volume
    if args.max_cells is not None:
        keep |= ncells > args.max_cells
    if keep.any():
        # Restore the original surface over every depression that exceeded a
        # limit. Reverting is exact and local: a kept basin's cells go back to
        # their input elevation, and the fill level of everything else was
        # fixed by its own spill path, so nothing upstream needs redoing.
        revert = np.isin(labels, np.flatnonzero(keep))
        filled[revert] = elev[revert]
        print("kept %d depression(s) unfilled (%d cells, %.4g m3 of fill not "
              "applied): over the limit"
              % (int(keep.sum()), int(revert.sum()), float(volume[keep].sum())))

    raised = filled > elev
    print("raised %d node(s) (%.2f%%), max fill %.3f m"
          % (raised.sum(), 100.0 * raised.sum() / npoints,
             (filled - elev).max() if raised.any() else 0.0))

    if args.check:
        print("interior sinks after: %d"
              % count_sinks(filled, indptr, indices, outlets, absorbing))

    out_key = args.out_key or args.elev_key
    arrays[out_key] = filled
    np.savez_compressed(args.out, **arrays)
    print("wrote %s (elevation in '%s', float64)" % (args.out, out_key))
    return 0


if __name__ == "__main__":
    sys.exit(main())
