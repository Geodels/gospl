"""
Post-processing: build a 3-D **wedge volume** of the goSPL stratigraphy for one
or more output steps, ready to open in ParaView.

goSPL records the stratigraphic pile in the distributed ``stratal.<step>.p<rank>
.h5`` files (per-layer ``stratZ``/``stratH``/``phiS`` and, optionally,
``stratHf``/``phiF`` for dual lithology and ``stratP`` for provenance), with the
triangular surface mesh in ``topology.p<rank>.h5``. This tool stacks that
triangular mesh across the recorded layers into a **triangular-prism (wedge)
volume**: each vertex of the surface mesh is replicated at every layer's
deposition elevation ``stratZ``, and each surface triangle becomes a column of
wedge cells between consecutive layers. The layer properties are written as
**cell data** on the wedges, so colouring a slab shows that layer's composition.

Two field modes (``--field``):

* ``lithology`` — coarse / fine thickness, fine fraction and the per-fraction
  porosities (needs the dual-lithology ``stratHf``/``phiF`` fields).
* ``provenance`` — per-source-class volume fraction and the dominant source
  (needs the provenance ``stratP`` field).

The basal **bedrock sentinel** layer (the ~1e6 m infinite reservoir goSPL keeps
under the recorded deposits) is auto-detected and skipped, so the volume spans
only the real recorded layers (override with ``--first-layer``).

Output mirrors goSPL's own layout — one HDF5 per input partition per step
(``<outdir>/<out>.<step>.p<p>.h5``) plus a single self-contained
``<outdir>/<out>.xdmf`` temporal collection referencing them. It runs serially
(looping over all partitions) or under ``mpirun`` (the partitions are split
across ranks — independent of the simulation's processor count); the output is
identical either way, so the result is trivially parallel and time-evolving.

Options
-------
``--h5dir DIR``
    goSPL output ``h5`` directory (**required**).
``--outdir DIR``
    output directory, created if needed (default ``strata``).
``--out PREFIX``
    output file prefix (default ``strata``); writes ``<outdir>/<out>.xdmf`` plus
    ``<outdir>/<out>.<step>.p<p>.h5``.
``--field {lithology,provenance}``
    cell-field set (default ``lithology``). ``lithology`` = coarse/fine
    thickness, fine fraction and per-fraction porosity (needs the dual-lithology
    ``stratHf``/``phiF`` fields); ``provenance`` = per-source-class volume
    fraction + dominant source (needs ``stratP``).
``--steps S0,S1,...``
    comma-separated output steps to convert (default: all found). Steps with
    fewer than 2 recorded layers are skipped.
``--tout YR`` (with optional ``--tstart YR``)
    set the ParaView ``Time`` to the simulation year ``tstart + step*tout``
    instead of the bare step index (default: step index).
``--first-layer N``
    first layer index to include (default: auto-skip the basal bedrock-sentinel
    layer, detected as thickness >= 1e5 m).
``--mesh-base BASE``
    mesh-output file base providing the current surface ``elev`` (default
    ``gospl``). Layer interfaces are placed at the current surface minus the
    cumulative compacted thickness, so eroded layers collapse instead of
    overhanging; falls back to the recorded ``stratZ`` if the mesh file is
    absent.
``--file-base BASE``
    stratal-file base name (default ``stratal``).

Runnable as the installed command ``gospl-strata-volume`` or, with no install
step, as ``python -m gospl.analyse.stratamesh``. Examples::

    # defaults -> writes ./strata/strata.xdmf (+ per-partition .h5)
    gospl-strata-volume --h5dir myrun/h5 --field lithology --steps 0,5,10

    # provenance source mix, ParaView time in years, run in parallel
    mpirun -np 4 gospl-strata-volume --h5dir myrun/h5 --outdir vol \\
        --field provenance --tout 100000
"""

import os
import glob
import argparse

import numpy as np

try:
    import h5py
except ImportError:  # pragma: no cover - h5py is required at runtime
    h5py = None

try:
    from mpi4py import MPI

    _COMM = MPI.COMM_WORLD
    _RANK = _COMM.Get_rank()
    _SIZE = _COMM.Get_size()
except ImportError:  # pragma: no cover - serial fallback
    _COMM = None
    _RANK = 0
    _SIZE = 1

# Bedrock-sentinel thickness goSPL uses for the infinite basal reservoir; any
# layer thicker than this cut-off is treated as basement and excluded.
_SENTINEL_CUTOFF = 1.0e5


def discover_steps(h5dir, file_base="stratal"):
    """Return ``(sorted step list, partition count)`` found in ``h5dir``."""
    steps = set()
    for f in glob.glob(os.path.join(h5dir, "%s.*.p*.h5" % file_base)):
        name = os.path.basename(f)
        try:
            steps.add(int(name.split(".")[1]))
        except (IndexError, ValueError):
            continue
    if not steps:
        raise FileNotFoundError(
            "no %s.<step>.p*.h5 files in %s" % (file_base, h5dir)
        )
    nparts = len(
        glob.glob(os.path.join(h5dir, "%s.%d.p*.h5" % (file_base, min(steps))))
    )
    return sorted(steps), nparts


def detect_first_layer(stratal_path):
    """
    Index of the first non-basement layer: skip leading layers whose thickness
    exceeds the sentinel cut-off (the ~1e6 m infinite-bedrock reservoir). The
    layer structure is global, so reading one partition is representative.
    """
    with h5py.File(stratal_path, "r") as f:
        H = np.asarray(f["stratH"])              # (nodes, nlayers)
    nlay = H.shape[1]
    lo = 0
    while lo < nlay and float(H[:, lo].max()) >= _SENTINEL_CUTOFF:
        lo += 1
    return lo


def build_partition(stratal_path, topology_path, lo, field, mesh_path=None):
    """
    Build the wedge volume for one partition / step.

    :arg stratal_path: ``stratal.<step>.p<p>.h5`` (per-layer fields).
    :arg topology_path: ``topology.p<p>.h5`` (``coords``, ``cells``).
    :arg lo: first layer to include (basement layers below are dropped).
    :arg field: ``"lithology"`` or ``"provenance"``.
    :arg mesh_path: optional ``<mesh>.<step>.p<p>.h5`` providing the current
        surface ``elev``. When given, the layer-interface elevations are
        reconstructed from the **current** surface minus the cumulative
        (compacted) thickness above each layer — so eroded layers (``stratH``
        ``= 0``) collapse to zero-thickness wedges instead of being drawn at
        their stale deposition elevation (which produced spurious cells above
        the eroded surface). Without it, the recorded ``stratZ`` is used.

    :return: dict with ``vertices`` (Nv, 3), ``wedges`` (Nw, 6), ``cells`` cell
        fields (each (Nw,)), ``frac`` (Nw, C) for provenance, and ``n_classes``.
        Returns ``None`` when the step has too few recorded layers to form a
        wedge.
    """
    with h5py.File(topology_path, "r") as f:
        coords = np.asarray(f["coords"], dtype=np.float64)       # (n, 3)
        tris = np.asarray(f["cells"], dtype=np.int64) - 1        # 0-indexed (n, 3)

    with h5py.File(stratal_path, "r") as f:
        stratZ = np.asarray(f["stratZ"], dtype=np.float64)       # (n, L)
        stratH = np.asarray(f["stratH"], dtype=np.float64)
        nlay = stratZ.shape[1]
        avail = set(f.keys())
        litho = field == "lithology"
        prov = field == "provenance"
        if litho:
            if "stratHf" not in avail or "phiF" not in avail:
                raise ValueError(
                    "field 'lithology' needs the dual-lithology 'stratHf'/'phiF' "
                    "fields, absent from %s (was the run dual lithology?)"
                    % stratal_path
                )
            stratHf = np.asarray(f["stratHf"], dtype=np.float64)
            phiS = np.asarray(f["phiS"], dtype=np.float64)
            phiF = np.asarray(f["phiF"], dtype=np.float64)
        if prov:
            if "stratP" not in avail:
                raise ValueError(
                    "field 'provenance' needs the 'stratP' field, absent from %s "
                    "(was the run provenance-enabled?)" % stratal_path
                )
            stratP = np.asarray(f["stratP"], dtype=np.float64)   # (n, L, C)

    n = coords.shape[0]
    m = tris.shape[0]
    nsurf = nlay - lo                                            # selected surfaces
    if nsurf < 2:
        return None                                              # no slab to draw
    nint = nsurf - 1                                             # wedge intervals

    # --- layer-interface elevations -------------------------------------------
    # Prefer the physically-current geometry: anchor at the present surface
    # (`elev`) and subtract the cumulative compacted thickness above each layer.
    # `erodeStrat` never updates `stratZ`, so an eroded layer keeps a stale
    # (higher) deposition elevation; using it would draw cells above the current
    # surface. The cumulative reconstruction collapses such layers (stratH == 0)
    # and stays monotonic. Falls back to the recorded `stratZ` if no mesh file.
    if mesh_path is not None and os.path.exists(mesh_path):
        with h5py.File(mesh_path, "r") as mf:
            topZ = np.asarray(mf["elev"], dtype=np.float64).reshape(-1)   # (n,)
        # above[:, l] = sum of thickness of the layers strictly above layer l.
        above = np.cumsum(stratH[:, ::-1], axis=1)[:, ::-1] - stratH
        surfZ = topZ[:, None] - above                            # top of layer l
    else:
        surfZ = stratZ

    # --- vertices: node replicated at each selected surface elevation ---------
    X = np.tile(coords[:, 0], nsurf)
    Y = np.tile(coords[:, 1], nsurf)
    Z = surfZ[:, lo:nlay].T.reshape(-1)                          # surface-major
    vertices = np.column_stack([X, Y, Z]).astype(np.float32)

    # --- wedges: triangle prism between consecutive surfaces ------------------
    si = np.arange(nint)
    bottom = tris[None, :, :] + (si[:, None, None] * n)          # (nint, m, 3)
    top = tris[None, :, :] + ((si + 1)[:, None, None] * n)
    wedges = np.concatenate([bottom, top], axis=2).reshape(-1, 6).astype(np.int32)

    # --- per-wedge (cell) fields = the slab's recorded layer property ---------
    # Interval i is the slab below surface i+1 == original layer (lo+i+1).
    lay = slice(lo + 1, nlay)                                    # layers lo+1..L-1

    def _cell(per_node_layer):
        # (n, nint) per-node-per-layer -> (Nw,) per-wedge, interval-major.
        return per_node_layer[tris].mean(axis=1).T.reshape(-1).astype(np.float32)

    out = {
        "vertices": vertices,
        "wedges": wedges,
        "cells": {},
        "frac": None,
        "n_classes": 0,
    }
    out["cells"]["thickness"] = _cell(stratH[:, lay])

    if litho:
        Hl = stratH[:, lay]
        Hf = stratHf[:, lay]
        with np.errstate(divide="ignore", invalid="ignore"):
            ff = np.where(Hl > 0.0, Hf / Hl, 0.0)
        out["cells"]["coarseThick"] = _cell(Hl - Hf)
        out["cells"]["fineThick"] = _cell(Hf)
        out["cells"]["fineFrac"] = _cell(ff)
        out["cells"]["phiCoarse"] = _cell(phiS[:, lay])
        out["cells"]["phiFine"] = _cell(phiF[:, lay])

    if prov:
        C = stratP.shape[2]
        P = stratP[:, lay, :]                                   # (n, nint, C)
        Hl = stratH[:, lay]
        with np.errstate(divide="ignore", invalid="ignore"):
            fnode = np.where(Hl[..., None] > 0.0, P / Hl[..., None], 0.0)
        # Per-wedge per-class fraction (interval-major to match `wedges`).
        fcell = fnode[tris].mean(axis=1)                        # (m, nint, C)
        frac = np.transpose(fcell, (1, 0, 2)).reshape(-1, C).astype(np.float32)
        tot = frac.sum(axis=1)
        dominant = np.where(tot > 0.0, frac.argmax(axis=1), -1).astype(np.int32)
        out["frac"] = frac
        out["n_classes"] = C
        out["cells"]["dominant"] = dominant

    return out


def write_partition_h5(path, data):
    """Write one partition's wedge volume to HDF5."""
    opts = {"compression": "gzip"}
    with h5py.File(path, "w") as f:
        f.create_dataset("vertices", data=data["vertices"], **opts)
        f.create_dataset("wedges", data=data["wedges"], **opts)
        for name, arr in data["cells"].items():
            f.create_dataset(name, data=arr, **opts)
        if data["frac"] is not None:
            f.create_dataset("frac", data=data["frac"], **opts)


def _attr_cell(name, h5file, path, nw, dtype="Float", prec=4):
    return (
        '         <Attribute Name="%s" AttributeType="Scalar" Center="Cell">\n'
        '          <DataItem Format="HDF" NumberType="%s" Precision="%d" '
        'Dimensions="%d 1">%s:/%s</DataItem>\n'
        "         </Attribute>\n" % (name, dtype, prec, nw, h5file, path)
    )


def _attr_frac_col(c, h5file, nw, ncl):
    # HyperSlab column c of the (nw, ncl) /frac array.
    return (
        '         <Attribute Name="src_class%d" AttributeType="Scalar" Center="Cell">\n'
        '          <DataItem ItemType="HyperSlab" Dimensions="%d 1">\n'
        '           <DataItem Dimensions="3 2" Format="XML">0 %d 1 1 %d 1</DataItem>\n'
        '           <DataItem Format="HDF" NumberType="Float" Precision="4" '
        'Dimensions="%d %d">%s:/frac</DataItem>\n'
        "          </DataItem>\n"
        "         </Attribute>\n" % (c, nw, c, nw, nw, ncl, h5file)
    )


def write_xdmf(out_prefix, meta, field, times=None):
    """
    Write a single self-contained ``<out>.xdmf`` (temporal collection): one time
    per step, each a spatial collection of the per-partition wedge grids.

    :arg meta: ``{step: [ {h5, nv, nw, cells:[names], ncl}, ... per partition ]}``.
    :arg times: optional ``{step: value}`` for the ParaView ``Time`` (e.g. the
        simulation year); defaults to the step index.
    """
    times = times or {}
    L = [
        '<?xml version="1.0" ?>\n<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n',
        '<Xdmf Version="2.0">\n <Domain>\n',
        '  <Grid Name="StrataVolume" GridType="Collection" CollectionType="Temporal">\n',
    ]
    for step in sorted(meta):
        L.append('   <Grid Name="step_%d" GridType="Collection" '
                 'CollectionType="Spatial">\n' % step)
        L.append('    <Time Value="%g"/>\n' % times.get(step, step))
        for blk in meta[step]:
            h5 = os.path.basename(blk["h5"])
            nv, nw = blk["nv"], blk["nw"]
            L.append('     <Grid Name="%s" GridType="Uniform">\n'
                     % blk["h5"].rsplit(".", 1)[0].rsplit(".", 1)[-1])
            L.append(
                '      <Topology TopologyType="Wedge" NumberOfElements="%d">\n'
                '       <DataItem Format="HDF" DataType="Int" Dimensions="%d 6">'
                "%s:/wedges</DataItem>\n"
                "      </Topology>\n" % (nw, nw, h5)
            )
            L.append(
                '      <Geometry GeometryType="XYZ">\n'
                '       <DataItem Format="HDF" NumberType="Float" Precision="4" '
                'Dimensions="%d 3">%s:/vertices</DataItem>\n'
                "      </Geometry>\n" % (nv, h5)
            )
            for name in blk["cells"]:
                dt = "Int" if name == "dominant" else "Float"
                L.append("   " + _attr_cell(name, h5, name, nw, dtype=dt))
            if field == "provenance":
                for c in range(blk["ncl"]):
                    L.append("   " + _attr_frac_col(c, h5, nw, blk["ncl"]))
            L.append("     </Grid>\n")
        L.append("   </Grid>\n")
    L.append("  </Grid>\n </Domain>\n</Xdmf>\n")
    with open(out_prefix + ".xdmf", "w") as f:
        f.writelines(L)


def main(argv=None):
    p = argparse.ArgumentParser(
        description="Build a 3-D wedge volume of the goSPL stratigraphy for "
        "ParaView (lithology or provenance)."
    )
    p.add_argument("--h5dir", required=True, help="goSPL output 'h5' directory")
    p.add_argument("--outdir", default="strata",
                   help="output directory, created if needed (default: strata)")
    p.add_argument("--out", default="strata",
                   help="output file prefix (default: strata)")
    p.add_argument("--field", default="lithology",
                   choices=["lithology", "provenance"],
                   help="cell field set to attach (default: lithology)")
    p.add_argument("--steps", default=None,
                   help="comma-separated steps (default: all found)")
    p.add_argument("--tout", type=float, default=None,
                   help="output interval in years; sets ParaView Time = "
                        "tstart + step*tout (default: use the step index)")
    p.add_argument("--tstart", type=float, default=0.0,
                   help="simulation start time in years for the --tout mapping")
    p.add_argument("--first-layer", type=int, default=None,
                   help="first layer index to include (default: auto-skip the "
                        "bedrock sentinel)")
    p.add_argument("--file-base", default="stratal",
                   help="stratal file base name (default: stratal)")
    p.add_argument("--mesh-base", default="gospl",
                   help="mesh-output file base providing the current surface "
                        "'elev' (default: gospl); used to place layer interfaces "
                        "by cumulative thickness so eroded layers don't overhang")
    args = p.parse_args(argv)

    if h5py is None:
        raise ImportError("h5py is required (pip install gospl[analysis])")

    # Rank 0 creates the output directory before any rank writes into it.
    if _RANK == 0:
        os.makedirs(args.outdir, exist_ok=True)
    if _COMM is not None and _SIZE > 1:
        _COMM.Barrier()

    all_steps, nparts = discover_steps(args.h5dir, args.file_base)
    steps = (
        [int(s) for s in args.steps.split(",")] if args.steps else all_steps
    )

    def sp(step, part):
        return os.path.join(args.h5dir, "%s.%d.p%d.h5" % (args.file_base, step, part))

    def tp(part):
        return os.path.join(args.h5dir, "topology.p%d.h5" % part)

    def mp(step, part):
        return os.path.join(
            args.h5dir, "%s.%d.p%d.h5" % (args.mesh_base, step, part)
        )

    local_meta = {}                      # step -> list of block dicts (this rank)
    for step in steps:
        lo = (
            args.first_layer
            if args.first_layer is not None
            else detect_first_layer(sp(step, 0))
        )
        blocks = []
        for part in range(nparts):
            if part % _SIZE != _RANK:    # split partitions across ranks
                continue
            data = build_partition(
                sp(step, part), tp(part), lo, args.field, mesh_path=mp(step, part)
            )
            if data is None:
                continue
            out_h5 = os.path.join(
                args.outdir, "%s.%d.p%d.h5" % (args.out, step, part)
            )
            write_partition_h5(out_h5, data)
            blocks.append({
                "h5": out_h5,
                "nv": len(data["vertices"]),
                "nw": len(data["wedges"]),
                "cells": list(data["cells"].keys()),
                "ncl": data["n_classes"],
            })
        if blocks:
            local_meta[step] = blocks

    # Gather metadata to rank 0 (which writes the single XDMF).
    if _COMM is not None and _SIZE > 1:
        gathered = _COMM.gather(local_meta, root=0)
    else:
        gathered = [local_meta]

    if _RANK == 0:
        meta = {}
        for part_meta in gathered:
            for step, blocks in part_meta.items():
                meta.setdefault(step, []).extend(blocks)
        for step in meta:                # stable partition order in the XDMF
            meta[step].sort(key=lambda b: b["h5"])
        # ParaView Time: the simulation year when --tout is given, else step.
        times = (
            {s: args.tstart + s * args.tout for s in meta}
            if args.tout is not None
            else None
        )
        out_prefix = os.path.join(args.outdir, args.out)
        write_xdmf(out_prefix, meta, args.field, times=times)
        print("wrote %s.xdmf (%d step(s), %d partition(s)) — open in ParaView"
              % (out_prefix, len(meta), nparts))
        skipped = [s for s in steps if s not in meta]
        if skipped:
            print("  (skipped step(s) %s — fewer than 2 recorded layers, so no "
                  "slab to draw)" % ", ".join(str(s) for s in skipped))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
