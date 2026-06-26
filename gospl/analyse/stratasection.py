"""
Post-processing: publication-quality **stratigraphic sections, wells and
Wheeler diagrams** from a goSPL run (matplotlib, vector PDF/SVG).

Reassembles the global stratigraphic pile (``stratZ``/``stratH`` and, when
present, ``stratHf``/``phiS``/``phiF`` and provenance ``stratP``) plus the
current surface ``elev``, drops the bedrock-sentinel layer, and builds the
**current** layer-interface elevations (``surface = elev``; each interface =
``elev − Σ thickness above``; ``basement = elev − Σ thickness``), so eroded
layers pinch out correctly. On that it offers four products:

* :func:`cross_section` — a vertical section along the **x**- or **y**-direction
  or an arbitrary **(x, y) path**: stratigraphy drawn as a coloured layer mesh
  (``color_by`` = deposition elevation / thickness / lithology / porosity / age
  / provenance) bounded by **thick surface and basement lines**, sampled at the
  mesh resolution.
* :func:`horizontal_slice` — a map at constant **elevation z**: which layer (and
  its property) is preserved at that depth (subcrop / structure slice).
* :func:`synthetic_well` — the full stratigraphic column at a point, as a
  borehole log coloured by any layer property.
* :func:`wheeler` — a chronostratigraphic (Wheeler) chart along a transect:
  distance vs deposition time, coloured by preserved thickness / composition,
  with **hiatuses and erosion shown as gaps**.

Each returns the matplotlib ``Figure``/``Axes`` for paper tuning. CLI
``gospl-section`` (or ``python -m gospl.analyse.stratasection``) with
``--kind cross|slice|well|wheeler``.
"""

import os
import argparse

import numpy as np


# Per-layer colour fields. Each maps to an (n_nodes, n_layers) array.
_COLOR_FIELDS = ("deposition", "thickness", "lithology", "coarse",
                 "porosity", "age", "provenance")
_CMAP = {
    "deposition": "terrain", "thickness": "viridis", "lithology": "YlGnBu",
    "coarse": "YlOrBr", "porosity": "cividis", "age": "Spectral",
    "provenance": "tab20",
}
_LABEL = {
    "deposition": "deposition elevation (m)", "thickness": "layer thickness (m)",
    "lithology": "fine fraction", "coarse": "coarse fraction",
    "porosity": "bulk porosity", "age": "layer (oldest=0)",
    "provenance": "dominant source class",
}
_SENTINEL = 1.0e5


def load_strata(h5dir, mesh, step=None, vkey="v", ckey="c", file_base="gospl",
                strat_base="stratal"):
    """
    Reassemble the global stratigraphic pile + current surface for one step.

    :return: a dict with the mesh (``coords``/``cells``/``x``/``y``,
        ``geographic``), ``elev`` (n,), per-layer ``stratZ``/``stratH`` and
        optional ``stratHf``/``phiS``/``phiF``/``stratP`` (n, L) with the bedrock
        sentinel dropped, ``interfaces`` (n, L+1) current interface elevations
        (basement..surface), ``nlayers`` and ``lo`` (first kept layer).
    """
    import glob
    import h5py
    from scipy.spatial import cKDTree

    data = np.load(mesh)
    coords = np.asarray(data[vkey], dtype=np.float64)
    cells = np.asarray(data[ckey], dtype=np.int64)
    npts = coords.shape[0]

    r = np.linalg.norm(coords, axis=1)
    geographic = bool(r.mean() > 1.0e5 and r.std() / max(r.mean(), 1.0) < 1.0e-3)
    if geographic:
        R = float(r.mean())
        x = np.degrees(np.arctan2(coords[:, 1], coords[:, 0]))
        y = np.degrees(np.arcsin(np.clip(coords[:, 2] / R, -1.0, 1.0)))
    else:
        x, y = coords[:, 0], coords[:, 1]

    if step is None:
        step = max(
            int(os.path.basename(f).split(".")[1])
            for f in glob.glob(os.path.join(h5dir, "%s.*.p*.h5" % strat_base))
        )
    tree = cKDTree(coords)

    # First pass: discover layer count C from the stratal file.
    sparts = sorted(glob.glob(os.path.join(h5dir, "%s.%d.p*.h5" % (strat_base, step))))
    if not sparts:
        raise FileNotFoundError("no %s.%d.p*.h5 in %s" % (strat_base, step, h5dir))
    with h5py.File(sparts[0], "r") as f:
        L = f["stratH"].shape[1]
        has_litho = "stratHf" in f and "phiF" in f
        has_prov = "stratP" in f
        nprov = f["stratP"].shape[2] if has_prov else 0

    fields = {k: np.full((npts, L), np.nan)
              for k in ("stratZ", "stratH", "phiS")}
    if has_litho:
        fields["stratHf"] = np.full((npts, L), np.nan)
        fields["phiF"] = np.full((npts, L), np.nan)
    prov = np.full((npts, L, nprov), np.nan) if has_prov else None
    elev = np.full(npts, np.nan)

    for pth in sparts:
        p = os.path.basename(pth).split(".p")[-1].split(".")[0]
        with h5py.File(os.path.join(h5dir, "topology.p%s.h5" % p), "r") as tf:
            lc = np.asarray(tf["coords"], dtype=np.float64)
        _, idx = tree.query(lc)
        with h5py.File(pth, "r") as f:
            for k in fields:
                if k in f:
                    fields[k][idx] = np.asarray(f[k])
            if has_prov:
                prov[idx] = np.asarray(f["stratP"])
        mpath = os.path.join(h5dir, "%s.%d.p%s.h5" % (file_base, step, p))
        if os.path.exists(mpath):
            with h5py.File(mpath, "r") as f:
                if "elev" in f:
                    elev[idx] = np.asarray(f["elev"])[:, 0]

    # Drop the bedrock sentinel (and any leading huge layer).
    lo = 0
    while lo < L and np.nanmax(fields["stratH"][:, lo]) >= _SENTINEL:
        lo += 1
    sl = slice(lo, L)
    stratH = np.nan_to_num(fields["stratH"][:, sl])
    out = {
        "coords": coords, "cells": cells, "x": x, "y": y,
        "geographic": geographic, "elev": elev,
        "stratZ": fields["stratZ"][:, sl], "stratH": stratH,
        "phiS": fields["phiS"][:, sl],
        "stratHf": fields["stratHf"][:, sl] if has_litho else None,
        "phiF": fields["phiF"][:, sl] if has_litho else None,
        "stratP": prov[:, sl, :] if has_prov else None,
        "nlayers": L - lo, "lo": lo, "step": step,
    }
    # Current interface elevations (basement .. surface).
    above = np.cumsum(stratH[:, ::-1], axis=1)[:, ::-1] - stratH   # above each layer
    top = elev[:, None] - above                                   # top of each layer
    out["interfaces"] = np.column_stack([top[:, :1] - stratH[:, :1], top])
    return out


def color_field(data, color_by, tstart=0.0, dt=1.0):
    """Per-(node, layer) array, colormap and label for a ``color_by`` choice."""
    H = data["stratH"]
    if color_by in ("lithology", "fine"):
        f = data["stratHf"]
        if f is None:
            raise ValueError("color_by='lithology' needs dual-lithology output")
        with np.errstate(divide="ignore", invalid="ignore"):
            arr = np.where(H > 0, f / H, np.nan)
    elif color_by == "coarse":
        f = data["stratHf"]
        if f is None:
            raise ValueError("color_by='coarse' needs dual-lithology output")
        with np.errstate(divide="ignore", invalid="ignore"):
            arr = np.where(H > 0, 1.0 - f / H, np.nan)
    elif color_by == "thickness":
        arr = np.where(H > 0, H, np.nan)
    elif color_by == "deposition":
        arr = data["stratZ"]
    elif color_by == "porosity":
        if data["phiF"] is not None:
            with np.errstate(divide="ignore", invalid="ignore"):
                ff = np.where(H > 0, data["stratHf"] / H, 0.0)
            arr = (1 - ff) * data["phiS"] + ff * data["phiF"]
        else:
            arr = data["phiS"]
    elif color_by == "age":
        arr = np.broadcast_to(
            tstart + (data["lo"] + np.arange(data["nlayers"])) * dt,
            H.shape,
        ).astype(float)
    elif color_by == "provenance":
        if data["stratP"] is None:
            raise ValueError("color_by='provenance' needs provenance output")
        arr = data["stratP"].argmax(axis=2).astype(float)
        arr[H <= 0] = np.nan
    else:
        raise ValueError("unknown color_by %r (choose from %s)"
                         % (color_by, ", ".join(_COLOR_FIELDS)))
    return arr, _CMAP.get(color_by, "viridis"), _LABEL.get(color_by, color_by)


def _transect(data, kind="x", at=None, path=None, npts=None):
    """
    Sample coordinates along a transect: ``kind='x'`` (constant y, spans x),
    ``kind='y'`` (constant x, spans y), or an explicit ``path`` of (x, y)
    waypoints. ``npts`` defaults to the mesh resolution. Returns (xs, ys, dist).
    """
    x, y = data["x"], data["y"]
    if npts is None:
        e = np.vstack([data["cells"][:, [0, 1]], data["cells"][:, [1, 2]]])
        res = float(np.median(np.linalg.norm(
            data["coords"][e[:, 0], :2] - data["coords"][e[:, 1], :2], axis=1)))
        if data["geographic"]:
            res = float(np.median(np.abs(np.diff(np.sort(np.unique(x))))) or 1.0)
    if path is not None:
        pts = np.asarray(path, dtype=float)
        seg = np.r_[0, np.cumsum(np.linalg.norm(np.diff(pts, axis=0), axis=1))]
        n = npts or max(2, int(seg[-1] / res) + 1)
        d = np.linspace(0, seg[-1], n)
        xs = np.interp(d, seg, pts[:, 0])
        ys = np.interp(d, seg, pts[:, 1])
        return xs, ys, d
    if kind == "x":
        at = 0.5 * (y.min() + y.max()) if at is None else at
        n = npts or max(2, int((x.max() - x.min()) / res) + 1)
        xs = np.linspace(x.min(), x.max(), n)
        ys = np.full(n, at)
    elif kind == "y":
        at = 0.5 * (x.min() + x.max()) if at is None else at
        n = npts or max(2, int((y.max() - y.min()) / res) + 1)
        ys = np.linspace(y.min(), y.max(), n)
        xs = np.full(n, at)
    else:
        raise ValueError("kind must be 'x', 'y', or pass path=")
    d = np.hypot(xs - xs[0], ys - ys[0])
    return xs, ys, d


def _interp_layers(data, arr, xs, ys):
    """Interpolate a per-(node, layer) array onto transect points -> (npts, L)."""
    import matplotlib.tri as mtri

    tri = mtri.Triangulation(data["x"], data["y"], data["cells"])
    out = np.full((xs.size, arr.shape[1]), np.nan)
    for k in range(arr.shape[1]):
        out[:, k] = mtri.LinearTriInterpolator(tri, arr[:, k])(xs, ys)
    return out


def cross_section(data, kind="x", at=None, path=None, color_by="lithology",
                  npts=None, vexag=1.0, sea_level=None, cmap=None, ax=None,
                  tstart=0.0, dt=1.0):
    """
    Vertical stratigraphic section along x / y / a path. Layers are drawn as a
    coloured ``pcolormesh`` (pinching where eroded) between the current
    interfaces, bounded by thick **surface** and **basement** lines.
    """
    import matplotlib.pyplot as plt

    xs, ys, dist = _transect(data, kind, at, path, npts)
    iface = _interp_layers(data, data["interfaces"], xs, ys)        # (npts, L+1)
    cfld, dcmap, label = color_field(data, color_by, tstart, dt)
    cvals = _interp_layers(data, cfld, xs, ys)                      # (npts, L)

    Y = (iface * vexag).T                                          # (L+1, npts)
    X = np.broadcast_to(dist, Y.shape)
    # flat shading: C is one smaller than the corner mesh in each axis.
    C = np.ma.masked_invalid(cvals.T[:, :-1])                      # (L, npts-1)

    if ax is None:
        _, ax = plt.subplots(figsize=(9, 4))
    mesh = ax.pcolormesh(X, Y, C, cmap=cmap or dcmap, shading="flat")
    surf = iface[:, -1] * vexag
    base = iface[:, 0] * vexag
    ax.plot(dist, surf, color="k", lw=2.2, zorder=5, label="surface")
    ax.plot(dist, base, color="0.25", lw=2.2, zorder=5, label="basement")
    if sea_level is not None:
        ax.axhline(sea_level * vexag, color="steelblue", lw=1.0, ls="--",
                   zorder=4)
    ax.figure.colorbar(mesh, ax=ax, label=label, shrink=0.8)
    ax.set_xlabel("distance along %s" % ("path" if path is not None else kind))
    ax.set_ylabel("elevation (m)%s" % ("" if vexag == 1 else " x%g" % vexag))
    ax.set_title("Stratigraphic cross-section (%s)" % color_by)
    return ax


def horizontal_slice(data, z, color_by="age", cmap=None, ax=None,
                     tstart=0.0, dt=1.0):
    """
    Map at constant elevation ``z``: colour each location by the property of the
    layer preserved at that elevation (subcrop / structure slice). Cells where z
    is above the surface or below the basement are blank.
    """
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri

    iface = data["interfaces"]                                    # (n, L+1)
    cfld, dcmap, label = color_field(data, color_by, tstart, dt)  # (n, L)
    # Which layer interval contains z, per node.
    below = iface <= z                                            # (n, L+1)
    k = below.sum(axis=1) - 1                                     # layer index
    inside = (k >= 0) & (k < data["nlayers"])
    val = np.full(data["x"].shape, np.nan)
    rows = np.where(inside)[0]
    val[rows] = cfld[rows, k[rows]]

    if ax is None:
        _, ax = plt.subplots(figsize=(7, 5))
    tri = mtri.Triangulation(data["x"], data["y"], data["cells"])
    tp = ax.tripcolor(tri, np.ma.masked_invalid(val), cmap=cmap or dcmap,
                      shading="gouraud")
    ax.figure.colorbar(tp, ax=ax, label=label, shrink=0.8)
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Layer %s preserved at z = %g m" % (color_by, z))
    return ax


def synthetic_well(data, x, y, color_by="lithology", cmap=None, ax=None,
                   tstart=0.0, dt=1.0, width=1.0):
    """
    Synthetic well at ``(x, y)``: the stratigraphic column drawn as stacked
    intervals (true current elevations) coloured by a layer property.
    """
    import matplotlib.pyplot as plt
    from matplotlib import cm, colors

    iface = _interp_layers(data, data["interfaces"],
                           np.array([x]), np.array([y]))[0]        # (L+1,)
    cfld, dcmap, label = color_field(data, color_by, tstart, dt)
    cval = _interp_layers(data, cfld, np.array([x]), np.array([y]))[0]  # (L,)

    finite = cval[np.isfinite(cval)]
    norm = colors.Normalize(finite.min(), finite.max()) if finite.size \
        else colors.Normalize(0, 1)
    cmp = plt.get_cmap(cmap or dcmap)

    if ax is None:
        _, ax = plt.subplots(figsize=(2.6, 6))
    for k in range(data["nlayers"]):
        z0, z1 = iface[k], iface[k + 1]
        if not np.isfinite(z0) or not np.isfinite(z1) or z1 - z0 <= 0:
            continue
        col = cmp(norm(cval[k])) if np.isfinite(cval[k]) else "white"
        ax.add_patch(plt.Rectangle((0, z0), width, z1 - z0, facecolor=col,
                                   edgecolor="0.5", lw=0.3))
    ax.axhline(iface[-1], color="k", lw=2.0)                      # surface
    ax.axhline(iface[0], color="0.25", lw=2.0)                    # basement
    ax.set_xlim(0, width)
    ax.set_ylim(np.nanmin(iface), np.nanmax(iface))
    ax.set_xticks([])
    ax.set_ylabel("elevation (m)")
    ax.figure.colorbar(cm.ScalarMappable(norm=norm, cmap=cmp), ax=ax,
                       label=label, shrink=0.8)
    ax.set_title("Well (%.0f, %.0f)" % (x, y))
    return ax


def wheeler(data, kind="x", at=None, path=None, color_by="thickness",
            npts=None, tstart=0.0, dt=1.0, sea_level=None, cmap=None, ax=None):
    """
    Wheeler (chronostratigraphic) diagram along a transect: distance (x) vs
    deposition time / layer (y), coloured by ``color_by``; **hiatuses and
    erosion are blank** (layers with zero preserved thickness along the
    transect). If ``sea_level`` is given, the **shoreline trajectory** is
    overlaid — the locus where each layer's deposition elevation (``stratZ``)
    crosses sea level — showing transgression/regression through time.
    """
    import matplotlib.pyplot as plt

    xs, ys, dist = _transect(data, kind, at, path, npts)
    thick = _interp_layers(data, data["stratH"], xs, ys)          # (npts, L)
    cfld, dcmap, label = color_field(data, color_by, tstart, dt)
    cvals = _interp_layers(data, cfld, xs, ys)                    # (npts, L)
    C = np.where(np.nan_to_num(thick) > 1.0e-6, cvals, np.nan)    # gap = hiatus

    times = tstart + (data["lo"] + np.arange(data["nlayers"] + 1)) * dt
    X = np.broadcast_to(dist, (times.size, dist.size))
    Y = np.broadcast_to(times[:, None], X.shape)

    if ax is None:
        _, ax = plt.subplots(figsize=(9, 4))
    mesh = ax.pcolormesh(X, Y, np.ma.masked_invalid(C.T[:, :-1]),
                         cmap=cmap or dcmap, shading="flat")
    if sea_level is not None:
        depoz = _interp_layers(data, data["stratZ"], xs, ys)      # (npts, L)
        tmid = 0.5 * (times[:-1] + times[1:])                     # layer centres
        ax.contour(dist, tmid, depoz.T, levels=[sea_level],
                   colors="steelblue", linewidths=1.6)
        ax.plot([], [], color="steelblue", lw=1.6, label="shoreline (sea level)")
        ax.legend(loc="upper right", fontsize=8)
    ax.figure.colorbar(mesh, ax=ax, label=label, shrink=0.8)
    ax.set_xlabel("distance along %s" % ("path" if path is not None else kind))
    ax.set_ylabel("deposition time" + ("" if dt == 1 else " (yr)"))
    ax.set_title("Wheeler diagram (%s; gaps = hiatus/erosion)" % color_by)
    return ax


def main(argv=None):
    p = argparse.ArgumentParser(
        description="Stratigraphic sections / wells / Wheeler diagrams from a "
        "goSPL run (matplotlib, vector output).")
    p.add_argument("--h5dir", required=True, help="goSPL output 'h5' directory")
    p.add_argument("--mesh", required=True,
                   help="global mesh FILE.npz[:vkey[:ckey]] (default keys v,c)")
    p.add_argument("--kind", default="cross",
                   choices=["cross", "slice", "well", "wheeler"])
    p.add_argument("--out", default="section.pdf", help="output figure (vector)")
    p.add_argument("--step", type=int, default=None)
    p.add_argument("--color-by", default="lithology", choices=list(_COLOR_FIELDS))
    p.add_argument("--along", default="x", help="cross/wheeler: 'x', 'y'")
    p.add_argument("--at", type=float, default=None, help="transect position")
    p.add_argument("--path", default=None,
                   help="x0,y0;x1,y1;... polyline (overrides --along)")
    p.add_argument("--z", type=float, default=None, help="slice: elevation")
    p.add_argument("--xy", default=None, help="well: 'x,y'")
    p.add_argument("--vexag", type=float, default=1.0)
    p.add_argument("--sea-level", type=float, default=None)
    p.add_argument("--strat-dt", type=float, default=1.0, help="layer interval (yr)")
    p.add_argument("--tstart", type=float, default=0.0)
    p.add_argument("--file-base", default="gospl")
    args = p.parse_args(argv)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    parts = args.mesh.split(":")
    mesh, vkey, ckey = (parts + ["v", "c"])[:1] + \
        [parts[1] if len(parts) > 1 else "v", parts[2] if len(parts) > 2 else "c"]
    data = load_strata(args.h5dir, mesh, step=args.step, vkey=vkey, ckey=ckey,
                       file_base=args.file_base)
    path = None
    if args.path:
        path = [[float(v) for v in pt.split(",")] for pt in args.path.split(";")]

    if args.kind == "cross":
        ax = cross_section(data, kind=args.along, at=args.at, path=path,
                           color_by=args.color_by, vexag=args.vexag,
                           sea_level=args.sea_level, tstart=args.tstart,
                           dt=args.strat_dt)
    elif args.kind == "slice":
        if args.z is None:
            raise SystemExit("--z is required for --kind slice")
        ax = horizontal_slice(data, args.z, color_by=args.color_by,
                              tstart=args.tstart, dt=args.strat_dt)
    elif args.kind == "well":
        if not args.xy:
            raise SystemExit("--xy 'x,y' is required for --kind well")
        wx, wy = (float(v) for v in args.xy.split(","))
        ax = synthetic_well(data, wx, wy, color_by=args.color_by,
                            tstart=args.tstart, dt=args.strat_dt)
    else:
        ax = wheeler(data, kind=args.along, at=args.at, path=path,
                     color_by=args.color_by, tstart=args.tstart, dt=args.strat_dt,
                     sea_level=args.sea_level)
    ax.figure.savefig(args.out, bbox_inches="tight", dpi=200)
    print("wrote %s (%s; %d layers)" % (args.out, args.kind, data["nlayers"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
