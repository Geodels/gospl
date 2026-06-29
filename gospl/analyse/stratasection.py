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
  borehole log coloured by any layer property (:func:`well_panel` draws several
  wells side by side on one figure with a shared colour scale).
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
                 "porosity", "age", "provenance", "facies")
_CMAP = {
    "deposition": "terrain", "thickness": "viridis", "lithology": "YlGnBu",
    "coarse": "YlOrBr", "porosity": "cividis", "age": "Spectral",
    "provenance": "tab20",
}
_LABEL = {
    "deposition": "deposition elevation (m)", "thickness": "layer thickness (m)",
    "lithology": "fine fraction", "coarse": "coarse fraction",
    "porosity": "bulk porosity", "age": "layer (oldest=0)",
    "provenance": "dominant source class", "facies": "depositional facies",
}
_SENTINEL = 1.0e5

# Depositional-facies defaults (tunable). The facies of a layer is set by the
# water depth at the time of deposition, depth = sea_level - deposition_elev
# (stratZ): subaerial (< 0) = fluvial/deltaic plain, then the increasing-depth
# marine facies. `_FACIES_DEPTHS` are the bin edges (m below sea level), giving
# len+1 facies; `_FACIES_COLORS`/`_FACIES_LABELS` map one entry per facies.
_FACIES_DEPTHS = (0.0, 20.0, 50.0, 75.0)
_FACIES_COLORS = ("limegreen", "darkkhaki", "sandybrown", "khaki", "c", "teal")
_FACIES_LABELS = ("fluvial / deltaic plain", "shoreface", "distal offshore",
                  "upper slope", "lower slope")


def facies_field(data, sea_level=None, depths=_FACIES_DEPTHS):
    """
    Per-(node, layer) **depositional-facies index** from the water depth at the
    time of deposition, ``depth = sea_level - stratZ`` (``stratZ`` is the
    recorded deposition elevation, i.e. the paleo-seabed). Binned by ``depths``
    (m below sea level), facies 0 is the shallowest (subaerial). Eroded/absent
    layers (``stratH <= 0``) are NaN.

    :arg sea_level: scalar, or a per-layer array (length ``nlayers``) for a
        time-varying paleo sea level; defaults to the step's sea level
        (``data['sea_level']``), else 0.
    :arg depths: increasing facies bin edges (m below sea level).
    """
    if sea_level is None:
        sea_level = data.get("sea_level")
        if sea_level is None:
            sea_level = 0.0
    sl = np.asarray(sea_level, dtype=float)
    if sl.ndim == 1:                                   # per-layer paleo sea level
        sl = sl[None, :]
    depth = sl - data["stratZ"]                        # >0 below sea level
    idx = np.digitize(depth, np.asarray(depths, dtype=float)).astype(float)
    idx[data["stratH"] <= 0] = np.nan
    return idx


def _facies_cmap_norm(depths, colors):
    """Discrete ListedColormap + BoundaryNorm for ``len(depths)+1`` facies."""
    from matplotlib.colors import ListedColormap, BoundaryNorm
    nfac = len(depths) + 1
    fcolors = list(colors)[:nfac]
    cmap = ListedColormap(fcolors)
    norm = BoundaryNorm(np.arange(-0.5, nfac + 0.5, 1.0), cmap.N)
    return cmap, norm, fcolors, nfac


def _facies_legend(ax, fcolors, labels, nfac, loc="lower left", extra=None):
    """Discrete facies legend (one patch per facies, plus any ``extra`` handles
    such as the shoreline line) on ``ax``."""
    from matplotlib.patches import Patch
    flabels = list(labels)[:nfac]
    handles = [Patch(facecolor=fcolors[i], label=flabels[i])
               for i in range(min(nfac, len(flabels)))]
    if extra:
        handles += list(extra)
    ax.legend(handles=handles, loc=loc, fontsize=8, framealpha=0.9,
              title="facies")


def _effective_dt(data, tstart, dt):
    """
    Resolve the per-layer time increment (yr). If ``dt`` is given, use it.
    Otherwise derive it from the step's **stratal display time**
    (``data['time']`` read from the .xmf) so the age / Wheeler time axis is the
    real simulation time: the surface layer maps to ``data['time']`` and the
    spacing is ``(time - tstart) / (lo + nlayers - 1)`` (the surface's global
    layer index). Falls back to 1.0 (unitless layer count) if no time is known.
    """
    if dt is not None:
        return dt
    T = data.get("time")
    if T is None:
        return 1.0
    denom = max(data["lo"] + data["nlayers"] - 1, 1)
    return (float(T) - tstart) / denom


def _km(ax, axis="x", geographic=False):
    """Display a distance axis in **km** (data stays in m). No-op for geographic
    (lon/lat) meshes, where the axis is in degrees."""
    if geographic:
        return "deg"
    from matplotlib.ticker import FuncFormatter
    fmt = FuncFormatter(lambda v, _: "%g" % (v / 1000.0))
    (ax.xaxis if axis == "x" else ax.yaxis).set_major_formatter(fmt)
    return "km"


def _ky(ax, axis="y"):
    """Display a time axis in **ky** (data stays in yr)."""
    from matplotlib.ticker import FuncFormatter
    fmt = FuncFormatter(lambda v, _: "%g" % (v / 1000.0))
    (ax.xaxis if axis == "x" else ax.yaxis).set_major_formatter(fmt)
    return "ky"


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
    # Sea level for this step, read from the simulation (the step's .xmf `sea`
    # constant); None if the file is absent. Used as the default sea level for
    # cross_section / wheeler so the shoreline datum matches the run without a
    # manual --sea-level.
    # Display time (years) of this step, also read from the .xmf, so the age /
    # Wheeler time axis is the real simulation time without a manual --strat-dt.
    try:
        from gospl.analyse.gridexport import _read_sealevel, _read_time
        out["sea_level"] = _read_sealevel(h5dir, file_base, step)
        out["time"] = _read_time(h5dir, file_base, step)
    except Exception:
        out["sea_level"] = None
        out["time"] = None
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
                  tstart=0.0, dt=None, figsize=None, layer_lines=0,
                  layer_line_kw=None, facies_depths=_FACIES_DEPTHS,
                  facies_colors=_FACIES_COLORS, facies_labels=_FACIES_LABELS,
                  legend_loc="lower left", title_fontsize=None, xlim=None,
                  ylim=None):
    """
    Vertical stratigraphic section along x / y / a path. Layers are drawn as a
    coloured ``pcolormesh`` (pinching where eroded) between the current
    interfaces, bounded by thick **surface** and **basement** lines, with a
    light-grey fill below the basement and a solid sea-level datum line drawn at
    the back (behind the coloured section).

    The y-axis shows the **true elevation** (m). Vertical exaggeration is applied
    as a real data aspect ratio (``vexag``), so the tick labels stay true; with
    ``vexag == 1`` the box is auto-scaled and ``figsize`` controls the look.

    :arg figsize: figure size (w, h) inches when a new figure is created.
    :arg vexag: vertical exaggeration (data aspect); 1 = auto-scaled (labels are
        always true elevation either way).
    :arg layer_lines: if > 0, overlay a thin interface line every N layers (e.g.
        2 = every other layer) on top of the surface/basement lines.
    :arg layer_line_kw: dict of matplotlib line kwargs for those layer lines
        (defaults: thin grey).
    :arg color_by: any of ``deposition``/``thickness``/``lithology``/``coarse``/
        ``porosity``/``age``/``provenance``/``facies``. ``facies`` classifies
        each layer by the water depth at deposition (see ``facies_field``) and
        draws a discrete colour + legend instead of a colour bar.
    :arg facies_depths / facies_colors / facies_labels: tunable facies bin edges
        (m below sea level), colours and legend labels.
    :arg sea_level: datum line + facies depth reference; defaults to the value
        read from the simulation for this step (``data["sea_level"]``). Omitted
        if the run did not record one.
    :arg xlim / ylim: explicit axis ranges ``(min, max)`` in data units (m). When
        either is given the view is clipped to exactly those bounds and the
        aspect is left **auto** (so ``figsize`` sets the proportions); ``vexag``
        only locks the data aspect when neither limit is pinned.
    """
    import matplotlib.pyplot as plt

    if sea_level is None:
        sea_level = data.get("sea_level")
    dt = _effective_dt(data, tstart, dt)
    facies = color_by == "facies"

    xs, ys, dist = _transect(data, kind, at, path, npts)
    iface = _interp_layers(data, data["interfaces"], xs, ys)        # (npts, L+1)
    if facies:
        cfld = facies_field(data, sea_level, facies_depths)
        label = _LABEL["facies"]
    else:
        cfld, dcmap, label = color_field(data, color_by, tstart, dt)
    cvals = _interp_layers(data, cfld, xs, ys)                      # (npts, L)

    Y = iface.T                                                    # true elevation
    X = np.broadcast_to(dist, Y.shape)
    # flat shading: C is one smaller than the corner mesh in each axis.
    C = np.ma.masked_invalid(cvals.T[:, :-1])                      # (L, npts-1)

    if ax is None:
        _, ax = plt.subplots(figsize=figsize or (9, 4))

    surf = iface[:, -1]
    base = iface[:, 0]

    # Sea-level datum: solid line at the BACK of the plot (low zorder) so the
    # coloured section is drawn over it and it only shows through gaps / above
    # the surface.
    if sea_level is not None and np.ndim(sea_level) == 0:
        ax.axhline(float(sea_level), color="steelblue", lw=1.0, zorder=0)

    # Light-grey basement: fill everything below the basement line. Drawn at the
    # back so the section sits on top of it; extend the y-limit down so it
    # reaches the bottom of the axes with no white gap.
    if np.isfinite(base).any():
        ylo = float(np.nanmin(base))
        rng = float(np.nanmax(surf) - ylo) if np.isfinite(surf).any() else 0.0
        ylo -= 0.03 * rng if rng > 0 else 1.0
        ax.fill_between(dist, ylo, base, color="0.85", zorder=0.5,
                        label="basement")

    if facies:
        fcmap, norm, fcolors, nfac = _facies_cmap_norm(facies_depths, facies_colors)
        # interpolated facies indices are continuous; snap to the nearest class.
        C = np.ma.masked_invalid(np.round(cvals.T[:, :-1]))
        mesh = ax.pcolormesh(X, Y, C, cmap=fcmap, norm=norm, shading="flat",
                             zorder=2)
        _facies_legend(ax, fcolors, facies_labels, nfac, loc=legend_loc)
    else:
        mesh = ax.pcolormesh(X, Y, C, cmap=cmap or dcmap, shading="flat",
                             zorder=2)
        ax.figure.colorbar(mesh, ax=ax, label=label, shrink=0.8)

    ax.plot(dist, surf, color="k", lw=2.2, zorder=5, label="surface")
    ax.plot(dist, base, color="0.25", lw=2.2, zorder=5)
    if layer_lines and layer_lines > 0:
        lkw = {"color": "0.4", "lw": 0.5, "alpha": 0.7, "zorder": 4}
        lkw.update(layer_line_kw or {})
        # interface columns are basement(0) .. surface(L); draw interior ones
        # every `layer_lines` layers (skip basement/surface, already drawn).
        for j in range(layer_lines, iface.shape[1] - 1, layer_lines):
            ax.plot(dist, iface[:, j], **lkw)

    # Axis limits / aspect. If the user pins xlim and/or ylim we honour them
    # EXACTLY (the view is clipped to those bounds) and leave the aspect auto so
    # `figsize` sets the proportions — a locked data aspect (`set_aspect`) would
    # otherwise re-adjust the limits at draw time and override them. With no
    # pinned limits, `vexag` locks a true vertical exaggeration (adjustable=
    # "datalim" keeps the box at `figsize`; the default "box" would resize it).
    pinned = xlim is not None or ylim is not None
    if vexag and vexag != 1 and not pinned:
        ax.set_aspect(vexag, adjustable="datalim")
    else:
        ax.set_aspect("auto")
    if xlim is not None:
        ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)
    u = _km(ax, "x", data["geographic"])
    ax.set_xlabel("distance along %s (%s)" % ("path" if path is not None else kind, u))
    ax.set_ylabel("elevation (m)")
    ax.set_title("Stratigraphic cross-section (%s)" % color_by,
                 fontsize=title_fontsize)
    return ax


def horizontal_slice(data, z, color_by="age", cmap=None, ax=None,
                     tstart=0.0, dt=None, figsize=None, title_fontsize=None):
    """
    Map at constant elevation ``z``: colour each location by the property of the
    layer preserved at that elevation (subcrop / structure slice). Cells where z
    is above the surface or below the basement are blank.
    """
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri

    dt = _effective_dt(data, tstart, dt)
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
        _, ax = plt.subplots(figsize=figsize or (7, 5))
    tri = mtri.Triangulation(data["x"], data["y"], data["cells"])
    tp = ax.tripcolor(tri, np.ma.masked_invalid(val), cmap=cmap or dcmap,
                      shading="gouraud")
    ax.figure.colorbar(tp, ax=ax, label=label, shrink=0.8)
    ax.set_aspect("equal")
    ux = _km(ax, "x", data["geographic"])
    _km(ax, "y", data["geographic"])
    ax.set_xlabel("x (%s)" % ux)
    ax.set_ylabel("y (%s)" % ux)
    ax.set_title("Layer %s preserved at z = %g m" % (color_by, z),
                 fontsize=title_fontsize)
    return ax


def synthetic_well(data, x, y, color_by="lithology", cmap=None, ax=None,
                   tstart=0.0, dt=None, width=1.0, figsize=None,
                   title_fontsize=None, cbar_orientation="horizontal",
                   vmin=None, vmax=None, colorbar=True, ylim=None):
    """
    Synthetic well at ``(x, y)``: the stratigraphic column drawn as stacked
    intervals (true current elevations) coloured by a layer property.

    To draw **several wells on one figure** with a shared colour scale, pass each
    an ``ax=`` (your own subplots) plus common ``vmin``/``vmax`` and
    ``colorbar=False``, then add a single colour bar — or just use
    :func:`well_panel`, which does this for a list of locations.

    :arg cbar_orientation: ``'horizontal'`` (default, colour bar across the base)
        or ``'vertical'``.
    :arg vmin / vmax: colour-scale limits; default is this well's own data range
        (pass shared values to make several wells comparable).
    :arg colorbar: draw the per-well colour bar (set ``False`` when sharing one).
    :arg ylim: elevation range ``(min, max)``; default is this well's own
        basement-to-surface span.
    :arg title_fontsize: font size for the title (matplotlib default if None).
    """
    import matplotlib.pyplot as plt
    from matplotlib import cm, colors

    dt = _effective_dt(data, tstart, dt)
    iface = _interp_layers(data, data["interfaces"],
                           np.array([x]), np.array([y]))[0]        # (L+1,)
    cfld, dcmap, label = color_field(data, color_by, tstart, dt)
    cval = _interp_layers(data, cfld, np.array([x]), np.array([y]))[0]  # (L,)

    finite = cval[np.isfinite(cval)]
    lo = vmin if vmin is not None else (finite.min() if finite.size else 0.0)
    hi = vmax if vmax is not None else (finite.max() if finite.size else 1.0)
    norm = colors.Normalize(lo, hi)
    cmp = plt.get_cmap(cmap or dcmap)

    if ax is None:
        _, ax = plt.subplots(figsize=figsize or (2.6, 6))
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
    if ylim is not None:
        ax.set_ylim(*ylim)
    else:
        ax.set_ylim(np.nanmin(iface), np.nanmax(iface))
    ax.set_xticks([])
    ax.set_ylabel("elevation (m)")
    if colorbar:
        sm = cm.ScalarMappable(norm=norm, cmap=cmp)
        if cbar_orientation == "horizontal":
            # A well is usually drawn tall & very narrow (e.g. figsize=(1, 8));
            # a default horizontal colour bar there is squished with overlapping
            # ticks. Use a thicker bar (low aspect), only a few ticks and a
            # small font so it stays readable.
            from matplotlib.ticker import MaxNLocator
            cb = ax.figure.colorbar(sm, ax=ax, location="bottom", fraction=0.08,
                                    pad=0.06, aspect=8)
            cb.ax.xaxis.set_major_locator(MaxNLocator(nbins=3))
            cb.ax.tick_params(labelsize=6)
            cb.set_label(label, fontsize=7)
        else:
            ax.figure.colorbar(sm, ax=ax, label=label, shrink=0.8)
    ax.set_title("Well (%.0f, %.0f)" % (x, y), fontsize=title_fontsize)
    return ax


def well_panel(data, locations, color_by="lithology", cmap=None, tstart=0.0,
               dt=None, width=1.0, figsize=None, title_fontsize=None,
               labels=None, sharey=True, vmin=None, vmax=None, ylim=None):
    """
    Draw **several synthetic wells side by side on one figure**, with a shared
    colour scale and a single colour bar at the base.

    :arg locations: list of ``(x, y)`` well positions.
    :arg labels: per-well titles; default ``"well 1"``, ``"well 2"``, … (pass a
        list to name them, one entry per location).
    :arg color_by: layer property coloured in every well (shared scale).
    :arg vmin / vmax: colour-scale limits; default is the field's finite range
        across all the wells.
    :arg ylim: elevation range ``(min, max)`` for every well; default is the
        **union of the wells' elevation bounds** (basement of the deepest to
        surface of the highest), so the panel spans the full range between the
        defined wells.
    :arg sharey: share the elevation axis across wells (default True).
    :return: ``(fig, axes)``.
    """
    import matplotlib.pyplot as plt
    from matplotlib import cm, colors

    dt = _effective_dt(data, tstart, dt)
    cfld, dcmap, label = color_field(data, color_by, tstart, dt)
    cmp = plt.get_cmap(cmap or dcmap)

    xs = np.array([p[0] for p in locations], dtype=float)
    ys = np.array([p[1] for p in locations], dtype=float)

    # Shared colour scale: the field's finite range across all wells, unless the
    # caller pins vmin/vmax.
    if vmin is None or vmax is None:
        vals = _interp_layers(data, cfld, xs, ys)                  # (nwell, L)
        fin = vals[np.isfinite(vals)]
        gmin, gmax = (float(fin.min()), float(fin.max())) if fin.size else (0.0, 1.0)
        vmin = gmin if vmin is None else vmin
        vmax = gmax if vmax is None else vmax

    # Default elevation range = union of all wells' interface elevations, so
    # every column is shown over the full elevation span between the wells.
    if ylim is None:
        ifc = _interp_layers(data, data["interfaces"], xs, ys)     # (nwell, L+1)
        if np.isfinite(ifc).any():
            ylim = (float(np.nanmin(ifc)), float(np.nanmax(ifc)))

    n = len(locations)
    fig, axes = plt.subplots(1, n, figsize=figsize or (1.4 * n + 0.6, 7),
                             sharey=sharey, squeeze=False)
    axes = axes[0]
    for i, (wx, wy) in enumerate(locations):
        synthetic_well(data, wx, wy, color_by=color_by, cmap=cmap, ax=axes[i],
                       tstart=tstart, dt=dt, width=width, vmin=vmin, vmax=vmax,
                       colorbar=False, title_fontsize=title_fontsize)
        # Default title "well N"; `labels` overrides (overwrites the per-well
        # "Well (x, y)" set by synthetic_well).
        ttl = labels[i] if labels is not None else "well %d" % (i + 1)
        axes[i].set_title(ttl, fontsize=title_fontsize)
        if ylim is not None:
            axes[i].set_ylim(*ylim)
        if i > 0 and sharey:
            axes[i].set_ylabel("")
    # one shared horizontal colour bar under all the wells
    sm = cm.ScalarMappable(norm=colors.Normalize(vmin, vmax), cmap=cmp)
    fig.colorbar(sm, ax=list(axes), location="bottom", fraction=0.05, pad=0.08,
                 aspect=40, label=label)
    return fig, axes


def wheeler(data, kind="x", at=None, path=None, color_by="thickness",
            npts=None, tstart=0.0, dt=None, sea_level=None, cmap=None, ax=None,
            figsize=None, facies_depths=_FACIES_DEPTHS,
            facies_colors=_FACIES_COLORS, facies_labels=_FACIES_LABELS,
            legend_loc="upper right", title_fontsize=None, xlim=None,
            ylim=None):
    """
    Wheeler (chronostratigraphic) diagram along a transect: distance (x) vs
    deposition **time** / layer (y), coloured by ``color_by``; **hiatuses and
    erosion are blank** (layers with zero preserved thickness along the
    transect).

    The **shoreline trajectory** is overlaid — the locus where the deposition
    water depth (``sea_level − stratZ``) crosses 0 (the subaerial↔marine
    boundary) — showing transgression/regression through time. ``sea_level``
    defaults to the value read from the simulation for this step
    (``data["sea_level"]``); pass a float, or a **per-layer array** (length
    ``nlayers``) for a time-varying paleo sea level, or it is omitted if the run
    did not record one.

    ``color_by="facies"`` colours each layer by the depositional facies (water
    depth at deposition; see :func:`facies_field`) with a discrete legend —
    consistent with the shoreline overlay. ``figsize`` and the ``facies_*`` lists
    are tunable.
    """
    import matplotlib.pyplot as plt

    if sea_level is None:
        sea_level = data.get("sea_level")
    dt = _effective_dt(data, tstart, dt)
    facies = color_by == "facies"

    xs, ys, dist = _transect(data, kind, at, path, npts)
    thick = _interp_layers(data, data["stratH"], xs, ys)          # (npts, L)
    if facies:
        cfld = facies_field(data, sea_level, facies_depths)
        label = _LABEL["facies"]
    else:
        cfld, dcmap, label = color_field(data, color_by, tstart, dt)
    cvals = _interp_layers(data, cfld, xs, ys)                    # (npts, L)
    C = np.where(np.nan_to_num(thick) > 1.0e-6, cvals, np.nan)    # gap = hiatus

    times = tstart + (data["lo"] + np.arange(data["nlayers"] + 1)) * dt
    X = np.broadcast_to(dist, (times.size, dist.size))
    Y = np.broadcast_to(times[:, None], X.shape)

    if ax is None:
        _, ax = plt.subplots(figsize=figsize or (9, 4))
    if facies:
        fcmap, norm, fcolors, nfac = _facies_cmap_norm(facies_depths, facies_colors)
        Cf = np.ma.masked_invalid(np.round(C.T[:, :-1]))
        mesh = ax.pcolormesh(X, Y, Cf, cmap=fcmap, norm=norm, shading="flat")
    else:
        mesh = ax.pcolormesh(X, Y, np.ma.masked_invalid(C.T[:, :-1]),
                             cmap=cmap or dcmap, shading="flat")
        ax.figure.colorbar(mesh, ax=ax, label=label, shrink=0.8)

    shoreline = None
    if sea_level is not None:
        # Shoreline position at EVERY time step: the distance along the transect
        # where the paleo-surface crosses sea level for each layer. Computed from
        # the RAW deposition elevation `stratZ` (recorded for every node and
        # layer, even where nothing was deposited), so it is drawn for every step
        # — not only where facies/deposits are present. One position per layer is
        # connected into a shoreline trajectory; it breaks only for steps whose
        # shoreline falls outside the transect.
        from matplotlib.lines import Line2D
        sZ = _interp_layers(data, data["stratZ"], xs, ys)         # (npts, L) raw
        sl = np.asarray(sea_level, dtype=float)
        depth = (sl[None, :] if sl.ndim == 1 else sl) - sZ        # >0 below sea
        tmid = 0.5 * (times[:-1] + times[1:])                     # layer centres
        shore = np.full(data["nlayers"], np.nan)
        for k in range(data["nlayers"]):
            dk = depth[:, k]
            ok = np.isfinite(dk)
            if ok.sum() < 2:
                continue
            di, xi = dk[ok], dist[ok]
            sub = di < 0.0                                        # subaerial side
            cr = np.where(sub[:-1] & ~sub[1:])[0]                 # land -> sea
            if cr.size == 0:
                cr = np.where(np.diff(np.sign(di)) != 0)[0]       # any crossing
            if cr.size:
                i = cr[0]
                d0, d1 = di[i], di[i + 1]
                f = d0 / (d0 - d1) if d1 != d0 else 0.0
                shore[k] = xi[i] + f * (xi[i + 1] - xi[i])
        ax.plot(shore, tmid, color="black", lw=2.0, zorder=6)
        shoreline = Line2D([], [], color="black", lw=2.0, label="shoreline")

    # Legend: facies patches (+ shoreline) for facies mode, else just the
    # shoreline. `legend_loc` is user-positionable.
    if facies:
        _facies_legend(ax, fcolors, facies_labels, nfac, loc=legend_loc,
                       extra=[shoreline] if shoreline is not None else None)
    elif shoreline is not None:
        ax.legend(handles=[shoreline], loc=legend_loc, fontsize=8)

    if xlim is not None:
        ax.set_xlim(*xlim)
    if ylim is not None:
        ax.set_ylim(*ylim)
    ux = _km(ax, "x", data["geographic"])
    ut = _ky(ax, "y")
    ax.set_xlabel("distance along %s (%s)" % ("path" if path is not None else kind, ux))
    ax.set_ylabel("deposition time (%s)" % ut)
    ax.set_title("Wheeler diagram (%s; gaps = hiatus/erosion)" % color_by,
                 fontsize=title_fontsize)
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
    p.add_argument("--vexag", type=float, default=1.0,
                   help="vertical exaggeration (true data aspect; labels stay "
                        "true elevation). 1 = auto-scaled")
    p.add_argument("--figsize", default=None,
                   help="figure size 'W,H' in inches (default 9,4)")
    p.add_argument("--layer-lines", type=int, default=0,
                   help="cross: overlay a thin interface line every N layers "
                        "(0 = off)")
    p.add_argument("--facies-depths", default=None,
                   help="facies bin edges (m below sea level), comma list "
                        "(default 0,20,50,75)")
    p.add_argument("--facies-colors", default=None,
                   help="facies colours, comma list (default "
                        "limegreen,darkkhaki,sandybrown,khaki,c,teal)")
    p.add_argument("--sea-level", type=float, default=None,
                   help="sea-level datum for the shoreline/section line and the "
                        "facies depth reference; default reads it from the "
                        "simulation (the step's xmf)")
    p.add_argument("--strat-dt", type=float, default=None,
                   help="layer interval (yr); default derives it from the "
                        "step's stratal display time (the .xmf Time)")
    p.add_argument("--tstart", type=float, default=0.0)
    p.add_argument("--legend-loc", default=None,
                   help="legend position (matplotlib loc, e.g. 'upper right', "
                        "'lower left')")
    p.add_argument("--xlim", default=None,
                   help="cross/wheeler: horizontal distance range 'MIN,MAX' "
                        "(same units as the transect, m; default: full extent)")
    p.add_argument("--ylim", default=None,
                   help="cross: elevation range 'MIN,MAX' (m); wheeler: time "
                        "range (yr). Pinning xlim/ylim leaves the aspect auto")
    p.add_argument("--title-fontsize", type=float, default=None,
                   help="title font size")
    p.add_argument("--tight", action="store_true",
                   help="crop the saved figure to its content (bbox_inches="
                        "'tight'); default keeps the requested --figsize")
    p.add_argument("--file-base", default="gospl")
    args = p.parse_args(argv)

    figsize = None
    if args.figsize:
        figsize = tuple(float(v) for v in args.figsize.split(","))
    fac_depths = _FACIES_DEPTHS
    if args.facies_depths:
        fac_depths = tuple(float(v) for v in args.facies_depths.split(","))
    fac_colors = _FACIES_COLORS
    if args.facies_colors:
        fac_colors = tuple(args.facies_colors.split(","))
    xlim = None
    if args.xlim:
        xlim = tuple(float(v) for v in args.xlim.split(","))
    ylim = None
    if args.ylim:
        ylim = tuple(float(v) for v in args.ylim.split(","))

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

    loc = args.legend_loc
    tf = args.title_fontsize
    if args.kind == "cross":
        ax = cross_section(data, kind=args.along, at=args.at, path=path,
                           color_by=args.color_by, vexag=args.vexag,
                           sea_level=args.sea_level, tstart=args.tstart,
                           dt=args.strat_dt, figsize=figsize,
                           layer_lines=args.layer_lines,
                           facies_depths=fac_depths, facies_colors=fac_colors,
                           legend_loc=loc or "lower left", title_fontsize=tf,
                           xlim=xlim, ylim=ylim)
    elif args.kind == "slice":
        if args.z is None:
            raise SystemExit("--z is required for --kind slice")
        ax = horizontal_slice(data, args.z, color_by=args.color_by,
                              tstart=args.tstart, dt=args.strat_dt,
                              figsize=figsize, title_fontsize=tf)
    elif args.kind == "well":
        if not args.xy:
            raise SystemExit("--xy 'x,y' is required for --kind well")
        wx, wy = (float(v) for v in args.xy.split(","))
        ax = synthetic_well(data, wx, wy, color_by=args.color_by,
                            tstart=args.tstart, dt=args.strat_dt,
                            figsize=figsize, title_fontsize=tf)
    else:
        ax = wheeler(data, kind=args.along, at=args.at, path=path,
                     color_by=args.color_by, tstart=args.tstart, dt=args.strat_dt,
                     sea_level=args.sea_level, figsize=figsize,
                     facies_depths=fac_depths, facies_colors=fac_colors,
                     legend_loc=loc or "upper right", title_fontsize=tf,
                     xlim=xlim, ylim=ylim)
    # Keep the requested --figsize by default; --tight crops to content.
    ax.figure.savefig(args.out, dpi=200,
                      bbox_inches="tight" if args.tight else None)
    print("wrote %s (%s; %d layers)" % (args.out, args.kind, data["nlayers"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
