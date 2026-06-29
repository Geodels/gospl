"""
Post-processing: extract, **per drainage basin**, the cell of maximum water
discharge and the cell of maximum sediment load from a gridded goSPL output —
i.e. each basin's **outflow point** (river mouth) and its flux.

This is the natural downstream step after :mod:`gospl.analyse.gridexport`: that
tool rasterises a goSPL surface to a CF-NetCDF grid carrying ``FA`` (flow
accumulation / water discharge), ``sedLoad`` (sediment load) and ``basin`` (the
drainage-basin id of every cell). Here, for every basin, we pick the single cell
where ``FA`` peaks (the water outlet) and the single cell where ``sedLoad`` peaks
(the sediment outlet), and write their ``lon``/``lat`` + value. Looped over a
time series, this gives the migrating river-mouth fluxes used in flux maps.

The water discharge defaults to ``FA``. Pass ``flow_var="fillFA"``
(``--flow-var fillFA``) to use the **depression-filled** accumulation instead,
which routes the trunk river *through* lakes / pits so a basin's outlet carries
its full upstream discharge even when the channel crosses a depression (raw
``FA`` can drop to zero inside a lake).

Consistent variable names
-------------------------
The defaults match the names :func:`gospl.analyse.gridexport.to_netcdf` writes,
so this tool reads a ``surface*.nc`` produced by ``gospl-grid`` **directly** — no
renaming step:

================  ====================  ===========================
quantity          gridexport name       legacy fallback (mapOutputs)
================  ====================  ===========================
water discharge   ``FA``                ``flowDischarge``
sediment load     ``sedLoad``           ``sedimentLoad``
basin id          ``basin``             ``basinID``
longitude         ``lon``               ``longitude``
latitude          ``lat``               ``latitude``
================  ====================  ===========================

Each variable falls back to the legacy name automatically if the gridexport one
is absent, so old ``fsdata*.nc`` files still work; explicit ``*_var`` arguments
override both.

Speed
-----
The per-basin maximum is a **grouped arg-max**, done once with a single
``lexsort`` over the subaerial cells (sort by ``(basin, value)``; the last cell
of each basin group is its maximum). This is :math:`O(N \\log N)` for the whole
grid, replacing the previous :math:`O(N_\\mathrm{basins} \\times N_\\mathrm{cells})`
scan (a full ``where(basin == k)`` per basin), which dominated runtime on a
global 0.1° grid (millions of cells × thousands of basins). It is fast enough to
run **serially** — the old MPI fan-out over basins is no longer needed.

Runnable as ``gospl-catchment`` (installed) or
``python -m gospl.analyse.catchment``::

    gospl-catchment -i inputSedFlow.csv -o flowsed

where the index CSV has two columns ``time,netcdf``::

    time,netcdf
    1,results/surface1.nc
    5,results/surface5.nc
    10,results/surface10.nc

writes ``flowsed/flow{time}.csv`` and ``flowsed/sed{time}.csv`` (columns
``basin,lon,lat,val``). In a notebook::

    from gospl.analyse.catchment import catchment_flux, basin_outflow
    out = catchment_flux("inputSedFlow.csv", "flowsed")   # batch -> CSVs
    o = basin_outflow("results/surface10.nc")             # one file
    flowdf, seddf = o["flow"], o["sed"]
"""

import os
import argparse

import numpy as np


# Variable-name aliases: the gridexport name first, the legacy mapOutputs name
# last. The first one present in the file is used (unless overridden). Water
# discharge defaults to the raw `FA`; `fillFA` (depression-filled accumulation)
# is kept as a fallback and can be selected explicitly with `--flow-var fillFA`.
_FLOW_ALIASES = ("FA", "fillFA", "flowDischarge")
_SED_ALIASES = ("sedLoad", "sedimentLoad")
_BASIN_ALIASES = ("basin", "basinID")
_LON_ALIASES = ("lon", "longitude")
_LAT_ALIASES = ("lat", "latitude")


def _resolve(present, requested, aliases, kind):
    """Pick the variable name to read: explicit request, else first alias present."""
    if requested is not None:
        if requested not in present:
            raise KeyError(
                "variable %r (%s) not in NetCDF; present: %s"
                % (requested, kind, sorted(present))
            )
        return requested
    for a in aliases:
        if a in present:
            return a
    raise KeyError(
        "no %s variable found (looked for %s) in NetCDF; present: %s"
        % (kind, list(aliases), sorted(present))
    )


def _read_grid(ncfile, flow_var, sed_var, basin_var, lon_var, lat_var):
    """
    Read the flow / sediment / basin fields and the lon-lat axes from a gridded
    NetCDF. ``ncfile`` may be a path, an open ``netCDF4.Dataset`` or an
    ``xarray.Dataset``. Returns ``(flow2d, sed2d, basin2d, lon1d, lat1d)``.
    """
    if isinstance(ncfile, (str, os.PathLike)):
        import netCDF4

        with netCDF4.Dataset(os.fspath(ncfile)) as ds:
            return _extract(ds, list(ds.variables), flow_var, sed_var,
                            basin_var, lon_var, lat_var)
    # already-open dataset: both netCDF4.Dataset and xarray.Dataset expose a
    # `.variables` mapping (xarray's includes coords) and index by name.
    return _extract(ncfile, list(ncfile.variables), flow_var, sed_var,
                    basin_var, lon_var, lat_var)


def _extract(ds, present, flow_var, sed_var, basin_var, lon_var, lat_var):
    def arr(name):
        v = ds[name]
        return np.asarray(v.values if hasattr(v, "values") else v[:])

    fv = _resolve(present, flow_var, _FLOW_ALIASES, "water discharge")
    sv = _resolve(present, sed_var, _SED_ALIASES, "sediment load")
    bv = _resolve(present, basin_var, _BASIN_ALIASES, "basin id")
    lov = _resolve(present, lon_var, _LON_ALIASES, "longitude")
    lav = _resolve(present, lat_var, _LAT_ALIASES, "latitude")
    return (arr(fv).astype(np.float64), arr(sv).astype(np.float64),
            arr(bv), arr(lov).astype(np.float64), arr(lav).astype(np.float64))


def _basin_argmax(basin, value, lon, lat, nb, min_cells):
    """
    For every basin id (``0..nb-1``), find the cell of MAXIMUM ``value`` and
    return its ``lon``/``lat``/``value``. Vectorised grouped arg-max: one stable
    ``lexsort`` by ``(basin, value)`` puts each basin's largest value last in its
    group, so the group boundaries give the per-basin maxima in one pass.

    Basins with at most ``min_cells`` cells are dropped (too small to be a
    meaningful catchment — matches the legacy ``len(ids) > 10`` filter). Returns
    a DataFrame ``basin,lon,lat,val`` (one row per kept basin).
    """
    import pandas as pd

    # Non-finite discharge can never be a maximum (NaN sorts last in lexsort, so
    # guard it explicitly); also lets us reject basins whose every cell is NaN.
    val = np.where(np.isfinite(value), value, -np.inf)
    order = np.lexsort((val, basin))            # basin asc, then value asc
    bs = basin[order]
    last = np.empty(bs.size, dtype=bool)
    last[-1] = True
    last[:-1] = bs[1:] != bs[:-1]               # last cell of each basin = its max
    sel = order[last]
    sb = basin[sel]

    counts = np.bincount(basin, minlength=nb)
    keep = (counts[sb] > min_cells) & np.isfinite(val[sel])
    sel, sb = sel[keep], sb[keep]

    df = pd.DataFrame({
        "basin": sb.astype(np.int64),
        "lon": lon[sel],
        "lat": lat[sel],
        "val": value[sel],
    })
    return df.sort_values("basin").reset_index(drop=True)


def basin_outflow(ncfile, min_cells=10, flow_var=None, sed_var=None,
                  basin_var=None, lon_var=None, lat_var=None):
    """
    Per-basin outflow points for one gridded NetCDF.

    :arg ncfile: path (or open ``netCDF4``/``xarray`` dataset) of a gridded
        surface — typically a ``gospl-grid`` ``surface*.nc``.
    :arg min_cells: basins with at most this many cells are ignored (default 10).
    :arg flow_var, sed_var, basin_var, lon_var, lat_var: override the variable
        names (defaults auto-detect the gridexport names, then the legacy ones).

    :return: ``{"flow": DataFrame, "sed": DataFrame}`` — each ``basin,lon,lat,val``
        (``val`` in m³/yr), one row per basin, the cell where that basin's water
        discharge (``flow``) / sediment load (``sed``) is largest.
    """
    flow, sed, basin, lon, lat = _read_grid(
        ncfile, flow_var, sed_var, basin_var, lon_var, lat_var)

    mlon, mlat = np.meshgrid(lon, lat)          # (nlat, nlon), matches the fields
    basin = basin.ravel()
    flow = flow.ravel()
    sed = sed.ravel()
    lonf = mlon.ravel()
    latf = mlat.ravel()

    # Subaerial cells carrying a basin id (marine / outside = -1 or NaN fill).
    bint = np.where(np.isfinite(basin.astype(np.float64)), basin, -1).astype(np.int64)
    sub = bint >= 0
    b = bint[sub]
    if b.size == 0:
        import pandas as pd
        empty = pd.DataFrame(columns=["basin", "lon", "lat", "val"])
        return {"flow": empty, "sed": empty.copy()}
    nb = int(b.max()) + 1

    flowdf = _basin_argmax(b, flow[sub], lonf[sub], latf[sub], nb, min_cells)
    seddf = _basin_argmax(b, sed[sub], lonf[sub], latf[sub], nb, min_cells)
    return {"flow": flowdf, "sed": seddf}


def catchment_flux(index, outdir=None, min_cells=10, verbose=True,
                   flow_var=None, sed_var=None, basin_var=None,
                   lon_var=None, lat_var=None):
    """
    Batch :func:`basin_outflow` over a time series of gridded NetCDFs.

    :arg index: a CSV path with columns ``time,netcdf`` (paths relative to the
        CSV's directory are resolved), a ``pandas.DataFrame`` with those columns,
        or an iterable of ``(time, ncfile)`` pairs.
    :arg outdir: if given, write ``outdir/flow{time}.csv`` and
        ``outdir/sed{time}.csv`` (columns ``basin,lon,lat,val``).

    :return: ``{time: {"flow": DataFrame, "sed": DataFrame}}``.
    """
    import pandas as pd
    from time import process_time

    base = ""
    if isinstance(index, (str, os.PathLike)):
        base = os.path.dirname(os.path.abspath(os.fspath(index)))
        rows = list(pd.read_csv(index)[["time", "netcdf"]].itertuples(
            index=False, name=None))
    elif isinstance(index, pd.DataFrame):
        rows = list(index[["time", "netcdf"]].itertuples(index=False, name=None))
    else:
        rows = list(index)

    if outdir is not None:
        os.makedirs(outdir, exist_ok=True)

    results = {}
    for time, ncfile in rows:
        t0 = process_time()
        path = os.fspath(ncfile)
        if base and not os.path.isabs(path) and not os.path.exists(path):
            cand = os.path.join(base, path)
            if os.path.exists(cand):
                path = cand
        if verbose:
            print("\nOpen output", path, flush=True)

        out = basin_outflow(path, min_cells=min_cells, flow_var=flow_var,
                            sed_var=sed_var, basin_var=basin_var,
                            lon_var=lon_var, lat_var=lat_var)
        results[time] = out

        if outdir is not None:
            step = int(time) if float(time).is_integer() else time
            out["flow"].to_csv(os.path.join(outdir, "flow%s.csv" % step),
                               index=False)
            out["sed"].to_csv(os.path.join(outdir, "sed%s.csv" % step),
                              index=False)
        if verbose:
            print("  +  %d flow / %d sediment outlets  (%.2f s)"
                  % (len(out["flow"]), len(out["sed"]), process_time() - t0),
                  flush=True)
    return results


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    p = argparse.ArgumentParser(
        description="Extract per-basin outflow (max water-discharge and "
        "max sediment-load cell) from gridded goSPL NetCDF output.",
        add_help=True,
    )
    p.add_argument("-i", "--input", required=True,
                   help="index CSV with columns time,netcdf")
    p.add_argument("-o", "--output", required=True,
                   help="output folder for flow{t}.csv / sed{t}.csv")
    p.add_argument("--min-cells", type=int, default=10,
                   help="ignore basins with <= this many cells (default 10)")
    p.add_argument("--flow-var", default=None,
                   help="water-discharge variable (default: FA; e.g. fillFA for "
                        "through-lake discharge; else flowDischarge)")
    p.add_argument("--sed-var", default=None,
                   help="sediment-load variable (default: sedLoad, else sedimentLoad)")
    p.add_argument("--basin-var", default=None,
                   help="basin-id variable (default: basin, else basinID)")
    p.add_argument("-q", "--quiet", action="store_true", help="no per-file log")
    args = p.parse_args(argv)

    catchment_flux(args.input, args.output, min_cells=args.min_cells,
                   verbose=not args.quiet, flow_var=args.flow_var,
                   sed_var=args.sed_var, basin_var=args.basin_var)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
