"""
Tests for the per-basin outflow extractor (``gospl.analyse.catchment``).

A tiny synthetic gridded NetCDF (lon/lat axes + FA/sedLoad/basin fields) with a
known maximum cell per basin is written, so the grouped arg-max, the variable
alias fallback, the min-cells filter and the CSV batch path are exercised
without a model run.
"""

import numpy as np
import pytest


def _write_grid(path, names, basin, fa, sed, lon, lat):
    """Write a gridded NetCDF using the given (lon, lat, fa, sed, basin) names."""
    netCDF4 = pytest.importorskip("netCDF4")
    lonn, latn, fan, sedn, basinn = names
    with netCDF4.Dataset(str(path), "w") as ds:
        ds.createDimension(latn, lat.size)
        ds.createDimension(lonn, lon.size)
        ds.createVariable(lonn, "f8", (lonn,))[:] = lon
        ds.createVariable(latn, "f8", (latn,))[:] = lat
        ds.createVariable(fan, "f8", (latn, lonn))[:, :] = fa
        ds.createVariable(sedn, "f8", (latn, lonn))[:, :] = sed
        ds.createVariable(basinn, "i4", (latn, lonn))[:, :] = basin


def _synthetic(tmp_path, names):
    """
    4x4 grid, two basins (0 left half, 1 right half) + a marine strip (-1).
    The FA / sedLoad maxima sit at known, *different* cells per basin.
    """
    lon = np.array([0.0, 1.0, 2.0, 3.0])
    lat = np.array([0.0, 1.0, 2.0, 3.0])
    basin = np.array([
        [0, 0, 1, 1],
        [0, 0, 1, 1],
        [0, 0, 1, 1],
        [-1, -1, -1, -1],     # marine row (excluded)
    ], dtype=np.int32)
    fa = np.zeros((4, 4))
    sed = np.zeros((4, 4))
    fa[1, 0] = 500.0          # basin 0 water max at (lon=0, lat=1)
    sed[2, 1] = 90.0          # basin 0 sed   max at (lon=1, lat=2)
    fa[0, 3] = 800.0          # basin 1 water max at (lon=3, lat=0)
    sed[1, 2] = 70.0          # basin 1 sed   max at (lon=2, lat=1)
    p = tmp_path / "grid.nc"
    _write_grid(p, names, basin, fa, sed, lon, lat)
    return str(p)


_GRIDEXPORT = ("lon", "lat", "FA", "sedLoad", "basin")
_LEGACY = ("longitude", "latitude", "flowDischarge", "sedimentLoad", "basinID")


@pytest.mark.parametrize("names", [_GRIDEXPORT, _LEGACY])
def test_basin_outflow_maxima_and_aliases(tmp_path, names):
    """Per-basin outflow picks the right max cell — for both name conventions."""
    cm = pytest.importorskip("gospl.analyse.catchment")
    ncf = _synthetic(tmp_path, names)

    out = cm.basin_outflow(ncf, min_cells=0)   # keep the small synthetic basins
    flow = out["flow"].set_index("basin")
    sed = out["sed"].set_index("basin")

    # Water-discharge maxima (lon, lat, val).
    assert (flow.loc[0, "lon"], flow.loc[0, "lat"], flow.loc[0, "val"]) == (0.0, 1.0, 500.0)
    assert (flow.loc[1, "lon"], flow.loc[1, "lat"], flow.loc[1, "val"]) == (3.0, 0.0, 800.0)
    # Sediment-load maxima.
    assert (sed.loc[0, "lon"], sed.loc[0, "lat"], sed.loc[0, "val"]) == (1.0, 2.0, 90.0)
    assert (sed.loc[1, "lon"], sed.loc[1, "lat"], sed.loc[1, "val"]) == (2.0, 1.0, 70.0)
    # Marine cells (basin -1) never appear.
    assert set(flow.index) == {0, 1}


def test_min_cells_filter(tmp_path):
    """Basins with <= min_cells cells are dropped (each synthetic basin has 6)."""
    cm = pytest.importorskip("gospl.analyse.catchment")
    ncf = _synthetic(tmp_path, _GRIDEXPORT)
    assert len(cm.basin_outflow(ncf, min_cells=5)["flow"]) == 2   # 6 > 5: kept
    assert len(cm.basin_outflow(ncf, min_cells=6)["flow"]) == 0   # 6 > 6 false: dropped


def test_catchment_flux_batch_csv(tmp_path):
    """catchment_flux writes flow{t}.csv / sed{t}.csv with the expected columns."""
    pd = pytest.importorskip("pandas")
    cm = pytest.importorskip("gospl.analyse.catchment")
    ncf = _synthetic(tmp_path, _GRIDEXPORT)

    index = tmp_path / "index.csv"
    pd.DataFrame({"time": [3], "netcdf": [ncf]}).to_csv(index, index=False)
    outdir = tmp_path / "flowsed"
    res = cm.catchment_flux(str(index), str(outdir), min_cells=0, verbose=False)

    assert set(res.keys()) == {3}
    fcsv = outdir / "flow3.csv"
    scsv = outdir / "sed3.csv"
    assert fcsv.exists() and scsv.exists()
    df = pd.read_csv(fcsv)
    assert list(df.columns) == ["basin", "lon", "lat", "val"]
    assert df.set_index("basin").loc[1, "val"] == 800.0
