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


def test_flow_var_fillfa_opt_in(tmp_path):
    """Water flux defaults to FA; flow_var='fillFA' opts into the filled field."""
    netCDF4 = pytest.importorskip("netCDF4")
    cm = pytest.importorskip("gospl.analyse.catchment")
    lon = np.array([0.0, 1.0, 2.0, 3.0])
    lat = np.array([0.0, 1.0, 2.0, 3.0])
    basin = np.zeros((4, 4), dtype=np.int32)
    fa = np.zeros((4, 4))
    fillfa = np.zeros((4, 4))
    fa[0, 0] = 100.0          # raw-FA max here
    fillfa[3, 3] = 999.0      # filled-FA (trunk-through-lake) max elsewhere
    p = tmp_path / "grid.nc"
    with netCDF4.Dataset(str(p), "w") as ds:
        ds.createDimension("lat", 4)
        ds.createDimension("lon", 4)
        ds.createVariable("lon", "f8", ("lon",))[:] = lon
        ds.createVariable("lat", "f8", ("lat",))[:] = lat
        ds.createVariable("FA", "f8", ("lat", "lon"))[:, :] = fa
        ds.createVariable("fillFA", "f8", ("lat", "lon"))[:, :] = fillfa
        ds.createVariable("sedLoad", "f8", ("lat", "lon"))[:, :] = np.zeros((4, 4))
        ds.createVariable("basin", "i4", ("lat", "lon"))[:, :] = basin

    # Default uses FA -> outlet at the raw-FA max (0,0), value 100.
    flow = cm.basin_outflow(str(p), min_cells=0)["flow"].set_index("basin")
    assert (flow.loc[0, "lon"], flow.loc[0, "lat"], flow.loc[0, "val"]) == (0.0, 0.0, 100.0)
    # Opt in to fillFA -> outlet at the filled-FA max (3,3), value 999.
    flow_fill = cm.basin_outflow(str(p), min_cells=0, flow_var="fillFA")["flow"].set_index("basin")
    assert (flow_fill.loc[0, "lon"], flow_fill.loc[0, "lat"], flow_fill.loc[0, "val"]) == (3.0, 3.0, 999.0)


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
    # No solute field in this grid -> no solute output, and no "solute" key.
    assert not (outdir / "solute3.csv").exists()
    assert "solute" not in res[3]


def _write_solute_grid(path):
    """4x4 two-basin grid + riverSolute total and two per-species fields."""
    netCDF4 = pytest.importorskip("netCDF4")
    lon = np.array([0.0, 1.0, 2.0, 3.0])
    lat = np.array([0.0, 1.0, 2.0, 3.0])
    basin = np.array([
        [0, 0, 1, 1],
        [0, 0, 1, 1],
        [0, 0, 1, 1],
        [-1, -1, -1, -1],
    ], dtype=np.int32)
    fa = np.zeros((4, 4))
    sed = np.zeros((4, 4))
    riv = np.zeros((4, 4))
    carb = np.zeros((4, 4))
    sil = np.zeros((4, 4))
    # Basin 0 solute outlet at (lon=1, lat=0): total 30 = 20 carbonate + 10 silica.
    riv[0, 1], carb[0, 1], sil[0, 1] = 30.0, 20.0, 10.0
    # Basin 1 solute outlet at (lon=2, lat=2): total 50 = 5 carbonate + 45 silica.
    riv[2, 2], carb[2, 2], sil[2, 2] = 50.0, 5.0, 45.0
    with netCDF4.Dataset(str(path), "w") as ds:
        ds.createDimension("lat", 4)
        ds.createDimension("lon", 4)
        ds.createVariable("lon", "f8", ("lon",))[:] = lon
        ds.createVariable("lat", "f8", ("lat",))[:] = lat
        for nm, g in (("FA", fa), ("sedLoad", sed), ("riverSolute", riv),
                      ("riverSolute_carbonate", carb), ("riverSolute_silica", sil)):
            ds.createVariable(nm, "f8", ("lat", "lon"))[:, :] = g
        ds.createVariable("basin", "i4", ("lat", "lon"))[:, :] = basin
    return str(path)


def test_basin_solute_flux_per_species(tmp_path):
    """
    basin_solute_flux picks each basin's solute outlet (max total flux) and reports
    the total plus each species' flux AT that cell; per-species columns sum to val.
    """
    pytest.importorskip("netCDF4")
    cm = pytest.importorskip("gospl.analyse.catchment")
    ncf = _write_solute_grid(tmp_path / "sol.nc")

    df = cm.basin_solute_flux(ncf, min_cells=0).set_index("basin")
    assert list(df.columns) == ["lon", "lat", "val", "carbonate", "silica"]
    # Basin 0 outlet + total + species breakdown.
    assert (df.loc[0, "lon"], df.loc[0, "lat"], df.loc[0, "val"]) == (1.0, 0.0, 30.0)
    assert (df.loc[0, "carbonate"], df.loc[0, "silica"]) == (20.0, 10.0)
    # Basin 1 outlet elsewhere, silica-dominated.
    assert (df.loc[1, "lon"], df.loc[1, "lat"], df.loc[1, "val"]) == (2.0, 2.0, 50.0)
    assert (df.loc[1, "carbonate"], df.loc[1, "silica"]) == (5.0, 45.0)
    # Species columns sum to the total flux at each outlet.
    assert np.allclose(df["carbonate"] + df["silica"], df["val"])


def test_basin_solute_flux_absent_and_alias(tmp_path):
    """No solute field -> basin_outflow omits 'solute' (KeyError swallowed);
    with soluteflux (no riverSolute) the fallback alias is used."""
    netCDF4 = pytest.importorskip("netCDF4")
    cm = pytest.importorskip("gospl.analyse.catchment")

    # (a) plain flow/sed grid: basin_outflow has no 'solute'; direct call raises.
    plain = _synthetic(tmp_path, _GRIDEXPORT)
    assert "solute" not in cm.basin_outflow(plain, min_cells=0)
    with pytest.raises(KeyError):
        cm.basin_solute_flux(plain, min_cells=0)

    # (b) only `soluteflux` present (no riverSolute) -> alias fallback works.
    lon = np.array([0.0, 1.0, 2.0, 3.0])
    lat = np.array([0.0, 1.0, 2.0, 3.0])
    basin = np.zeros((4, 4), dtype=np.int32)
    sf = np.zeros((4, 4)); sf[2, 2] = 12.0
    p = tmp_path / "sf.nc"
    with netCDF4.Dataset(str(p), "w") as ds:
        ds.createDimension("lat", 4); ds.createDimension("lon", 4)
        ds.createVariable("lon", "f8", ("lon",))[:] = lon
        ds.createVariable("lat", "f8", ("lat",))[:] = lat
        ds.createVariable("soluteflux", "f8", ("lat", "lon"))[:, :] = sf
        ds.createVariable("basin", "i4", ("lat", "lon"))[:, :] = basin
    df = cm.basin_solute_flux(str(p), min_cells=0).set_index("basin")
    assert (df.loc[0, "lon"], df.loc[0, "lat"], df.loc[0, "val"]) == (2.0, 2.0, 12.0)
    # No per-species fields present -> only the base columns.
    assert list(df.columns) == ["lon", "lat", "val"]
