"""
Tests for the NetCDF grid / river-profile post-processor
(``gospl.analyse.gridexport``).

A tiny lattice mesh with a known tilted-plane elevation is written in goSPL's
output layout (mesh ``.npz`` + ``topology.p0.h5`` + ``gospl.0.p0.h5``), so the
rasterisation, raster D8 hydrology, basin/chi fields and river extraction are
exercised end-to-end without a model run.
"""

import os

import numpy as np
import pytest


def _synthetic_run(tmp_path, nx=10, ny=8, dx=1000.0):
    """Lattice mesh tilted down toward x=0; flow drains west."""
    h5py = pytest.importorskip("h5py")
    xs = np.arange(nx) * dx
    ys = np.arange(ny) * dx
    X, Y = np.meshgrid(xs, ys)
    v = np.column_stack([X.ravel(), Y.ravel(), np.zeros(nx * ny)])
    z = 0.02 * X.ravel() + 10.0 + 1e-3 * Y.ravel()        # west-low plane
    cells = []
    for j in range(ny - 1):
        for i in range(nx - 1):
            n0 = j * nx + i
            cells += [[n0, n0 + 1, n0 + nx + 1], [n0, n0 + nx + 1, n0 + nx]]
    cells = np.array(cells, dtype=np.int64)

    mesh = tmp_path / "mesh.npz"
    np.savez(str(mesh), v=v, c=cells, z=z)
    d = tmp_path / "h5"
    d.mkdir()
    with h5py.File(str(d / "topology.p0.h5"), "w") as f:
        f["coords"] = v
        f["cells"] = cells + 1                            # goSPL is 1-indexed
    with h5py.File(str(d / "gospl.0.p0.h5"), "w") as f:
        f["elev"] = z[:, None]
        f["erodep"] = np.zeros((nx * ny, 1))
        f["FA"] = np.ones((nx * ny, 1))
    return str(d), str(mesh), dx


def test_gridexport_hydrology(tmp_path):
    """grid_export produces basins/chi/area that respect the tilted surface."""
    gx = pytest.importorskip("gospl.analyse.gridexport")
    pytest.importorskip("matplotlib")
    pytest.importorskip("scipy")
    h5dir, mesh, dx = _synthetic_run(tmp_path)

    g = gx.grid_export(h5dir, mesh, spacing=dx)
    mask = g["mask"]
    # Sea level (hydrology base level) is exposed; defaults to 0 with no xmf.
    assert isinstance(g["base_level"], float) and g["base_level"] == 0.0
    assert gx.grid_export(h5dir, mesh, spacing=dx, base_level=3.5)["base_level"] == 3.5
    # Variable metadata table (units + definition), incl. flexure / rainfall.
    assert gx._VAR_META["elev"] == ("m", "surface elevation")
    assert gx._VAR_META["flexIso"][0] == "m" and "isostatic" in gx._VAR_META["flexIso"][1]
    assert gx._VAR_META["rain"][0] == "m/yr"
    assert g["elev"].shape == (g["y"].size, g["x"].size)
    # Drainage area grows toward the western (low) outlet.
    area = g["drainage_area"]
    assert np.nanmean(area[:, 1]) > np.nanmean(area[:, -2])
    # Every subaerial cell is assigned to a basin; chi finite and >= 0.
    assert np.all(g["basin"][mask] >= 0)
    assert np.isfinite(g["chi"][mask]).all()
    assert np.nanmin(g["chi"][mask]) >= -1e-9
    # chi increases away from the outlet (east interior > west).
    assert np.nanmean(g["chi"][:, -2]) > np.nanmean(g["chi"][:, 1])


def test_gridexport_basin_rivers_and_plots(tmp_path):
    """basin_rivers extracts a main stem + tributaries; the plots render."""
    gx = pytest.importorskip("gospl.analyse.gridexport")
    mpl = pytest.importorskip("matplotlib")
    mpl.use("Agg")
    h5dir, mesh, dx = _synthetic_run(tmp_path)

    g = gx.grid_export(h5dir, mesh, spacing=dx)
    assert "filled" in g                              # hydrologically-filled DEM
    riv = gx.basin_rivers(g, area_threshold=2.0 * dx * dx)
    assert riv["main_stem"] is not None
    ms = riv["main_stem"]
    for k in ("x", "y", "dist", "elev", "chi", "area", "elev_filled"):
        assert k in ms and len(ms[k]) >= 2
    # Main-stem elevation is monotonic along the channel (tilted plane).
    assert np.all(np.diff(np.sort(ms["elev"])) >= -1e-6)
    # The FILLED profile is strictly monotonic along the stem (packed
    # source -> outlet, so it decreases: each cell drains to a lower receiver),
    # unlike the raw elevation which can have small reversals.
    assert np.all(np.diff(ms["elev_filled"]) <= 1e-9)

    ax = gx.plot_long_profile(riv)                    # raw
    assert ax is not None
    axf = gx.plot_long_profile(riv, which="filled")   # conditioned (monotonic)
    assert axf is not None
    # plot_basin_map draws the sea-level coastline (proxy in the legend).
    ax2 = gx.plot_basin_map(g, riv, sea_level=50.0, figsize=(6, 4))
    assert ax2 is not None
    labels = [t.get_text() for t in ax2.get_legend().get_texts()]
    assert any("sea level" in s for s in labels)


def test_gridexport_geographic(tmp_path):
    """A spherical mesh is auto-detected and gridded in lon/lat (degrees)."""
    gx = pytest.importorskip("gospl.analyse.gridexport")
    pytest.importorskip("matplotlib")
    h5py = pytest.importorskip("h5py")
    R = 6.371e6
    nlon, nlat = 14, 12
    lons = np.radians(np.linspace(-20, 20, nlon))
    lats = np.radians(np.linspace(-20, 20, nlat))
    LON, LAT = np.meshgrid(lons, lats)
    X = R * np.cos(LAT) * np.cos(LON)
    Y = R * np.cos(LAT) * np.sin(LON)
    Z = R * np.sin(LAT)
    v = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    z = 500.0 * (LAT.ravel() - lats.min()) / (lats.max() - lats.min())  # rises N
    cells = []
    for j in range(nlat - 1):
        for i in range(nlon - 1):
            n0 = j * nlon + i
            cells += [[n0, n0 + 1, n0 + nlon + 1], [n0, n0 + nlon + 1, n0 + nlon]]
    cells = np.array(cells, dtype=np.int64)

    mesh = tmp_path / "sphere.npz"
    np.savez(str(mesh), v=v, c=cells, z=z)
    d = tmp_path / "h5"
    d.mkdir()
    with h5py.File(str(d / "topology.p0.h5"), "w") as f:
        f["coords"] = v
        f["cells"] = cells + 1
    with h5py.File(str(d / "gospl.0.p0.h5"), "w") as f:
        f["elev"] = z[:, None]

    g = gx.grid_export(str(d), str(mesh), spacing=4.0)     # 4-degree grid
    assert g["geographic"] is True
    assert -25.0 < g["x"].min() and g["x"].max() < 25.0    # lon in degrees
    assert -25.0 < g["y"].min() and g["y"].max() < 25.0    # lat in degrees
    assert np.any(g["basin"][g["mask"]] >= 0)              # basins delineated
    assert np.isfinite(g["chi"][g["mask"] & (g["basin"] >= 0)]).all()


def test_gridexport_netcdf(tmp_path):
    """to_netcdf writes a CF grid with the coordinate + field variables."""
    gx = pytest.importorskip("gospl.analyse.gridexport")
    nc = pytest.importorskip("netCDF4")
    h5dir, mesh, dx = _synthetic_run(tmp_path)

    g = gx.grid_export(h5dir, mesh, spacing=dx, base_level=5.0)
    assert g["base_level"] == 5.0                       # exposed in the result
    out = str(tmp_path / "surface.nc")
    gx.to_netcdf(g, out, time=1000.0)
    assert os.path.exists(out)
    with nc.Dataset(out) as ds:
        assert "x" in ds.variables and "y" in ds.variables
        for v in ("elev", "basin", "chi", "drainage_area"):
            assert v in ds.variables
        assert ds.variables["elev"].dimensions == ("y", "x")
        # Sea level recorded as a scalar variable + global attribute.
        assert "sea_level" in ds.variables
        assert float(ds.variables["sea_level"][...]) == 5.0
        assert float(ds.sea_level) == 5.0
        # Per-variable metadata (units + definition) is attached.
        assert ds.variables["elev"].units == "m"
        assert "drainage" in ds.variables["drainage_area"].long_name
        assert ds.variables["chi"].units == "m"
