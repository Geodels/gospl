"""
Tests for the stratigraphic section / well / Wheeler toolkit
(``gospl.analyse.stratasection``).

A small flat dual-lithology pile (with an eroded top layer = a hiatus) is
written in goSPL's output layout, so the reader, current-interface geometry,
colour fields and the four plot products are exercised without a model run.
"""

import numpy as np
import pytest


def _synthetic(tmp_path, nx=12, ny=10, dx=1000.0, L=5):
    h5py = pytest.importorskip("h5py")
    X, Y = np.meshgrid(np.arange(nx) * dx, np.arange(ny) * dx)
    n = nx * ny
    v = np.column_stack([X.ravel(), Y.ravel(), np.zeros(n)])
    cells = []
    for j in range(ny - 1):
        for i in range(nx - 1):
            a = j * nx + i
            cells += [[a, a + 1, a + nx + 1], [a, a + nx + 1, a + nx]]
    cells = np.array(cells, dtype=np.int64)
    elev = (0.01 * X - 50.0 + 1e-3 * Y).ravel()
    stratH = np.zeros((n, L))
    stratH[:, 0] = 1.0e6                                   # bedrock sentinel
    stratH[:, 1], stratH[:, 2], stratH[:, 3] = 20.0, 15.0, 25.0
    stratH[:, 4] = np.where(X.ravel() > 0.6 * X.max(), 0.0, 18.0)   # eroded east
    stratZ = np.zeros((n, L))
    for k in range(1, L):
        stratZ[:, k] = -40.0 + 10.0 * k + 1e-3 * Y.ravel()
    stratHf = np.zeros((n, L))
    for k in range(1, L):
        stratHf[:, k] = stratH[:, k] * (0.2 + 0.1 * k)
    d = tmp_path / "h5"
    d.mkdir()
    np.savez(str(tmp_path / "m.npz"), v=v, c=cells)
    with h5py.File(str(d / "topology.p0.h5"), "w") as f:
        f["coords"] = v
        f["cells"] = cells + 1
    with h5py.File(str(d / "stratal.0.p0.h5"), "w") as f:
        f["stratZ"] = stratZ
        f["stratH"] = stratH
        f["phiS"] = np.full((n, L), 0.49)
        f["stratHf"] = stratHf
        f["phiF"] = np.full((n, L), 0.63)
    with h5py.File(str(d / "gospl.0.p0.h5"), "w") as f:
        f["elev"] = elev[:, None]
    return str(d), str(tmp_path / "m.npz"), elev


def test_stratasection_load_and_interfaces(tmp_path):
    ss = pytest.importorskip("gospl.analyse.stratasection")
    pytest.importorskip("scipy")
    h5dir, mesh, elev = _synthetic(tmp_path)
    d = ss.load_strata(h5dir, mesh)
    assert d["nlayers"] == 4 and d["lo"] == 1          # sentinel dropped
    assert d["interfaces"].shape == (elev.size, 5)     # L+1
    assert np.allclose(d["interfaces"][:, -1], elev, equal_nan=True)   # surface
    # basement = surface - total recorded thickness
    tot = np.nansum(d["stratH"], axis=1)
    assert np.allclose(d["interfaces"][:, 0], elev - tot, equal_nan=True)


def test_stratasection_color_fields(tmp_path):
    ss = pytest.importorskip("gospl.analyse.stratasection")
    h5dir, mesh, _ = _synthetic(tmp_path)
    d = ss.load_strata(h5dir, mesh)
    for cb in ("deposition", "thickness", "lithology", "coarse", "porosity",
               "age"):
        arr, cmap, label = ss.color_field(d, cb)
        assert arr.shape == d["stratH"].shape
    # lithology fine fraction is in [0, 1] where present.
    fine, _, _ = ss.color_field(d, "lithology")
    fin = fine[np.isfinite(fine)]
    assert (fin >= -1e-9).all() and (fin <= 1 + 1e-9).all()


def test_stratasection_products_render(tmp_path):
    ss = pytest.importorskip("gospl.analyse.stratasection")
    mpl = pytest.importorskip("matplotlib")
    mpl.use("Agg")
    h5dir, mesh, _ = _synthetic(tmp_path)
    d = ss.load_strata(h5dir, mesh)

    ax = ss.cross_section(d, kind="x", color_by="lithology", vexag=20,
                          sea_level=0.0)
    assert ax.collections                                 # the layer mesh
    assert ss.cross_section(d, path=[[0, 0], [11000, 9000]],
                            color_by="thickness") is not None
    assert ss.horizontal_slice(d, z=-10.0, color_by="age") is not None
    assert ss.synthetic_well(d, 3000.0, 4000.0, color_by="lithology") is not None
    # Wheeler with the shoreline (sea-level) overlay.
    axw = ss.wheeler(d, kind="x", color_by="lithology", sea_level=0.0,
                     dt=5000.0)
    assert axw is not None


def test_stratasection_cli(tmp_path):
    ss = pytest.importorskip("gospl.analyse.stratasection")
    pytest.importorskip("matplotlib")
    h5dir, mesh, _ = _synthetic(tmp_path)
    out = str(tmp_path / "sec.pdf")
    rc = ss.main(["--h5dir", h5dir, "--mesh", mesh, "--kind", "wheeler",
                  "--color-by", "thickness", "--sea-level", "0", "--out", out])
    assert rc == 0
    import os
    assert os.path.exists(out)
