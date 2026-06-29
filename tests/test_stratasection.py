"""
Tests for the stratigraphic section / well / Wheeler toolkit
(``gospl.analyse.stratasection``).

A small flat dual-lithology pile (with an eroded top layer = a hiatus) is
written in goSPL's output layout, so the reader, current-interface geometry,
colour fields and the four plot products are exercised without a model run.
"""

import numpy as np
import pytest


def _synthetic(tmp_path, nx=12, ny=10, dx=1000.0, L=5, sea=None, time=None):
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
    if sea is not None or time is not None:
        # Minimal xmf carrying the sea level (a `sea` Function constant) and/or
        # the display Time the way goSPL encodes them, so _read_sealevel /
        # _read_time can recover them.
        xd = tmp_path / "xmf"
        xd.mkdir()
        txt = '<Grid>\n'
        if time is not None:
            txt += '  <Time Value="%r"/>\n' % float(time)
        if sea is not None:
            txt += ('  <Attribute Name="sea">\n'
                    '    <DataItem ItemType="Function" Function="$0 * 1.0 + %r">\n'
                    '    </DataItem>\n  </Attribute>\n' % float(sea))
        txt += '</Grid>\n'
        (xd / "gospl0.xmf").write_text(txt)
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


def test_stratasection_sealevel_from_sim(tmp_path):
    """load_strata reads the step's sea level; cross_section/wheeler default to
    it (drawing the datum line) without an explicit --sea-level."""
    ss = pytest.importorskip("gospl.analyse.stratasection")
    pytest.importorskip("scipy")
    pytest.importorskip("matplotlib")
    import numpy as np
    h5dir, mesh, _ = _synthetic(tmp_path, sea=-30.0)
    d = ss.load_strata(h5dir, mesh)
    assert d["sea_level"] == -30.0                         # read from the sim

    # cross_section with no explicit sea_level draws the datum at the sim value.
    ax = ss.cross_section(d, kind="x", color_by="thickness")
    assert any(np.allclose(line.get_ydata(), -30.0) for line in ax.lines)

    # No xmf -> sea_level is None and no datum line is drawn.
    nosea = tmp_path / "nosea"
    nosea.mkdir()
    h5dir2, mesh2, _ = _synthetic(nosea)
    d2 = ss.load_strata(h5dir2, mesh2)
    assert d2["sea_level"] is None
    ax2 = ss.cross_section(d2, kind="x", color_by="thickness")
    assert not any(np.allclose(line.get_ydata(), -30.0) for line in ax2.lines)


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


def test_stratasection_facies(tmp_path):
    """facies_field bins by deposition water depth; cross_section renders the
    discrete facies with a legend and true-elevation axis."""
    ss = pytest.importorskip("gospl.analyse.stratasection")
    pytest.importorskip("scipy")
    mpl = pytest.importorskip("matplotlib")
    mpl.use("Agg")
    import numpy as np
    h5dir, mesh, _ = _synthetic(tmp_path, sea=0.0)
    d = ss.load_strata(h5dir, mesh)

    # Deposition elevations in _synthetic span -30..0 (stratZ = -40+10k); with
    # sea level 0, depths 0..40 -> facies 0..2 present. Indices are integer in
    # [0, len(depths)] where layers exist; NaN where eroded.
    fac = ss.facies_field(d, sea_level=0.0)
    assert fac.shape == d["stratH"].shape
    fin = fac[np.isfinite(fac)]
    assert fin.min() >= 0 and fin.max() <= len(ss._FACIES_DEPTHS)
    assert np.allclose(fin, np.round(fin))                 # integer classes
    # Eroded layer (east column, layer 4) is NaN.
    assert np.isnan(fac[d["stratH"] <= 0]).all()

    # Custom thresholds change the classification.
    fac2 = ss.facies_field(d, sea_level=0.0, depths=(10.0,))
    assert np.nanmax(fac2) <= 1

    # Render: true elevation, custom figsize, layer lines every 2, vexag aspect.
    ax = ss.cross_section(d, kind="x", color_by="facies", sea_level=0.0,
                          figsize=(12, 3), layer_lines=2, vexag=10)
    assert ax.collections                                  # facies mesh
    assert ax.get_legend() is not None                     # facies legend
    assert ax.get_ylabel() == "elevation (m)"              # true elevation label
    assert ax.get_xlabel().endswith("(km)")                # distance shown in km
    assert ax.get_aspect() == 10                           # vexag applied as aspect
    assert ax.get_adjustable() == "datalim"                # so figsize is honoured
    assert tuple(ax.figure.get_size_inches()) == (12, 3)   # figsize honoured
    # Light-grey basement fill + solid sea-level line at the back.
    from matplotlib.collections import PolyCollection
    assert any(isinstance(c, PolyCollection) for c in ax.collections)
    sea = [ln for ln in ax.lines if ln.get_color() == "steelblue"]
    assert sea and sea[0].get_linestyle() == "-" and sea[0].get_zorder() < 1
    # Pinned xlim/ylim are honoured exactly; aspect left auto so figsize rules
    # (a locked vexag aspect would otherwise override the limits).
    axp = ss.cross_section(d, kind="x", color_by="thickness", vexag=50,
                           xlim=(0.0, 8000.0), ylim=(-100.0, 50.0))
    assert axp.get_xlim() == (0.0, 8000.0) and axp.get_ylim() == (-100.0, 50.0)
    assert axp.get_aspect() == "auto"

    # Wheeler also colours by facies (discrete legend) + honours figsize, and
    # the shoreline is in the legend (positionable via legend_loc).
    axw = ss.wheeler(d, kind="x", color_by="facies", sea_level=0.0,
                     figsize=(10, 5), legend_loc="upper left")
    assert axw.collections and axw.get_legend() is not None
    assert tuple(axw.figure.get_size_inches()) == (10, 5)
    assert axw.get_ylabel() == "deposition time (ky)"      # time on y-axis, ky
    assert axw.get_xlabel().endswith("(km)")               # distance in km
    lt = [t.get_text() for t in axw.get_legend().get_texts()]
    assert "shoreline" in lt                                # sea-level line shown
    # Shoreline position is computed at EVERY step (not only where facies are
    # present): the black trajectory line has one point per layer.
    sl = [ln for ln in axw.lines if ln.get_color() == "black"]
    assert sl and len(sl[0].get_xdata()) == d["nlayers"]
    # Settable horizontal-distance + time ranges.
    axx = ss.wheeler(d, kind="x", color_by="thickness", xlim=(1000.0, 5000.0),
                     ylim=(0.0, 1.0e5))
    assert axx.get_xlim() == (1000.0, 5000.0)
    assert axx.get_ylim() == (0.0, 1.0e5)

    # synthetic_well renders with a horizontal colour bar + custom title size.
    axwl = ss.synthetic_well(d, 3000.0, 4000.0, color_by="lithology",
                             title_fontsize=7, ylim=(-150.0, 80.0))
    assert axwl is not None
    assert axwl.title.get_fontsize() == 7
    assert axwl.get_ylim() == (-150.0, 80.0)               # settable elevation range

    # Multiple wells on one figure with a shared colour scale + single colour bar.
    fig, axes = ss.well_panel(d, [(3000.0, 4000.0), (6000.0, 4000.0),
                                  (9000.0, 4000.0)], color_by="porosity")
    assert len(axes) == 3
    # exactly one shared colour bar for the whole panel: 3 well axes + 1 cbar
    assert len(fig.axes) == 4
    # the wells share the elevation axis
    assert axes[0].get_ylim() == axes[1].get_ylim()
    # default titles are "well 1", "well 2", ...
    assert [a.get_title() for a in axes] == ["well 1", "well 2", "well 3"]

    # Default elevation range = union of the wells' interface elevations.
    locs = [(3000.0, 4000.0), (6000.0, 4000.0), (9000.0, 4000.0)]
    import numpy as np
    ifc = ss._interp_layers(d, d["interfaces"],
                            np.array([p[0] for p in locs]),
                            np.array([p[1] for p in locs]))
    fig2, ax2 = ss.well_panel(d, locs, color_by="porosity")
    lo, hi = ax2[0].get_ylim()
    assert abs(lo - float(np.nanmin(ifc))) < 1.0 and abs(hi - float(np.nanmax(ifc))) < 1.0
    # explicit vmin/vmax + ylim + custom labels are honoured
    fig3, ax3 = ss.well_panel(d, locs, color_by="porosity", vmin=0.0, vmax=0.6,
                              ylim=(-200.0, 100.0), labels=["A", "B", "C"])
    assert ax3[0].get_ylim() == (-200.0, 100.0)
    assert [a.get_title() for a in ax3] == ["A", "B", "C"]


def test_stratasection_time_from_sim(tmp_path):
    """The age / Wheeler time axis derives the layer interval from the step's
    stratal display time (the .xmf Time) when --strat-dt is not given."""
    ss = pytest.importorskip("gospl.analyse.stratasection")
    pytest.importorskip("scipy")
    pytest.importorskip("matplotlib")
    import numpy as np
    h5dir, mesh, _ = _synthetic(tmp_path, time=4.0e5)      # 4 layers, T=400 kyr
    d = ss.load_strata(h5dir, mesh)
    assert d["time"] == 4.0e5

    # age field spans real years up to the display time (surface = T).
    age, _, _ = ss.color_field(d, "age", dt=ss._effective_dt(d, 0.0, None))
    assert np.isclose(np.nanmax(age), 4.0e5)
    # Wheeler with no explicit dt uses the derived spacing -> y extent ~ T.
    axw = ss.wheeler(d, kind="x", color_by="thickness")
    assert axw.get_ylim()[1] >= 3.5e5


def test_stratasection_cli(tmp_path):
    ss = pytest.importorskip("gospl.analyse.stratasection")
    pytest.importorskip("matplotlib")
    h5dir, mesh, _ = _synthetic(tmp_path)
    out = str(tmp_path / "sec.pdf")
    rc = ss.main(["--h5dir", h5dir, "--mesh", mesh, "--kind", "wheeler",
                  "--color-by", "thickness", "--sea-level", "0", "--out", out])
    assert rc == 0


def test_stratasection_cli_figsize(tmp_path):
    """--figsize is honoured in the saved output (no forced tight bbox)."""
    ss = pytest.importorskip("gospl.analyse.stratasection")
    plt = pytest.importorskip("matplotlib.pyplot")
    import numpy as np
    h5dir, mesh, _ = _synthetic(tmp_path)
    out = str(tmp_path / "cs.png")
    rc = ss.main(["--h5dir", h5dir, "--mesh", mesh, "--kind", "cross",
                  "--color-by", "thickness", "--figsize", "8,3", "--out", out])
    assert rc == 0
    img = plt.imread(out)                      # (h, w, 4); dpi=200 -> 1600x600
    assert abs(img.shape[1] - 1600) <= 2 and abs(img.shape[0] - 600) <= 2
    import os
    assert os.path.exists(out)
