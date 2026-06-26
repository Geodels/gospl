"""
Tests for the stratigraphy wedge-volume post-processor
(``gospl.analyse.stratamesh``).

Synthetic, fast and model-free: tiny ``stratal``/``topology`` HDF5 files are
written by hand (mirroring goSPL's output layout) so the geometry, connectivity,
field reduction and XDMF wiring are exercised without a full run.
"""

import os
import xml.etree.ElementTree as ET

import numpy as np
import pytest


def _write_dataset(tmp_path, nparts=2):
    """
    Write a tiny multi-partition stratal dataset: 4 nodes / 2 triangles per
    partition, 3 layers where layer 0 is the bedrock sentinel (1e6 m thick).
    Dual-lithology + provenance fields included.
    """
    h5py = pytest.importorskip("h5py")
    d = tmp_path / "h5"
    d.mkdir()
    n, m, L, C = 4, 2, 3, 2
    coords = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]], float)
    cells = np.array([[1, 2, 3], [2, 4, 3]], int)        # 1-indexed (gospl)

    stratZ = np.zeros((n, L))
    stratZ[:, 1] = -50.0
    stratZ[:, 2] = -20.0
    stratH = np.zeros((n, L))
    stratH[:, 0] = 1.0e6                                 # sentinel
    stratH[:, 1] = 50.0
    stratH[:, 2] = 30.0
    stratHf = np.zeros((n, L))
    stratHf[:, 1] = 20.0                                 # fine frac 0.4
    stratHf[:, 2] = 9.0                                  # fine frac 0.3
    phiS = np.full((n, L), 0.49)
    phiF = np.full((n, L), 0.63)
    stratP = np.zeros((n, L, C))
    stratP[:, :, 0] = stratH * 0.6
    stratP[:, :, 1] = stratH * 0.4

    for p in range(nparts):
        with h5py.File(str(d / ("topology.p%d.h5" % p)), "w") as f:
            f["coords"] = coords
            f["cells"] = cells
        with h5py.File(str(d / ("stratal.0.p%d.h5" % p)), "w") as f:
            f["stratZ"] = stratZ
            f["stratH"] = stratH
            f["phiS"] = phiS
            f["stratHf"] = stratHf
            f["phiF"] = phiF
            f["stratP"] = stratP
    return str(d), n, m, L, C, nparts


def test_stratamesh_basement_and_geometry(tmp_path):
    """detect_first_layer skips the sentinel; build_partition makes the wedge."""
    sm = pytest.importorskip("gospl.analyse.stratamesh")
    h5dir, n, m, L, C, _ = _write_dataset(tmp_path, nparts=1)

    lo = sm.detect_first_layer(os.path.join(h5dir, "stratal.0.p0.h5"))
    assert lo == 1, "the 1e6 m bedrock sentinel (layer 0) must be skipped"

    d = sm.build_partition(
        os.path.join(h5dir, "stratal.0.p0.h5"),
        os.path.join(h5dir, "topology.p0.h5"),
        lo, "lithology",
    )
    nsurf = L - lo                                        # 2 surfaces
    assert d["vertices"].shape == (n * nsurf, 3)
    assert d["wedges"].shape == (m * (nsurf - 1), 6)
    # Each wedge is 3 bottom-surface nodes + the same 3 one surface up (+n).
    bottom, top = d["wedges"][:, :3], d["wedges"][:, 3:]
    assert np.array_equal(top, bottom + n)
    # Vertex Z is the per-layer deposition elevation (surface-major): the first
    # selected surface (layer lo=1) is the wedge bottom, the next its top.
    assert np.allclose(d["vertices"][0:n, 2], -50.0)         # surface lo=1 (bottom)
    assert np.allclose(d["vertices"][n : 2 * n, 2], -20.0)   # surface lo+1=2 (top)
    # Cell fields reduce the slab (layer 2): thickness 30, fine fraction 9/30.
    assert np.allclose(d["cells"]["thickness"], 30.0)
    assert np.allclose(d["cells"]["fineFrac"], 0.3)
    assert np.allclose(d["cells"]["coarseThick"], 21.0)
    assert np.allclose(d["cells"]["phiCoarse"], 0.49)


def test_stratamesh_provenance_fields(tmp_path):
    """Provenance mode yields per-class fractions + a dominant-source field."""
    sm = pytest.importorskip("gospl.analyse.stratamesh")
    h5dir, n, m, L, C, _ = _write_dataset(tmp_path, nparts=1)
    d = sm.build_partition(
        os.path.join(h5dir, "stratal.0.p0.h5"),
        os.path.join(h5dir, "topology.p0.h5"),
        1, "provenance",
    )
    assert d["n_classes"] == C
    assert d["frac"].shape == (m, C)                     # 1 interval, m wedges
    assert np.allclose(d["frac"], [0.6, 0.4])
    assert np.all(d["cells"]["dominant"] == 0)           # class 0 dominates
    # Per-layer porosity is attached alongside the provenance composition.
    assert np.allclose(d["cells"]["porosity"], 0.49)
    assert np.allclose(d["cells"]["phiFine"], 0.63)      # dual run -> phiF too
    # Per-cell thickness + elevation (mid-height of the wedge: layer 2 spans
    # surfaces stratZ -50 .. -20 -> mid -35; thickness 30).
    assert np.allclose(d["cells"]["thickness"], 30.0)
    assert np.allclose(d["cells"]["elevation"], -35.0)
    # Per-cell recorded-layer index (only layer 2 here).
    assert np.all(d["cells"]["layer"] == 2)


def test_stratamesh_provenance_inherits_underlying(tmp_path):
    """
    An eroded / pinched-out layer (stratH == 0 -> no source composition) takes
    the composition of the cell directly below it in the column, so it is never
    a no-source cell (dominant == -1).
    """
    sm = pytest.importorskip("gospl.analyse.stratamesh")
    h5py = pytest.importorskip("h5py")
    d = tmp_path
    coords = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]], float)
    cells = np.array([[1, 2, 3], [2, 4, 3]], int)            # 1-indexed
    n, C = 4, 2
    L = 5                                                    # sentinel + 4 layers
    stratZ = np.zeros((n, L))
    stratZ[:, 1] = -80.0
    stratZ[:, 2] = -60.0
    stratZ[:, 3] = -40.0
    stratZ[:, 4] = -20.0
    stratH = np.zeros((n, L))
    stratH[:, 0] = 1.0e6                                     # basal sentinel
    stratH[:, 1:4] = 20.0                                    # layers 1-3 present
    stratH[:, 4] = 0.0                                       # top layer eroded
    # Wedge intervals are layers lo+1..L-1 = 2,3,4. layer 2 -> class 0,
    # layer 3 -> class 1, layer 4 eroded -> should inherit layer 3's class 1.
    stratP = np.zeros((n, L, C))
    stratP[:, 1, 0] = stratH[:, 1]
    stratP[:, 2, 0] = stratH[:, 2]                           # interval 0 -> class 0
    stratP[:, 3, 1] = stratH[:, 3]                           # interval 1 -> class 1
    with h5py.File(str(d / "topology.p0.h5"), "w") as f:
        f["coords"] = coords
        f["cells"] = cells
    with h5py.File(str(d / "stratal.0.p0.h5"), "w") as f:
        f["stratZ"] = stratZ
        f["stratH"] = stratH
        f["phiS"] = np.full((n, L), 0.49)
        f["stratP"] = stratP

    out = sm.build_partition(
        str(d / "stratal.0.p0.h5"), str(d / "topology.p0.h5"), 1, "provenance"
    )
    m = cells.shape[0]
    dom = out["cells"]["dominant"]
    assert dom.shape == (3 * m,)                             # 3 intervals
    assert not np.any(dom == -1)                             # no no-source cells
    assert np.all(dom[0:m] == 0)                             # interval 0 = class 0
    assert np.all(dom[m : 2 * m] == 1)                       # interval 1 = class 1
    # Top interval (eroded layer 4) inherits the underlying layer 3 -> class 1.
    assert np.all(dom[2 * m : 3 * m] == 1)
    # Per-cell recorded-layer index: intervals map to original layers 2, 3, 4.
    lyr = out["cells"]["layer"]
    assert np.all(lyr[0:m] == 2) and np.all(lyr[m : 2 * m] == 3)
    assert np.all(lyr[2 * m : 3 * m] == 4)


def test_stratamesh_eroded_no_overhang(tmp_path):
    """
    With a mesh ``elev`` file, an eroded top layer (``stratH == 0`` but a stale,
    higher ``stratZ``) must NOT produce cells above the current surface: the
    interfaces are rebuilt from the current surface minus cumulative thickness,
    so the eroded layer collapses. The recorded-``stratZ`` fallback shows the
    overhang the fix removes.
    """
    sm = pytest.importorskip("gospl.analyse.stratamesh")
    h5py = pytest.importorskip("h5py")
    d = tmp_path / "h5"
    d.mkdir()
    n, m, L = 4, 2, 4
    coords = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]], float)
    cells = np.array([[1, 2, 3], [2, 4, 3]], int)
    # layer 0 sentinel; 1,2 real; 3 eroded (H=0) but with a stale stratZ = +100.
    stratZ = np.zeros((n, L))
    stratZ[:, 1] = -50.0
    stratZ[:, 2] = -20.0
    stratZ[:, 3] = 100.0                                 # stale (eroded away)
    stratH = np.zeros((n, L))
    stratH[:, 0] = 1.0e6
    stratH[:, 1] = 50.0
    stratH[:, 2] = 30.0
    stratH[:, 3] = 0.0                                   # eroded
    with h5py.File(str(d / "topology.p0.h5"), "w") as f:
        f["coords"], f["cells"] = coords, cells
    with h5py.File(str(d / "stratal.0.p0.h5"), "w") as f:
        f["stratZ"], f["stratH"] = stratZ, stratH
        f["phiS"] = np.full((n, L), 0.49)
        f["stratHf"] = stratH * 0.4
        f["phiF"] = np.full((n, L), 0.63)
    with h5py.File(str(d / "gospl.0.p0.h5"), "w") as f:
        f["elev"] = np.full((n, 1), -20.0)               # current surface

    sp = str(d / "stratal.0.p0.h5")
    tp = str(d / "topology.p0.h5")
    mesh = str(d / "gospl.0.p0.h5")

    # Legacy (no mesh): the stale +100 surface is drawn -> overhang above -20.
    legacy = sm.build_partition(sp, tp, 1, "lithology")
    assert legacy["vertices"][:, 2].max() > 90.0

    # elev-anchored: nothing above the current surface; eroded layer collapsed.
    fixed = sm.build_partition(sp, tp, 1, "lithology", mesh_path=mesh)
    assert np.isclose(fixed["vertices"][:, 2].max(), -20.0)
    assert fixed["vertices"][:, 2].min() <= -50.0 + 1e-6   # base preserved


def test_stratamesh_missing_field_errors(tmp_path):
    """Requesting lithology/provenance fields absent from the file complains."""
    sm = pytest.importorskip("gospl.analyse.stratamesh")
    h5py = pytest.importorskip("h5py")
    d = tmp_path / "h5"
    d.mkdir()
    n, L = 4, 2
    with h5py.File(str(d / "topology.p0.h5"), "w") as f:
        f["coords"] = np.zeros((n, 3))
        f["cells"] = np.array([[1, 2, 3]])
    with h5py.File(str(d / "stratal.0.p0.h5"), "w") as f:    # single-fraction only
        f["stratZ"] = np.zeros((n, L))
        f["stratH"] = np.full((n, L), 10.0)
        f["phiS"] = np.full((n, L), 0.49)
    with pytest.raises(ValueError, match="stratHf"):
        sm.build_partition(
            str(d / "stratal.0.p0.h5"), str(d / "topology.p0.h5"), 0, "lithology"
        )
    with pytest.raises(ValueError, match="stratP"):
        sm.build_partition(
            str(d / "stratal.0.p0.h5"), str(d / "topology.p0.h5"), 0, "provenance"
        )


def test_stratamesh_cli_end_to_end(tmp_path):
    """
    The CLI creates the output directory, writes a per-partition HDF5 per step,
    and a single well-formed XDMF referencing every partition block with Wedge
    topology; ``--tout`` maps the step index to a simulation time.
    """
    sm = pytest.importorskip("gospl.analyse.stratamesh")
    h5py = pytest.importorskip("h5py")
    h5dir, n, m, L, C, nparts = _write_dataset(tmp_path, nparts=2)
    outdir = str(tmp_path / "vol")                       # created by the tool

    rc = sm.main([
        "--h5dir", h5dir, "--outdir", outdir, "--field", "provenance",
        "--steps", "0", "--tout", "50", "--tstart", "1000",
    ])
    assert rc == 0

    # Default --out prefix is 'strata'; everything lands inside --outdir.
    xdmf = os.path.join(outdir, "strata.xdmf")
    assert os.path.isdir(outdir) and os.path.exists(xdmf)
    ET.parse(xdmf)                                       # raises on malformed XML
    txt = open(xdmf).read()
    assert txt.count('TopologyType="Wedge"') == nparts
    assert "src_class0" in txt and "src_class1" in txt and "dominant" in txt
    # --tout/--tstart: ParaView Time = tstart + step*tout = 1000 + 0*50 = 1000.
    assert 'Time Value="1000"' in txt
    for p in range(nparts):
        ph5 = os.path.join(outdir, "strata.0.p%d.h5" % p)
        assert os.path.exists(ph5)
        with h5py.File(ph5, "r") as f:
            assert f["vertices"].shape == (n * (L - 1), 3)   # lo=1 -> 2 surfaces
            assert f["wedges"].shape == (m * (L - 2), 6)
            assert f["frac"].shape[1] == C
