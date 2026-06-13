"""
Tests for the sediment-provenance post-processor (gospl.analyse.provenance):
source-to-sink attribution by fraction routing on synthetic meshes, plus the
adjacency / receiver / downhill-edge helpers. Pure NumPy — no GIS/HDF5 needed.
"""

import numpy as np
import pytest

prov = pytest.importorskip("gospl.analyse.provenance")


def _chain():
    """A 5-node line 0->1->2->3->4 (node 4 the outlet), spacing 1, decreasing
    elevation so each node drains to the next."""
    coords = np.array([[float(i), 0.0] for i in range(5)])
    elev = np.array([4.0, 3.0, 2.0, 1.0, 0.0])
    indptr = np.array([0, 1, 3, 5, 7, 8])
    indices = np.array([1, 0, 2, 1, 3, 2, 4, 3])
    return coords, elev, (indptr, indices)


def test_build_neighbours_and_receivers():
    # Two triangles sharing an edge: vertices 0-1-2 and 1-2-3.
    cells = np.array([[0, 1, 2], [1, 2, 3]])
    indptr, indices = prov.build_neighbours(cells, 4)
    nbr = {i: set(indices[indptr[i]:indptr[i + 1]].tolist()) for i in range(4)}
    assert nbr[0] == {1, 2}
    assert nbr[1] == {0, 2, 3}
    assert nbr[2] == {0, 1, 3}
    assert nbr[3] == {1, 2}

    elev = np.array([3.0, 2.0, 2.5, 0.0])
    rcv = prov.steepest_receivers(elev, indptr, indices)
    assert rcv[0] == 1            # lowest lower neighbour of 0 is 1
    assert rcv[1] == 3            # 1 -> 3 (0.0)
    assert rcv[3] == 3            # outlet (no lower neighbour)


def test_provenance_mixing_distance_conservation():
    coords, elev, nbr = _chain()
    source_class = np.array([0, 1, 0, 0, 0])      # node0 class 0, node1 class 1
    t = prov.ProvenanceTracker(
        coords, source_class, n_classes=2, neighbours=nbr,
        basin_id=np.array([-1, -1, -1, 0, -1]),    # node 3 is the sink basin
        cu_weight=np.array([1.0, 0.0]),            # class 0 is Cu-fertile
    )
    t.step(elev, np.zeros(5))                      # prime (erodep_prev = 0)
    # Erode 2 m at node 0 and 1 m at node 1; deposit 3 m at node 3.
    t.step(elev, np.array([-2.0, -1.0, 0.0, 3.0, 0.0]))

    assert np.isclose(t.dep.sum(), 3.0)            # all eroded volume deposited
    assert np.isclose(t.exported.sum(), 0.0)
    frac = t.pixel_fractions()[3]                  # 2/3 class-0, 1/3 class-1
    assert np.allclose(frac, [2.0 / 3.0, 1.0 / 3.0])
    pct = t.basin_percentages()[0]
    assert np.allclose(pct, [200.0 / 3.0, 100.0 / 3.0])
    assert np.isclose(t.mean_distance()[3], 8.0 / 3.0)   # (2*3 + 1*2)/3
    assert np.isclose(t.cu_fraction()[3], 2.0 / 3.0)


def test_provenance_recycling():
    """Sediment deposited then re-eroded carries its mixed provenance onward."""
    coords, elev, nbr = _chain()
    source_class = np.array([0, 1, 2, 2, 2])
    t = prov.ProvenanceTracker(coords, source_class, n_classes=3, neighbours=nbr)
    t.step(elev, np.zeros(5))
    # Erode node0 (class0) + node1 (class1); deposit the mix at node 3.
    t.step(elev, np.array([-2.0, -1.0, 0.0, 3.0, 0.0]))
    assert np.allclose(t.pixel_fractions()[3], [2.0 / 3.0, 1.0 / 3.0, 0.0])

    # Re-erode node 3's pile (Δerodep[3] = -3) and deposit at node 4.
    t.step(elev, np.array([-2.0, -1.0, 0.0, 0.0, 3.0]))
    assert t.pile[3].sum() < 1.0e-9, "pile not re-entrained"
    # The recycled deposit at node 4 keeps the 2:1 class0:class1 ratio (no
    # class-2 bedrock added, since only the pile eroded).
    assert np.allclose(t.dep[4], [2.0, 1.0]) if t.dep.shape[1] == 2 else \
        np.allclose(t.dep[4], [2.0, 1.0, 0.0])


def _fan():
    """node 0 (source) drains to two lower neighbours 1 and 2 at equal distance
    but different slope; both are sinks (deposit everything)."""
    coords = np.array([[0.0, 0.0], [1.0, 0.0], [-1.0, 0.0], [0.0, -1.0]])
    elev = np.array([2.0, 1.0, 0.0, -1.0])         # slope 0->1 = 1, 0->2 = 2
    indptr = np.array([0, 2, 4, 6, 8])
    indices = np.array([1, 2, 0, 3, 0, 3, 1, 2])
    return coords, elev, (indptr, indices)


def test_mfd_splits_by_slope():
    coords, elev, nbr = _fan()
    src = np.array([0, 0, 0, 0])
    # MFD: node 0's 3 m^3 splits by slope (1 vs 2) -> 1/3 to node1, 2/3 to node2.
    t = prov.ProvenanceTracker(coords, src, n_classes=1, neighbours=nbr, routing="mfd")
    t.step(elev, np.zeros(4))
    t.step(elev, np.array([-3.0, 10.0, 10.0, 0.0]))   # erode 3 at 0; sinks at 1,2
    assert np.isclose(t.dep.sum(), 3.0)               # conserved
    assert np.isclose(t.dep[1, 0], 1.0)               # slope-1 branch
    assert np.isclose(t.dep[2, 0], 2.0)               # slope-2 branch (steeper)

    # SFD: everything goes to the single steepest neighbour (node 2).
    t2 = prov.ProvenanceTracker(coords, src, n_classes=1, neighbours=nbr, routing="sfd")
    t2.step(elev, np.zeros(4))
    t2.step(elev, np.array([-3.0, 10.0, 10.0, 0.0]))
    assert np.isclose(t2.dep[2, 0], 3.0)
    assert np.isclose(t2.dep[1, 0], 0.0)


def test_source_class_validation():
    """source_class must be a full per-vertex array in [0, n_classes)."""
    coords, _, nbr = _chain()
    # -1 (unclassified) is rejected — assign a background class instead.
    with pytest.raises(ValueError):
        prov.ProvenanceTracker(coords, np.array([-1, 0, 0, 0, 0]), 2, neighbours=nbr)
    # out-of-range class is rejected.
    with pytest.raises(ValueError):
        prov.ProvenanceTracker(coords, np.array([0, 0, 0, 0, 5]), 2, neighbours=nbr)
    # wrong length is rejected.
    with pytest.raises(ValueError):
        prov.ProvenanceTracker(coords, np.array([0, 1]), 2, neighbours=nbr)
    # neither cells nor neighbours -> error.
    with pytest.raises(ValueError):
        prov.ProvenanceTracker(coords, np.zeros(5, int), 1)


def test_gospl_output_partition_reassembly(tmp_path):
    """GosplOutput stitches goSPL's per-partition HDF5 onto the global mesh
    ordering (shared ghost nodes agree)."""
    h5py = pytest.importorskip("h5py")
    pytest.importorskip("scipy")

    # Global mesh: 4 vertices in a line.
    gcoords = np.array([[0.0, 0.0], [1.0, 0.0], [2.0, 0.0], [3.0, 0.0]])
    d = tmp_path
    # Partition 0 owns global {0,1,2}; partition 1 owns {2,3} (node 2 shared).
    with h5py.File(d / "topology.p0.h5", "w") as f:
        f["coords"] = gcoords[[0, 1, 2]]
    with h5py.File(d / "topology.p1.h5", "w") as f:
        f["coords"] = gcoords[[2, 3]]
    with h5py.File(d / "gospl.0.p0.h5", "w") as f:
        f["erodep"] = np.array([[10.0], [20.0], [30.0]])
    with h5py.File(d / "gospl.0.p1.h5", "w") as f:
        f["erodep"] = np.array([[30.0], [40.0]])   # node 2 agrees (30)

    reader = prov.GosplOutput(str(d), gcoords)
    ed = reader.field(0, "erodep")
    assert np.allclose(ed, [10.0, 20.0, 30.0, 40.0])   # reassembled in global order
