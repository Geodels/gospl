"""
scripts/fill_mesh_pits.py (input-mesh pre-filling).

Protects: AGENTS.md > Input-mesh contract.

Split out of the former tests/test_regression.py. Run this group alone with
`pytest -m flow`; see tests/README.md for the marker list.
"""

from __future__ import annotations

import numpy as np
import pytest

# Skip the whole module (rather than erroring at collection) when the
# goSPL runtime stack is not installed.
pytest.importorskip(
    "gospl.tools.inputparser",
    reason="goSPL runtime deps (petsc4py / ruamel.yaml / scipy) not installed",
)

pytestmark = [pytest.mark.flow]


def _load_fill_mesh_pits():
    """Import `scripts/fill_mesh_pits.py` by path (it is a script, not a module)."""
    import importlib.util
    import pathlib

    path = pathlib.Path(__file__).resolve().parents[1] / "scripts" / "fill_mesh_pits.py"
    spec = importlib.util.spec_from_file_location("fill_mesh_pits", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_flat_resolve_is_sink_free_and_reshapes_drainage():
    """
    Protects: `scripts/fill_mesh_pits.py --flat-resolve` (Garbrecht & Martz
    flat resolution) must leave NO cell without a lower neighbour, must never
    lower the input, and must actually change the drainage across a filled
    flat (otherwise it is an expensive no-op).

    Why it exists: Priority-Flood + epsilon drains a filled flat along the
    order the flood reached each cell, so flow radiates from the spill point
    in geometric bands. Reshaping that gradient is easy to get subtly wrong,
    and both ways were hit while writing it:

    1. Capping a flat's rise by its own drop protects the step DOWN to what it
       drains into, but not the case where the flat is the LOWER member of a
       stacked pair and climbs into the cell above it. The cap has to be
       symmetric: half the smallest gap across ANY boundary edge, so two
       adjacent flats each rise by at most half of their mutual gap.
    2. Resolving on a ZERO-increment fill is unrecoverable: a filled basin can
       come out exactly level with the ground it drains through, so no cell has
       a lower neighbour and the way out is simply not in the surface. The
       flood order carries it, so the fill must keep a token increment (the
       fortran backend gets this free from `nearest()`).

    Each mistake silently reintroduced sinks into a tool whose whole purpose
    is removing them, which is why this asserts the sink count directly.
    """
    pytest.importorskip("scipy", reason="flat resolution needs scipy")
    from scipy.spatial import Delaunay

    fmp = _load_fill_mesh_pits()

    # A tilted plane with a wide flat-bottomed basin punched into it.
    n = 40
    gx, gy = np.meshgrid(np.arange(n, dtype=float), np.arange(n, dtype=float))
    pts = np.column_stack([gx.ravel(), gy.ravel()])
    cells = Delaunay(pts).simplices.astype(np.int32)
    elev = 100.0 + 0.5 * pts[:, 0] + 0.2 * pts[:, 1]
    elev[np.hypot(pts[:, 0] - 20, pts[:, 1] - 20) < 9] = 96.0
    npoints = len(elev)

    indptr, indices = fmp.neighbour_csr(npoints, cells)
    outlets = np.zeros(npoints, dtype=bool)
    outlets[fmp.hull_nodes(cells)] = True

    flood = fmp.fill_python(elev, indptr, indices, outlets, 1.0e-4, False)
    token = fmp.fill_python(elev, indptr, indices, outlets, 1.0e-9, False)
    resolved, nflats, ncells = fmp.resolve_flats(
        token, indptr, indices, 1.0e-4, 1.0e-6
    )

    assert nflats > 0 and ncells > 0, "the fixture produced no flat to resolve."
    assert (resolved >= elev).all(), "flat resolution lowered the input surface."
    assert (resolved >= token).all(), "flat resolution lowered the filled surface."

    # The plain epsilon fill is sink-free; resolution must not undo that.
    flood_sinks = fmp.count_sinks(flood, indptr, indices, outlets, False)
    resolved_sinks = fmp.count_sinks(resolved, indptr, indices, outlets, False)
    assert flood_sinks == 0, "fixture is wrong: the plain fill left sinks."
    assert resolved_sinks == 0, (
        f"flat resolution reintroduced {resolved_sinks} sink(s); every flat "
        f"cell must keep a strictly lower neighbour."
    )

    # ...and it must actually re-route the flow. NOTE the assertion is that the
    # directions CHANGE, not that the flow gets more concentrated: on a large
    # flat the away-from-the-rim term does visibly channelise it (peak drainage
    # area +73% on a 610-cell flat), but on a small one the effect washes out,
    # so a concentration threshold would be a flaky test rather than a real
    # guarantee.
    def receivers(surface):
        rcv = np.full(npoints, -1)
        for node in range(npoints):
            nbrs = indices[indptr[node]:indptr[node + 1]]
            low = nbrs[np.argmin(surface[nbrs])]
            if surface[low] < surface[node]:
                rcv[node] = low
        return rcv

    flat = flood > elev
    changed = int((receivers(flood)[flat] != receivers(resolved)[flat]).sum())
    assert changed > 0.1 * flat.sum(), (
        f"only {changed} of {int(flat.sum())} flat cells changed flow "
        f"direction: the resolved gradient is not reaching the surface, so "
        f"the flood-order pattern is still what drains the flat."
    )


def test_fill_mesh_pits_flow_accumulation():
    """
    Protects: `scripts/fill_mesh_pits.py --flow-accum`, the QC view used to
    check a conditioned mesh before it is handed to a simulation.

    Three properties, each of which a broken sweep would violate silently and
    a plot would not reveal: every cell is counted exactly once at exactly one
    terminal (so the totals over terminals equal the mesh size), accumulation
    never decreases downstream, and the Numba and pure-Python kernels agree
    (they are the same source compiled two ways, the `auto` contract the
    analysis tools use, so a divergence means the njit path is wrong).
    """
    pytest.importorskip("scipy", reason="needs scipy for the test mesh")
    from scipy.spatial import Delaunay

    fmp = _load_fill_mesh_pits()

    n = 30
    gx, gy = np.meshgrid(np.arange(n, dtype=float), np.arange(n, dtype=float))
    pts = np.column_stack([gx.ravel(), gy.ravel()])
    cells = Delaunay(pts).simplices.astype(np.int32)
    elev = 100.0 + 0.5 * pts[:, 0] + 0.2 * pts[:, 1]
    elev[np.hypot(pts[:, 0] - 15, pts[:, 1] - 15) < 6] = 97.0
    npoints = len(elev)

    indptr, indices = fmp.neighbour_csr(npoints, cells)
    outlets = np.zeros(npoints, dtype=bool)
    outlets[fmp.hull_nodes(cells)] = True
    filled = fmp.fill_python(elev, indptr, indices, outlets, 1.0e-4, False)

    accum, receiver = fmp.flow_accumulation(
        filled, indptr, indices, outlets, method="python"
    )

    terminal = receiver < 0
    assert accum[terminal].sum() == npoints, (
        f"accumulation over the terminal cells is {accum[terminal].sum():.0f}, "
        f"not the {npoints} cells of the mesh: flow is being lost or "
        f"double-counted."
    )
    assert (accum >= 1.0).all(), "a cell accumulates less than its own area."
    downstream = receiver >= 0
    assert (accum[receiver[downstream]] >= accum[downstream]).all(), (
        "accumulation decreases downstream, so the sweep is not running in "
        "topological order."
    )

    numba = pytest.importorskip("numba", reason="njit path needs numba")
    del numba
    fast, _ = fmp.flow_accumulation(
        filled, indptr, indices, outlets, method="numba"
    )
    assert np.array_equal(fast, accum), (
        "the numba and python accumulation kernels disagree."
    )
