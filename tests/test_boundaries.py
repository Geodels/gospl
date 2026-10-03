"""
Flat-model boundary conditions (o/f/w/c) and closed-sink deposition.

Protects: AGENTS.md > Fixed (open-edge runaway, closed-sink needle).

Split out of the former tests/test_regression.py. Run this group alone with
`pytest -m boundaries`; see tests/README.md for the marker list.
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

pytestmark = [pytest.mark.boundaries]


def test_closed_sink_deposition_scope(incising_model, flat_wall_model, minimal_model):
    """
    Protects: `_closedDepo` (sedplex) is gated on `_domainHasOutlet` so a
    flow-terminal stray sink (`lsink & pitID<0`, e.g. a near-flat low pocket the
    flat-routing can't connect to its outlet) is deposited in place ONLY on a
    fully-CLOSED domain (sphere / all-wall) — to conserve mass — and is otherwise
    let LEAVE the domain on an OPEN domain.

    Silent failure prevented: PR #473's closed-sink closure fired on open domains
    too, needling all such sediment onto one near-flat edge node (a ~9 km spike
    on the glacial soil example), which then stalled the SPL/soil SNES. The fix
    keeps conservation on closed domains (sphere / all-wall) while open domains
    drain it out. `_domainHasOutlet` must therefore be False exactly when there
    is no draining outlet.
    """
    # Open flat domain (incising.yml is '0000') -> has draining outlets.
    assert incising_model._domainHasOutlet is True
    assert len(incising_model.outletIDs) > 0
    # All-wall flat box -> no outlet -> closed (conserves via _closedDepo).
    assert flat_wall_model._domainHasOutlet is False
    assert len(flat_wall_model.outletIDs) == 0
    # Global sphere -> no borders at all -> closed.
    assert minimal_model._domainHasOutlet is False


def test_wall_boundary_conservation(flat_wall_model):
    """
    Protects: WALL (closed) boundaries contain flow AND sediment for a SINGLE
    step. With every edge a wall and the sea far below the domain, the flat model
    is a closed box, so the sediment budget must balance over a step: eroded
    material that reaches the closed terminal sinks is deposited in place (the
    `_closedDepo` closure) rather than draining out — total deposited == total
    eroded, surface volume conserved.

    Scope: the fixture runs ONE step. Full conservation over MANY steps in a
    *fully* closed domain is NOT guaranteed — once the basins saturate, goSPL's
    spill-based cascade cannot aggrade a closed basin (it assumes excess always
    spills toward an outlet). That is a documented architectural limitation; the
    common case (walls plus at least one open/fixed edge, or marine) conserves.
    """
    from mpi4py import MPI

    model = flat_wall_model
    assert model.flatModel
    assert (model.south, model.east, model.north, model.west) == (3, 3, 3, 3), \
        "all edges should be wall (flag 3)"
    assert len(model.outletIDs) == 0, "all-wall domain should have no draining outlets"

    model.runProcesses()

    owned = model.inIDs == 1
    cumED = model.cumEDLocal.getArray()[owned]
    la = model.larea[owned]

    def _sum(x):
        return MPI.COMM_WORLD.allreduce(float(x), op=MPI.SUM)

    dV = _sum((cumED * la).sum())                 # net volume change (should ~0)
    activity = _sum((np.abs(cumED) * la).sum())   # volume actually redistributed
    assert activity > 0.0, "no erosion/deposition occurred — test is vacuous"
    rel = abs(dV) / activity
    assert rel < 1.0e-4, (
        f"wall boundary leaks: |dV|/activity = {rel:.2e} (sediment not conserved)"
    )


def test_open_edge_non_rising(incising_model):
    """
    Protects: open ('o') boundary edges are a NON-RISING free-outflow base
    level (`unstructuredmesh._drainOpenEdges`).

    Silent failure prevented: historically an open edge was reset every step to
    the neighbour AVERAGE (`getbc`/`fitedges`), so on an AGGRADING downstream
    edge the edge rose with the interior, lost base-level control, and the
    near-edge plain ponded sediment into a closed-sink deposition runaway
    (tens-of-km spurious peaks on the low, downstream open edge of an
    escarpment-retreat model). The fix sets each open node to
    `min(current elevation, min(interior neighbours))`: it may follow the domain
    DOWN (incision / subsidence) but must never be pushed UP by an aggrading
    interior. This test pins both directions directly (fast + deterministic;
    the full runaway is mesh-resolution dependent and impractical to reproduce
    in a tiny fixture).
    """
    model = incising_model
    assert model.flatModel
    # incising.yml is all-open ('0000'); collect every open-edge node.
    open_nodes = np.unique(np.concatenate([
        pts for flag, pts in (
            (model.north, model.northPts), (model.east, model.eastPts),
            (model.south, model.southPts), (model.west, model.westPts),
        ) if flag == 0 and pts is not None and len(pts) > 0
    ]))
    assert len(open_nodes) > 0, "incising fixture should have open edges"
    is_open = np.zeros(model.lpoints, dtype=bool)
    is_open[open_nodes] = True

    base = model.hLocal.getArray().copy()
    edge0 = base[open_nodes].copy()

    # 1) AGGRADATION: raise the whole interior far above the edge. The edge must
    #    NOT rise — min(current, min(neighbours)) stays at the current value.
    agg = base.copy()
    agg[~is_open] += 1000.0
    out = model._drainOpenEdges(agg.copy())
    assert np.all(out[open_nodes] <= edge0 + 1.0e-9), (
        "open edge rose with an aggrading interior (lost base-level control): "
        f"max rise = {float((out[open_nodes] - edge0).max()):.3e} m"
    )
    assert np.allclose(out[open_nodes], edge0, atol=1.0e-9), (
        "open edge value drifted under pure aggradation"
    )

    # 2) INCISION: drop the whole interior below the edge. The edge MUST follow
    #    down (free outflow tracks the lowering domain).
    inc = base.copy()
    inc[~is_open] -= 500.0
    out2 = model._drainOpenEdges(inc.copy())
    # every open node with an interior neighbour must have dropped
    assert np.all(out2[open_nodes] <= edge0 + 1.0e-9), "incision case raised the edge"
    assert float((edge0 - out2[open_nodes]).max()) > 1.0, (
        "open edge did not track the incising interior down"
    )
    assert np.all(np.isfinite(out2)), "open-edge fill produced non-finite values"


def test_open_edge_marine_tracks(incising_model):
    """
    Protects: the sea-level gate in `_drainOpenEdges`. BELOW sea level an open
    edge is a marine depocenter, not a subaerial outflow, so the min-neighbour
    reset is dropped there — the edge HOLDS its current elevation (like fixed)
    instead of being pulled down to its deepest interior neighbour. The
    pull-down deepens the edge relative to the aggrading basin and manufactures
    accommodation (a spurious sediment wedge against the lateral marine edges).
    The subaerial running-min (test_open_edge_non_rising) is untouched, so the
    escarpment base-level control is preserved.
    """
    model = incising_model
    assert model.flatModel
    open_nodes = np.unique(np.concatenate([
        pts for flag, pts in (
            (model.north, model.northPts), (model.east, model.eastPts),
            (model.south, model.southPts), (model.west, model.westPts),
        ) if flag == 0 and pts is not None and len(pts) > 0
    ]))
    assert len(open_nodes) > 0
    is_open = np.zeros(model.lpoints, dtype=bool)
    is_open[open_nodes] = True

    sl = model.sealevel
    # Whole domain below sea level (marine); the interior is DEEPER than the
    # edge, so the old min-neighbour reset would pull the edge down to sl-400.
    arr = np.full(model.lpoints, sl - 200.0)      # edge at sl-200
    arr[~is_open] = sl - 400.0                     # interior deeper
    out = model._drainOpenEdges(arr.copy())

    # Marine edge HOLDS its elevation (sl-200); it is NOT dragged down to the
    # deeper interior (sl-400), which would manufacture accommodation.
    assert np.allclose(out[open_nodes], sl - 200.0, atol=1.0e-6), (
        "marine open edge was pulled down to its deepest neighbour (would "
        "manufacture accommodation against the basin)"
    )
    assert np.all(out[open_nodes] > sl - 400.0 + 1.0)
    assert np.all(np.isfinite(out))


def test_cyclic_boundary(cyclic_cyl_model):
    """
    Protects: the cyclic (periodic) flow/sediment boundary option (`bc: '0c0c'`
    → E/W cyclic, N/S open) on a cylinder mesh. Checks that:
      - the cyclic edges parse to the cyclic flag (2), open edges to 0;
      - the FV neighbour graph wraps across the seam (the cylinder's seam cells
        link the two periodic edges) — this is what makes flow/sediment route
        across the boundary;
      - the cyclic seam is NOT given the open-outflow sentinel (it is not a
        drain — only the genuinely open N/S edges are);
      - the model runs end-to-end with finite flow accumulation.
    """
    model = cyclic_cyl_model
    assert model.flatModel
    # bc '0c0c' = [S=0 open, E=c cyclic→2, N=0 open, W=c cyclic→2]
    assert (model.east, model.west) == (2, 2), "E/W should be cyclic (2)"
    assert (model.south, model.north) == (0, 0), "N/S should be open (0)"

    # Cross-seam wrap: nodes on the θ≈0 seam (x≈xmax) have FV neighbours on the
    # θ≈2π side (x≈xmax, opposite embedding-z sign) — i.e. the graph is periodic.
    xmax = model.lcoords[:, 0].max()
    seam = np.where(np.isclose(model.lcoords[:, 0], xmax, atol=1.0))[0]
    ng = model.FVmesh_ngbID
    cross = 0
    for s in seam:
        nbrs = ng[s][ng[s] >= 0]
        for n in nbrs:
            if (abs(model.lcoords[n, 0] - xmax) < 200.0
                    and abs(model.lcoords[n, 2] + model.lcoords[s, 2]) < 200.0):
                cross += 1
    assert cross > 0, "no cross-seam FV links — mesh is not periodic / not wrapping"

    model.runProcesses()

    # The cyclic seam must not be drained like an open edge (no deep sentinel).
    h = model.hLocal.getArray()
    assert (h[seam] > -1.0e6).all(), "cyclic seam was forced to the open sentinel"

    fa = model.FAL.getArray()
    assert np.isfinite(fa).all() and fa.max() > 0.0, "flow accumulation invalid"


def test_cyclic_advection(cyclic_advect_model):
    """
    Protects: horizontal advection across a cyclic (periodic) seam. A cyclic 2D
    model runs on a cylinder, so a flat (vx, vy) displacement must be remapped
    onto the cylinder tangent — the periodic component becoming motion AROUND
    the seam (`_cylinderVelocity`), and the seam nodes kept out of the advection
    Dirichlet (`advectBorders`). Checks that:
      - the cyclic boundary is detected and the seam is excluded from the
        advection border set (advectBorders ⊂ idBorders);
      - the velocity transform is exactly tangent to the cylinder (zero radial
        component) and preserves the around-seam speed and the axial component;
      - an elevation bump started just inside the seam is advected ACROSS it
        (its peak wraps from +θ past ±π to the far side) and the far side gains
        mass — i.e. material genuinely crosses the periodic boundary;
      - mass is essentially conserved (the seam does not leak) and stays finite.
    """
    model = cyclic_advect_model
    assert model.flatModel and model.cyclicBC, "cyclic model not detected"
    # E/W periodic (bc 'ococ'): the seam must be dropped from the advection
    # Dirichlet, so advectBorders is a strict subset of idBorders.
    assert model.cyclicPts is not None and len(model.cyclicPts) > 0
    assert len(model.advectBorders) < len(model.idBorders)
    assert not np.isin(model.advectBorders, model.cyclicPts).any()

    # --- velocity transform: flat (vx, vy) -> cylinder tangent ---
    vx = 0.25
    hdisp = np.zeros((model.lpoints, 3))
    hdisp[:, 0] = vx          # around-seam (periodic-x) component
    vel = model._cylinderVelocity(hdisp)
    co = model.lcoords
    rad = np.sqrt(co[:, 0] ** 2 + co[:, 2] ** 2)
    rad[rad == 0.0] = 1.0
    normal = np.c_[co[:, 0] / rad, np.zeros(model.lpoints), co[:, 2] / rad]
    assert np.abs(np.sum(vel * normal, axis=1)).max() < 1.0e-9, \
        "transformed velocity is not tangent to the cylinder"
    assert np.abs(vel[:, 1]).max() < 1.0e-12, "spurious axial velocity"
    np.testing.assert_allclose(np.linalg.norm(vel, axis=1), vx, rtol=1e-9)

    # --- advect the bump across the seam ---
    theta = np.arctan2(model.lcoords[:, 2], model.lcoords[:, 0])
    base = -500.0
    z0 = model.hLocal.getArray().copy()
    far = theta < -1.5                       # far side of the seam (θ ≈ -π)
    peak0 = theta[np.argmax(z0)]
    mass0 = float(np.sum((z0 - base) * model.larea))
    farmass0 = float(np.sum((z0[far] - base) * model.larea[far]))
    assert peak0 > 2.0, "bump should start just inside the +θ seam"

    model.runProcesses()

    z1 = model.hLocal.getArray().copy()
    peak1 = theta[np.argmax(z1)]
    mass1 = float(np.sum((z1 - base) * model.larea))
    farmass1 = float(np.sum((z1[far] - base) * model.larea[far]))

    assert np.isfinite(z1).all(), "advected elevation is not finite"
    # The peak wrapped from the +θ side across the seam to the far (−θ) side.
    assert peak1 < 0.0, "bump peak did not cross the periodic seam"
    # Material accumulated on the far side of the seam.
    assert farmass1 > 3.0 * farmass0, "no mass transported across the seam"
    # The seam does not leak: total mass is conserved (small numerical
    # diffusion of the sharp bump is expected).
    assert abs(mass1 - mass0) / abs(mass0) < 0.1, "mass not conserved across seam"
