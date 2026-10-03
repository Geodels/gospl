"""
Flexural isostasy: flat FEM biharmonic and global DH-grid placement.

Protects: AGENTS.md > _extraFlex, Fixed (domain radius hang).

Split out of the former tests/test_regression.py. Run this group alone with
`pytest -m flexure`; see tests/README.md for the marker list.
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

pytestmark = [pytest.mark.flexure]


def test_flex_fem_2d_solver(flat_fem_flex_model):
    """
    Protects: the opt-in parallel FV biharmonic flexure solver for FLAT models
    (`flexure: method: fem`) runs end-to-end on the DMPlex (no regular
    grid) and produces a finite, non-trivial flexural field with subsidence
    under deposition.
    """
    model = flat_fem_flex_model
    assert model.flexOn and model.flex_method == "fem"

    model.runProcesses()
    flx = model.localFlex
    assert np.isfinite(flx).all(), "FEM-2D flexural field non-finite"
    # The incising fixture is net-erosional (unloading → isostatic rebound), so
    # the deflection is positive; only require a finite, non-trivial response.
    assert np.abs(flx).max() > 0.0, "FEM-2D flexure produced no deflection"


def test_flex_fem_2d_physical(flat_fem_flex_model):
    """
    Protects: the FEM-2D solver is physically sensible. A smooth deposition cap
    in the middle of the domain must produce subsidence under the load (negative
    deflection), a finite field, and a magnitude bounded by the local-isostatic
    limit q/Δρg = dz·ρ_s/ρ_a (a stiff plate deflects LESS than full isostasy).
    """
    from mpi4py import MPI

    model = flat_fem_flex_model
    assert model.flex_method == "fem"

    xy = model.lcoords[:, :2]
    cx, cy, amp, wid = 8000.0, 8000.0, 200.0, 1500.0
    load = amp * np.exp(-(((xy[:, 0] - cx) ** 2 + (xy[:, 1] - cy) ** 2)
                          / (2 * wid ** 2)))
    w = model._cmptFlexFEM(load)             # returns -deflection (neg = subsidence)
    assert np.isfinite(w).all(), "FEM-2D field non-finite"

    owned = model.inIDs == 1
    wmin = -MPI.COMM_WORLD.allreduce(
        float(-(w[owned].min())) if owned.sum() else -1e30, op=MPI.MAX
    )
    wabs = MPI.COMM_WORLD.allreduce(
        float(np.abs(w[owned]).max()) if owned.sum() else 0.0, op=MPI.MAX
    )
    assert wmin < 0.0, "no subsidence under a deposition load"
    iso = amp * model.flex_rhos / model.flex_rhoa     # local-isostatic deflection
    assert 0.0 < wabs < 1.5 * iso, (
        f"deflection {wabs:.1f} m outside (0, 1.5x isostatic={1.5 * iso:.1f})"
    )


def test_flex_fem_2d_clamped(flat_fem_flex_model):
    """
    Protects: the 0Displacement0Slope (clamped) boundary of the FEM-2D solver.
    All four sides are clamped and a BROAD load deflects the plate to the edges;
    w must be pinned to ~0 on the clamped edge nodes (relative to the interior
    peak), with a non-trivial interior deflection.
    """
    from mpi4py import MPI

    model = flat_fem_flex_model
    model.flex_bcN = model.flex_bcS = "0Displacement0Slope"
    model.flex_bcE = model.flex_bcW = "0Displacement0Slope"

    xy = model.lcoords[:, :2]
    cx, cy, wid = 8000.0, 8000.0, 4000.0
    load = 300.0 * np.exp(-(((xy[:, 0] - cx) ** 2 + (xy[:, 1] - cy) ** 2)
                            / (2 * wid ** 2)))
    w = model._cmptFlexFEM(load)

    def _max(x):
        return MPI.COMM_WORLD.allreduce(float(x), op=MPI.MAX)

    owned = model.inIDs == 1
    edge = np.unique(np.concatenate([model.southPts, model.northPts,
                                     model.eastPts, model.westPts]))
    edge = edge[model.inIDs[edge] == 1]
    w_edge = _max(float(np.abs(w[edge]).max()) if len(edge) else 0.0)
    w_peak = _max(float(np.abs(w[owned]).max()) if owned.sum() else 0.0)
    assert w_peak > 0.0, "no deflection"
    assert w_edge < 1.0e-3 * w_peak, (
        f"clamped edge not pinned: |w|_edge={w_edge:.3e} vs peak {w_peak:.3e}"
    )


def test_dh_flexure_grid_placed_at_mesh_radius():
    """
    Protects: `addprocess._buildDHGrid` must place the Driscoll-Healy flexure
    grid on the sphere the MESH occupies, not on the declared `domain: radius`.

    Silent failure prevented — and it is not silent so much as fatal: that grid
    exists only to supply DIRECTIONS for the mesh <-> grid interpolation, which
    is a 3-D nearest-neighbour + inverse-distance match against `self.mCoords`
    via `cKDTree`. If the two shells have different radii, every query point
    sits a constant radial offset from the whole mesh, all `k` neighbour
    distances become near-degenerate, and the KD-tree loses all pruning. On a
    real 5.9M-node mesh, declaring a Mars radius (3389.5 km) against an
    Earth-radius mesh (6371.2 km) took the grid query from ~2 s to ~4 h. Because
    `_buildDHGrid` runs on RANK 0 ONLY inside `Model.__init__` while every other
    rank waits at the next collective, that presents as a hang during
    initialisation rather than as a slow step.

    The declared radius is still what sets the flexural wavelength, so it must
    stay in the elastic-operator eigenvalues `dh_P_l ~ 1 / radius**4`.
    """
    addprocess = pytest.importorskip(
        "gospl.tools.addprocess",
        reason="goSPL runtime deps not installed",
    )

    # Bare instance: `_buildDHGrid` only reads mCoords / radius / flex_res_deg /
    # rgrd_interp, so skip the whole GridProcess bootstrap (the `__new__` stub
    # pattern used elsewhere in this file).
    grid = addprocess.GridProcess.__new__(addprocess.GridProcess)
    meshR = 6371220.0
    rng = np.random.default_rng(0)
    xyz = rng.normal(size=(4000, 3))
    xyz *= (meshR / np.linalg.norm(xyz, axis=1))[:, None]
    grid.mCoords = xyz
    grid.flex_res_deg = 10.0          # tiny DH grid: 18 x 36
    grid.rgrd_interp = 4

    # A radius that agrees with the mesh, and one that does not (Mars).
    grid.radius = meshR
    grid._buildDHGrid()
    w_match = grid.dhWeights.copy()
    pl_match = grid.dh_P_l.copy()

    grid.radius = 3389500.0
    grid._buildDHGrid()

    # NOTE the neighbour IDS are not a useful assertion: for query points on a
    # concentric shell the squared distance is monotonic in the angle, so the
    # k-nearest SET and its order survive any radial offset. What breaks is the
    # inverse-distance WEIGHTS (`1/d**2`), which flatten to near-uniform once
    # every neighbour sits ~offset away, degrading the mesh->DH interpolation to
    # a plain k-point average; and, far more seriously, the cKDTree query time,
    # because near-degenerate distances defeat all of the tree's pruning.
    assert np.allclose(grid.dhWeights, w_match, rtol=1.0e-12, atol=0.0), (
        "the mesh->DH inverse-distance weights changed when only the DECLARED "
        "radius changed, so the grid is being placed at `self.radius` instead "
        "of the mesh radius."
    )
    # The interpolation must stay LOCAL: neighbour distances of order the mesh
    # spacing, not of order the radius mismatch.
    dist = 1.0 / np.sqrt(grid.dhWeights[grid.dhWeights > 0.0])
    spacing = np.sqrt(4.0 * np.pi * meshR ** 2 / len(xyz))
    assert dist.max() < 10.0 * spacing, (
        f"largest mesh->DH neighbour distance {dist.max():.4g} m is far beyond "
        f"the mesh spacing {spacing:.4g} m: the DH grid is on a different shell "
        f"from the mesh."
    )
    # atol=0: these eigenvalues are ~1e-26, so the default absolute
    # tolerance would call any two of them equal.
    assert not np.allclose(grid.dh_P_l, pl_match, atol=0.0), (
        "`dh_P_l` did not change with the declared radius. The flexural "
        "eigenvalues are physical and MUST scale as 1/radius**4 — only the "
        "interpolation geometry should follow the mesh."
    )
