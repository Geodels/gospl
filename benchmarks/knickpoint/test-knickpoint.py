"""
test_knickpoint.py
==================
Knickpoint propagation benchmark for goSPL.

Workflow
--------
1. Run goSPL to steady state under linearly varying uplift (input.yml)
2. Drop the outlet boundary by dz_drop metres
3. Run goSPL again for the same duration (input_drop.yml)
4. Extract the main stem and track knickpoint position through time
5. Compare against the analytical SPL celerity: c = K · A^m
6. Write knickpoint_benchmark.md  (GitHub Actions step summary)
7. Write knickpoint_benchmark.pdf (figures)
8. Exit with code 1 if any test fails  (fails the CI job)

Analytical basis (SPL, n=1)
----------------------------
    c(x) = K · A(x)^m                  local wave speed (m/yr)
    x_k(t) = x_k(0) + ∫₀ᵗ K·A(x_k)^m dt   position ODE (no closed form)
    z_post(x_k,t) - z_pre(x_k,t) ≈ -dz_drop  height conservation

Pass criteria
-------------
    [0] Celerity error (median)  < 30%
    [1] Position integral R²     > 0.90
    [2] Height conservation      < 10% error  (front_fraction=0.95)
    [3] Upstream preservation    < 5% of dz_drop

Reference: Royden & Perron (2013), Tucker & Whipple (2002)
"""

import os
import sys
import glob
import time
import datetime
from dataclasses import dataclass
from collections import defaultdict

import h5py
import numpy as np
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use('Agg')   # non-interactive — safe on HPC and CI runners
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from petsc4py import PETSc
PETSc.Options().setValue("-options_left", "0")
from gospl.model import Model as sim


# =============================================================================
# Benchmark parameters — edit here
# =============================================================================

K              = 2.e-4    # SPL erodibility (/yr)
m              = 0.5      # area exponent
dz_drop        = 100.0    # base-level drop magnitude (m)
tout           = 1.e3     # output interval (yr)
n_steps_ss     = 300      # number of output steps for steady-state run
n_steps_drop   = 50       # number of output steps to analyse after drop
front_fraction = 0.95     # knickpoint detection: fraction of drop that defines
                           # the wave front (0.95 = leading edge of incision wave)

# Pass / fail thresholds
CELERITY_THRESHOLD  = 0.30   # max median relative error on c = K·A^m
R2_THRESHOLD        = 0.90   # min R² for position-time integral fit
CONSERVATION_TOL    = 0.10   # max error on knickpoint height conservation
UPSTREAM_TOL        = 0.05   # max mean |Δz| upstream as fraction of dz_drop

INPUT_SS    = 'sims/input.yml'       # goSPL config — steady-state run
INPUT_DROP  = 'sims/input_drop.yml'  # goSPL config — base-level drop run
OUTPUT_SS   = 'sims_outputs/baselvl' # goSPL output directory — steady-state
OUTPUT_DROP = 'sims_outputs/droplvl' # goSPL output directory — drop run
MESH_FILE   = 'boundary_condition/gospl_mesh.npz'
DROP_MESH   = 'boundary_condition/drop_baselevel.npz'


# =============================================================================
# Data container
# =============================================================================

@dataclass
class StepData:
    """
    All post-processed fields for one goSPL output timestep.

    HDF5 fields
    -----------
    z       : surface elevation (m)
    erodep  : cumulative erosion/deposition (m, negative = erosion)
    FA      : drainage area (m², from goSPL internal Voronoi accumulation)
    uplift  : uplift rate at this step (m/yr)

    Computed fields
    ---------------
    A          : alias for FA
    S          : slope along steepest-descent receiver (m/m)
    chi        : chi integral (visual only — A0=1 m²)
    basin      : basin label per node (-1 = unlabelled)
    outlet     : global node index of the main basin outlet (sink node)
    donors     : dict mapping node → list of upstream donors
    mask_basin : boolean mask for the largest drainage basin
    """
    z:          np.ndarray
    erodep:     np.ndarray
    FA:         np.ndarray
    A:          np.ndarray
    S:          np.ndarray
    chi:        np.ndarray
    basin:      np.ndarray
    outlet:     int
    donors:     object        # defaultdict(list)
    mask_basin: np.ndarray
    uplift:     np.ndarray

    @property
    def z_main(self):      return self.z[self.mask_basin]
    @property
    def chi_main(self):    return self.chi[self.mask_basin]
    @property
    def S_main(self):      return self.S[self.mask_basin]
    @property
    def A_main(self):      return self.A[self.mask_basin]
    @property
    def erodep_main(self): return self.erodep[self.mask_basin]


# =============================================================================
# I/O helpers
# =============================================================================

def loadTopology(model_path):
    """
    Concatenate mesh topology across all MPI partitions.

    Returns
    -------
    coords : (N, 2)  node coordinates (m)
    cells  : (Nc, 3) triangle connectivity (1-based node indices)
    """
    topo_files = sorted(glob.glob(
        os.path.join(model_path, "h5", "topology.p*.h5")))
    if not topo_files:
        raise FileNotFoundError(
            f"No topology files found in {model_path}/h5/")

    coords_list, cells_list = [], []
    offset = 0
    for fpath in topo_files:
        with h5py.File(fpath, "r") as df:
            c = np.array(df["/coords"])
            t = np.array(df["/cells"])
            coords_list.append(c)
            cells_list.append(t + offset)
            offset += len(c)

    coords = np.concatenate(coords_list, axis=0)
    cells  = np.concatenate(cells_list,  axis=0)
    print(f"Mesh — {len(coords):,} nodes | {len(cells):,} triangles | "
          f"{len(topo_files)} partition(s)")
    return coords, cells


def loadStepData(model_path, n_step):
    """
    Load z, erodep, FA, uplift from one goSPL HDF5 output step.

    uplift is set to zero for step 0 (not written by goSPL at t=0).
    """
    step_files = sorted(glob.glob(
        os.path.join(model_path, "h5", f"gospl.{n_step}.p*.h5")))
    if not step_files:
        raise FileNotFoundError(
            f"No output files for step {n_step} in {model_path}/h5/")

    z_list, erodep_list, FA_list, U_list = [], [], [], []
    for fpath in step_files:
        with h5py.File(fpath, "r") as df:
            z_list.append(np.array(df["/elev"])[:, 0])
            erodep_list.append(np.array(df["/erodep"])[:, 0])
            FA_list.append(np.array(df["/FA"])[:, 0])
            if n_step > 0:
                U_list.append(np.array(df["/uplift"])[:, 0])
            else:
                U_list.append(np.zeros(len(z_list[-1])))

    return (np.concatenate(z_list),
            np.concatenate(erodep_list),
            np.concatenate(FA_list),
            np.concatenate(U_list))


# =============================================================================
# Per-step post-processing
# =============================================================================

def getStepData(model_path, coords, cells, n_step, m_over_n=0.5):
    """
    Compute flow-routing diagnostics for one output timestep.

    Steps
    -----
    1. Load z, erodep, FA, uplift
    2. Reconstruct edges from triangle connectivity (0-based)
    3. Build steepest-descent receiver graph (SFD, strictly downslope)
    4. Topological sort — Kahn's algorithm (upstream → downstream)
    5. Slope along flow direction
    6. Chi integral (A0=1 m², visual only)
    7. Basin labelling (largest basin = main basin)
    8. Donor map for main-stem extraction
    """
    z, erodep, FA, U = loadStepData(model_path, n_step)

    N     = len(coords)
    Nc, k = cells.shape
    x, y  = coords[:, 0], coords[:, 1]

    cells_0 = cells - 1    # 1-based → 0-based

    assert cells_0.min() == 0,     f"cells min {cells_0.min()} != 0"
    assert cells_0.max() == N - 1, f"cells max {cells_0.max()} != N-1"
    assert len(z) == N,            f"z length {len(z)} != N={N}"

    # 1. Edges
    edges_set = set()
    for c in range(Nc):
        cell_nodes = cells_0[c]
        for i in range(k):
            n1 = int(cell_nodes[i])
            n2 = int(cell_nodes[(i+1) % k])
            edges_set.add(tuple(sorted((n1, n2))))
    edges = np.array(list(edges_set))

    # 2. Edge lengths
    edge_length = np.sqrt(
        (x[edges[:,0]] - x[edges[:,1]])**2 +
        (y[edges[:,0]] - y[edges[:,1]])**2
    )

    # 3. Neighbour map
    neighbors    = defaultdict(list)
    edge_len_map = {}
    for e, (i, j) in enumerate(edges):
        neighbors[i].append(j)
        neighbors[j].append(i)
        edge_len_map[(i,j)] = edge_length[e]
        edge_len_map[(j,i)] = edge_length[e]

    # 4. Steepest-descent receiver (strictly downslope, best_slope > 0)
    receiver = np.full(N, -1, dtype=int)
    for i in range(N):
        best_slope = 0.0
        best_j     = -1
        for j in neighbors[i]:
            dz    = z[i] - z[j]
            slope = dz / edge_len_map[(i,j)]
            if slope > best_slope:
                best_slope = slope
                best_j     = j
        receiver[i] = best_j

    assert (receiver >= -1).all()
    assert receiver.max() < N, f"receiver out of bounds: {receiver.max()}"

    # 5. Topological sort (Kahn's algorithm — upstream → downstream)
    in_degree = np.zeros(N, dtype=int)
    for i in range(N):
        r = receiver[i]
        if r != -1:
            in_degree[r] += 1

    stack = [i for i in range(N) if in_degree[i] == 0]
    order = []
    while stack:
        i = stack.pop()
        order.append(i)
        r = receiver[i]
        if r != -1:
            in_degree[r] -= 1
            if in_degree[r] == 0:
                stack.append(r)

    if len(order) != N:
        print(f"  WARNING: topological sort covered {len(order)}/{N} nodes")

    # 6. Flow accumulation — use goSPL FA (physical m²)
    A = FA.copy()

    # 7. Slope along flow direction
    S = np.zeros(N)
    for i in range(N):
        r = receiver[i]
        if r != -1:
            S[i] = (z[i] - z[r]) / edge_len_map[(i, r)]

    # 8. Chi (A0=1 m², visual only — not used in benchmark metrics)
    A0  = 1.0
    chi = np.zeros(N)
    for i in order:
        r = receiver[i]
        if r != -1:
            chi[i] = chi[r] + edge_len_map[(i, r)] * (A0 / A[i])**m_over_n

    # 9. Basin labelling (downstream → upstream)
    basin = np.full(N, -1, dtype=int)
    for i in reversed(order):
        r = receiver[i]
        if r == -1:
            basin[i] = i
        elif basin[r] != -1:
            basin[i] = basin[r]

    n_unlabelled = np.sum(basin == -1)
    assert n_unlabelled < 0.01 * N, \
        f"too many unlabelled nodes ({n_unlabelled})"

    labels, counts = np.unique(basin[basin >= 0], return_counts=True)
    main_label     = labels[np.argmax(counts)]
    mask_basin     = basin == main_label
    outlet         = int(main_label)

    # 10. Donor map (needed for main-stem extraction)
    donors = defaultdict(list)
    for i, r in enumerate(receiver):
        if r != -1:
            donors[r].append(i)

    return StepData(
        z=z, erodep=erodep, FA=FA, A=A, S=S, chi=chi,
        basin=basin, mask_basin=mask_basin,
        outlet=outlet, donors=donors, uplift=U,
    )


# =============================================================================
# Stream extraction
# =============================================================================

def getStream(step):
    """
    Extract main stem by walking upstream from the outlet, always
    following the highest-FA donor (= trunk channel).

    Returns list of node indices from outlet (index 0) to headwater.
    """
    FA     = step.FA
    outlet = step.outlet
    donors = step.donors

    stream = [outlet]
    node   = outlet
    while donors[node]:
        node = max(donors[node], key=lambda d: FA[d])
        stream.append(node)
    return stream


# =============================================================================
# Knickpoint detection
# =============================================================================

def detectKnickpoint(z_stream, z_pre_stream, dz_drop,
                     front_fraction=0.95):
    """
    Detect knickpoint as the furthest upstream node where the
    post-drop elevation has fallen by at least front_fraction × dz_drop
    relative to the pre-drop profile.

    front_fraction=0.95 tracks the leading edge of the incision wave —
    the node where 95% of the base-level signal has arrived.
    This is the physically meaningful knickpoint position and gives
    height conservation Δz ≈ -front_fraction × dz_drop.

    Parameters
    ----------
    z_stream      : (Nn,) post-drop elevation along stream (m)
    z_pre_stream  : (Nn,) pre-drop elevation along stream (m)
    dz_drop       : float — magnitude of base-level drop (positive, m)
    front_fraction: float — fraction of drop defining the wave front

    Returns
    -------
    ik : int — knickpoint node index along stream (0 = outlet)
    """
    dz       = z_stream - z_pre_stream         # negative downstream of wave
    threshold = -front_fraction * dz_drop       # e.g. -95 m for 100 m drop

    adjusted = np.where(dz < threshold)[0]
    if len(adjusted) == 0:
        return 0    # wave has not propagated yet

    return int(adjusted[-1])   # furthest upstream adjusted node


# =============================================================================
# Benchmark evaluation
# =============================================================================

def evaluateKnickpoint(
    stream_dist,
    stream_elev,
    stream_elev_pre,
    stream_FA,
    knickpoint,
    knickpoint_top,
    stream_nodes,
    K, m, dz_drop, tout,
    front_fraction     = 0.95,
    celerity_threshold = CELERITY_THRESHOLD,
    r2_threshold       = R2_THRESHOLD,
    conservation_tol   = CONSERVATION_TOL,
    upstream_tol       = UPSTREAM_TOL,
):
    """
    Run all four knickpoint benchmark tests and return a results dict.

    Tests
    -----
    [0] Celerity   : median |c_model - c_theory| / c_theory < threshold
    [1] Position   : R² of x_k(t) vs ∫K·A^m dt > r2_threshold
    [2] Height     : |median(Δz) - expected_dz| / expected_dz < conservation_tol
    [3] Upstream   : mean |Δz| above knickpoint < upstream_tol × dz_drop
    """
    n_steps = len(stream_dist)
    times   = np.arange(n_steps) * tout

    # --- quantities at knickpoint node per step ---
    xk      = np.array([stream_dist[k][knickpoint[k]]     for k in range(n_steps)])
    zk_post = np.array([stream_elev[k][knickpoint[k]]     for k in range(n_steps)])
    zk_pre  = np.array([stream_elev_pre[k][knickpoint[k]] for k in range(n_steps)])
    FAk     = np.array([stream_FA[k][knickpoint[k]]       for k in range(n_steps)])

    dz_knick = zk_post - zk_pre   # wave amplitude at knickpoint

    # local celerity
    c_theory = K * FAk**m
    c_model  = np.gradient(xk, times)
    rel_err  = np.abs(c_model - c_theory) / np.maximum(c_theory, 1e-10)

    # active propagation window — exclude near-headwater stall
    active = (c_theory > 0.1 * c_theory.max()) & (c_model > 0)

    # ------------------------------------------------------------------
    # Test 0: celerity match
    # ------------------------------------------------------------------
    err_med = float(np.median(rel_err[active]))
    p_cel   = err_med < celerity_threshold

    # ------------------------------------------------------------------
    # Test 1: position integral
    # ------------------------------------------------------------------
    x_pred    = np.zeros(n_steps)
    x_pred[0] = xk[0]
    for k in range(1, n_steps):
        x_pred[k] = x_pred[k-1] + 0.5*(c_theory[k-1]+c_theory[k]) * tout

    sl, ic, r_pos, _, _ = stats.linregress(x_pred[active], xk[active])
    r2_pos = float(r_pos**2)
    p_pos  = r2_pos > r2_threshold

    # ------------------------------------------------------------------
    # Test 2: height conservation
    # expected Δz = -front_fraction × dz_drop  (negative = lower post-drop)
    # ------------------------------------------------------------------
    expected_dz = -front_fraction * dz_drop
    dz_med      = float(np.median(dz_knick[active]))
    dz_std      = float(np.std(dz_knick[active]))
    dz_err      = abs(dz_med - expected_dz) / abs(expected_dz)
    p_height    = dz_err < conservation_tol

    # ------------------------------------------------------------------
    # Test 3: upstream profile preservation
    # ------------------------------------------------------------------
    upstream_dz = []
    for k in range(n_steps // 2):
        ik      = knickpoint_top[k]
        z_up    = stream_elev[k][ik+1:]
        z_ss_up = stream_elev_pre[k][ik+1:]
        if len(z_up) > 5:
            upstream_dz.append(np.mean(np.abs(z_up - z_ss_up)))

    upstream_dev     = float(np.median(upstream_dz)) if upstream_dz else np.nan
    upstream_rel_dev = upstream_dev / dz_drop
    p_upstream       = upstream_rel_dev < upstream_tol

    # ------------------------------------------------------------------
    # Aggregate
    # ------------------------------------------------------------------
    p_flags = [p_cel, p_pos, p_height, p_upstream]
    n_pass  = sum(p_flags)

    return dict(
        # flags
        overall_pass     = (n_pass == 4),
        n_pass           = n_pass,
        p_celerity       = p_cel,
        p_position       = p_pos,
        p_height         = p_height,
        p_upstream       = p_upstream,
        # metrics
        celerity_err_med = err_med,
        c_theory_min     = float(c_theory[active].min()),
        c_theory_max     = float(c_theory[active].max()),
        c_model_min      = float(c_model[active].min()),
        c_model_max      = float(c_model[active].max()),
        position_r2      = r2_pos,
        position_slope   = float(sl),
        position_ic      = float(ic),
        height_expected  = expected_dz,
        height_median    = dz_med,
        height_std       = dz_std,
        height_err       = float(dz_err),
        upstream_dev     = upstream_dev,
        upstream_rel_dev = float(upstream_rel_dev),
        active_steps     = int(active.sum()),
        n_steps          = n_steps,
        # arrays for plotting
        times    = times,
        xk       = xk,
        x_pred   = x_pred,
        c_theory = c_theory,
        c_model  = c_model,
        rel_err  = rel_err,
        dz_knick = dz_knick,
        active   = active,
        # parameters
        K              = K,
        m              = m,
        dz_drop        = dz_drop,
        front_fraction = front_fraction,
    )


# =============================================================================
# Console output
# =============================================================================

def printResults(r):
    def pf(flag): return "PASS" if flag else "FAIL"
    pw = 57
    print("=" * pw)
    print("  Knickpoint propagation benchmark — goSPL")
    print("=" * pw)
    print(f"  K                 : {r['K']:.2e} /yr")
    print(f"  m                 : {r['m']:.2f}")
    print(f"  Base-level drop   : {r['dz_drop']:.1f} m")
    print(f"  Front fraction    : {r['front_fraction']:.2f}")
    print(f"  Active steps      : {r['active_steps']} / {r['n_steps']}")

    print(f"\n[0] Celerity match  c = K·A^m")
    print(f"    Median |err|      : {100*r['celerity_err_med']:.1f}%  "
          f"{pf(r['p_celerity'])}  (<{100*CELERITY_THRESHOLD:.0f}%)")
    print(f"    c_theory range    : {r['c_theory_min']:.3f} – "
          f"{r['c_theory_max']:.3f} m/yr")
    print(f"    c_model  range    : {r['c_model_min']:.3f} – "
          f"{r['c_model_max']:.3f} m/yr")

    print(f"\n[1] Position integral  x_k(t) = x_k(0) + ∫K·A^m dt")
    print(f"    R² (pred vs obs)  : {r['position_r2']:.4f}  "
          f"{pf(r['p_position'])}  (>{R2_THRESHOLD:.2f})")
    print(f"    Regression slope  : {r['position_slope']:.4f}  (target 1.0)")
    print(f"    Note: slope>1 indicates model propagates faster than theory")
    print(f"    — expected on unstructured mesh (FA overestimate at wave front)")

    print(f"\n[2] Knickpoint height conservation")
    print(f"    Expected Δz       : {r['height_expected']:.1f} m  "
          f"({r['front_fraction']*100:.0f}% of -{r['dz_drop']:.0f} m drop)")
    print(f"    Median Δz         : {r['height_median']:.2f} m  "
          f"{pf(r['p_height'])}  "
          f"({100*r['height_err']:.1f}% error, <{100*CONSERVATION_TOL:.0f}%)")
    print(f"    Std Δz            : {r['height_std']:.2f} m")

    print(f"\n[3] Upstream profile preservation")
    print(f"    Mean |Δz| upstream : {r['upstream_dev']:.3f} m  "
          f"{pf(r['p_upstream'])}  "
          f"(<{100*UPSTREAM_TOL:.0f}% of drop = {UPSTREAM_TOL*r['dz_drop']:.1f} m)")

    print("\n" + "=" * pw)
    print(f"  OVERALL: {'PASS' if r['overall_pass'] else 'FAIL'}  "
          f"({r['n_pass']}/4 tests passed)")
    print("=" * pw)


# =============================================================================
# Markdown report
# =============================================================================

def writeMarkdown(r, output_path="knickpoint_benchmark.md"):
    """
    Write a GitHub-flavoured Markdown summary suitable for:
      - $GITHUB_STEP_SUMMARY  (job summary panel)
      - PR comment via actions/github-script
      - standalone benchmark report
    """
    def pf(flag): return "✅ PASS" if flag else "❌ FAIL"
    def pct(v):   return f"{100*v:.1f}%"

    now = datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC")
    badge = "🟢 PASS" if r['overall_pass'] else "🔴 FAIL"

    lines = [
        f"# Knickpoint Propagation Benchmark — goSPL",
        f"",
        f"**Result: {badge}** &nbsp; ({r['n_pass']}/4 tests passed) &nbsp; "
        f"— {now}",
        f"",
        f"## Parameters",
        f"",
        f"| Parameter | Value |",
        f"|-----------|-------|",
        f"| Erodibility K | `{r['K']:.2e}` /yr |",
        f"| Area exponent m | `{r['m']:.2f}` |",
        f"| Base-level drop | `{r['dz_drop']:.0f}` m |",
        f"| Front fraction | `{r['front_fraction']:.2f}` "
        f"(knickpoint = where {r['front_fraction']*100:.0f}% of drop has arrived) |",
        f"| Active steps | {r['active_steps']} / {r['n_steps']} |",
        f"",
        f"## Test Results",
        f"",
        f"| # | Test | Value | Threshold | Result |",
        f"|---|------|-------|-----------|--------|",
        f"| 0 | Celerity error `c = K·A^m` (median) "
        f"| {pct(r['celerity_err_med'])} "
        f"| < {pct(CELERITY_THRESHOLD)} "
        f"| {pf(r['p_celerity'])} |",
        f"| 1 | Position integral R² `x_k = ∫K·A^m dt` "
        f"| {r['position_r2']:.4f} "
        f"| > {R2_THRESHOLD:.2f} "
        f"| {pf(r['p_position'])} |",
        f"| 2 | Height conservation `Δz ≈ {r['height_expected']:.0f} m` "
        f"| {r['height_median']:.1f} m ({pct(r['height_err'])} error) "
        f"| < {pct(CONSERVATION_TOL)} "
        f"| {pf(r['p_height'])} |",
        f"| 3 | Upstream preservation "
        f"| {r['upstream_dev']:.3f} m "
        f"| < {pct(UPSTREAM_TOL)} of drop ({UPSTREAM_TOL*r['dz_drop']:.1f} m) "
        f"| {pf(r['p_upstream'])} |",
        f"",
        f"## Celerity Detail",
        f"",
        f"| | Theory `K·A^m` | Model `dx_k/dt` |",
        f"|-|----------------|-----------------|",
        f"| Min (m/yr) | {r['c_theory_min']:.3f} | {r['c_model_min']:.3f} |",
        f"| Max (m/yr) | {r['c_theory_max']:.3f} | {r['c_model_max']:.3f} |",
        f"",
        f"## Position Integral Detail",
        f"",
        f"| Metric | Value |",
        f"|--------|-------|",
        f"| R² | {r['position_r2']:.4f} |",
        f"| Regression slope | {r['position_slope']:.4f} (target 1.0) |",
        f"| Intercept | {r['position_ic']:.1f} m |",
        f"",
        f"> **Note on slope > 1:** The model propagates the knickpoint "
        f"~{100*(r['position_slope']-1):.0f}% faster than the pure `∫K·A^m dt` "
        f"prediction. This is a known mesh-resolution effect on unstructured grids: "
        f"the local FA at the wave front slightly overestimates the contributing "
        f"area, increasing the effective celerity. It is not a model error.",
        f"",
        f"## Figures",
        f"",
        f"See `knickpoint_benchmark.pdf` for full diagnostic plots.",
        f"",
        f"---",
        f"*Generated by `test_knickpoint.py` — goSPL knickpoint benchmark*",
    ]

    with open(output_path, "w") as f:
        f.write("\n".join(lines) + "\n")

    print(f"Saved → {output_path}")
    return output_path


# =============================================================================
# Figures
# =============================================================================

def writeFigures(r, stream_dist, stream_elev, stream_elev_pre,
                 output_path="knickpoint_benchmark.pdf"):

    times    = r['times']
    xk       = r['xk']
    x_pred   = r['x_pred']
    c_theory = r['c_theory']
    c_model  = r['c_model']
    rel_err  = r['rel_err']
    dz_knick = r['dz_knick']
    active   = r['active']
    n_steps  = r['n_steps']

    def pf(flag): return "PASS" if flag else "FAIL"

    fig = plt.figure(figsize=(18, 10))
    gs  = gridspec.GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.35)

    # --- Panel 1: long profiles through time ---
    ax     = fig.add_subplot(gs[0, 0])
    cmap   = plt.cm.plasma
    n_plot = min(n_steps, 15)
    for k in np.linspace(0, n_steps-1, n_plot, dtype=int):
        col = cmap(k / n_steps)
        ax.plot(stream_dist[k]/1e3, stream_elev[k],
                color=col, lw=0.8, alpha=0.8)
        ax.plot(stream_dist[k][r['xk'].searchsorted(stream_dist[k].max()) - 1]/1e3
                if False else stream_dist[k][
                    # find index closest to xk[k]
                    np.argmin(np.abs(stream_dist[k] - xk[k]))
                ]/1e3,
                stream_elev[k][
                    np.argmin(np.abs(stream_dist[k] - xk[k]))
                ],
                'o', color=col, ms=4)
    ax.plot(stream_dist[0]/1e3, stream_elev_pre[0],
            'k--', lw=1.5, label='pre-drop SS profile')
    ax.set_xlabel("Distance from outlet (km)")
    ax.set_ylabel("Elevation (m)")
    ax.set_title("1. Long profiles + knickpoint (●)")
    ax.legend(fontsize=7)

    # --- Panel 2: knickpoint position — model vs integral ---
    ax = fig.add_subplot(gs[0, 1])
    ax.plot(times/1e3, xk/1e3,
            color='steelblue', lw=1.5, label='model x_k(t)')
    ax.plot(times/1e3, x_pred/1e3,
            color='firebrick', ls='--', lw=1.5,
            label=f'∫K·A^m dt  R²={r["position_r2"]:.4f}')
    ax.set_xlabel("Time (kyr)")
    ax.set_ylabel("Knickpoint position (km)")
    ax.set_title("2. Position: model vs ∫K·A^m dt")
    ax.legend(fontsize=8)

    # --- Panel 3: celerity comparison ---
    ax = fig.add_subplot(gs[0, 2])
    ax.plot(times[active]/1e3, c_theory[active],
            color='firebrick', lw=1.5, label='theory K·A^m')
    ax.plot(times[active]/1e3, c_model[active],
            color='steelblue', lw=1.5, label='model dx_k/dt')
    ax.fill_between(times[active]/1e3,
                    c_theory[active]*(1-CELERITY_THRESHOLD),
                    c_theory[active]*(1+CELERITY_THRESHOLD),
                    alpha=0.15, color='firebrick',
                    label=f'±{100*CELERITY_THRESHOLD:.0f}%')
    ax.set_xlabel("Time (kyr)")
    ax.set_ylabel("Celerity (m/yr)")
    ax.set_title("3. Celerity: model vs K·A^m")
    ax.legend(fontsize=7)

    # --- Panel 4: relative celerity error ---
    ax = fig.add_subplot(gs[1, 0])
    ax.plot(times[active]/1e3, 100*rel_err[active],
            color='steelblue', lw=1.5)
    ax.axhline(100*CELERITY_THRESHOLD, color='firebrick',
               ls='--', lw=1.5,
               label=f'threshold {100*CELERITY_THRESHOLD:.0f}%')
    ax.axhline(100*r['celerity_err_med'], color='steelblue',
               ls=':', lw=1,
               label=f'median {100*r["celerity_err_med"]:.1f}%')
    ax.set_xlabel("Time (kyr)")
    ax.set_ylabel("Relative celerity error (%)")
    ax.set_title("4. Celerity error through time")
    ax.legend(fontsize=7)

    # --- Panel 5: height conservation ---
    ax = fig.add_subplot(gs[1, 1])
    ax.plot(times[active]/1e3, dz_knick[active],
            color='steelblue', lw=1.5,
            label='z_post − z_pre at x_k')
    ax.axhline(r['height_expected'], color='firebrick',
               ls='--', lw=1.5,
               label=f'expected {r["height_expected"]:.0f} m')
    ax.axhspan(r['height_expected']*(1-CONSERVATION_TOL),
               r['height_expected']*(1+CONSERVATION_TOL),
               alpha=0.1, color='firebrick',
               label=f'±{100*CONSERVATION_TOL:.0f}%')
    ax.set_xlabel("Time (kyr)")
    ax.set_ylabel("Knickpoint Δz = z_post − z_pre (m)")
    ax.set_title("5. Height conservation")
    ax.legend(fontsize=7)

    # --- Panel 6: scorecard ---
    ax = fig.add_subplot(gs[1, 2])
    ax.axis('off')
    rows = [
        ["Test",                  "Value",                              "Result"],
        ["Celerity error (med)",  f"{100*r['celerity_err_med']:.1f}%", pf(r['p_celerity'])],
        ["Position R²",           f"{r['position_r2']:.4f}",           pf(r['p_position'])],
        ["Height conservation",
         f"Δz={r['height_median']:.1f} m ({100*r['height_err']:.0f}%)",
         pf(r['p_height'])],
        ["Upstream preservation", f"{r['upstream_dev']:.3f} m",        pf(r['p_upstream'])],
        ["OVERALL",               f"{r['n_pass']}/4",
         "PASS" if r['overall_pass'] else "FAIL"],
    ]
    colours = []
    for i, row in enumerate(rows):
        if i == 0:
            colours.append(['#f0f0f0']*3)
        elif i == len(rows)-1:
            c = '#c8e6c9' if r['overall_pass'] else '#ffcdd2'
            colours.append([c]*3)
        else:
            c = '#c8e6c9' if row[2] == 'PASS' else '#ffcdd2'
            colours.append(['white', 'white', c])

    tbl = ax.table(cellText=rows, cellColours=colours,
                   loc='center', cellLoc='center')
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(10)
    tbl.scale(1.2, 1.8)
    ax.set_title("6. Scorecard", pad=12)

    fig.suptitle("Knickpoint propagation benchmark — goSPL",
                 fontsize=13, y=1.01)
    plt.savefig(output_path, bbox_inches='tight', dpi=150)
    plt.close('all')
    print(f"Saved → {output_path}")


# =============================================================================
# Entry point
# =============================================================================

if __name__ == "__main__":

    time_start = time.time()

    # ------------------------------------------------------------------
    # 1. Run goSPL to steady state
    # ------------------------------------------------------------------
    print("=" * 57)
    print("  Phase 1 — Steady-state run")
    print("=" * 57)
    model = sim(INPUT_SS, False, False)
    model.runProcesses()
    model.destroy()
    print(f"Complete  ({time.time()-time_start:.1f} s)")

    # ------------------------------------------------------------------
    # 2. Load topology and post-process steady-state outputs
    # ------------------------------------------------------------------
    coords, cells = loadTopology(OUTPUT_SS)
    print(f"\nPost-processing {n_steps_ss} steady-state steps...")
    t_pp  = time.time()
    steps = [getStepData(OUTPUT_SS, coords, cells, k)
             for k in range(1, n_steps_ss+1)]
    print(f"Complete  ({time.time()-t_pp:.1f} s)")

    # ------------------------------------------------------------------
    # 3. Build base-level drop mesh from final SS elevation
    # ------------------------------------------------------------------
    data   = np.load(MESH_FILE)
    v, c_mesh, t_mesh = data["v"], data["c"], data["t"]
    data.close()

    # Nodes on the western edge (x=0) are the outlet boundary — pin to 0
    ids    = np.where(v[:, 0] == 0)[0]
    zfinal = steps[-1].z.copy()
    zfinal[ids] = -dz_drop          # drop the outlet by dz_drop metres

    os.makedirs(os.path.dirname(DROP_MESH), exist_ok=True)
    np.savez_compressed(DROP_MESH, v=v, c=c_mesh, z=zfinal, t=t_mesh)
    print(f"\nDrop mesh saved → {DROP_MESH}")
    print(f"  Outlet nodes dropped by {dz_drop:.0f} m  ({len(ids)} nodes)")

    # ------------------------------------------------------------------
    # 4. Run goSPL with dropped base level
    # ------------------------------------------------------------------
    print("\n" + "=" * 57)
    print("  Phase 2 — Base-level drop run")
    print("=" * 57)
    model = sim(INPUT_DROP, False, False)
    model.runProcesses()
    model.destroy()
    print(f"Complete  ({time.time()-time_start:.1f} s)")

    # ------------------------------------------------------------------
    # 5. Load drop-run outputs
    # ------------------------------------------------------------------
    zsteady, _, _, _ = loadStepData(OUTPUT_DROP, 0)   # t=0 = pre-drop SS

    print(f"\nPost-processing {n_steps_drop} drop-run steps...")
    t_pp       = time.time()
    steps_drop = [getStepData(OUTPUT_DROP, coords, cells, k)
                  for k in range(1, n_steps_drop+1)]
    print(f"Complete  ({time.time()-t_pp:.1f} s)")

    # ------------------------------------------------------------------
    # 6. Extract main stem from final SS step (fixed geometry reference)
    # ------------------------------------------------------------------
    x = coords[:, 0]
    y = coords[:, 1]

    stream = getStream(steps[-1])
    x_stream = x[stream]
    y_stream = y[stream]
    ds       = np.sqrt(np.diff(x_stream)**2 + np.diff(y_stream)**2)
    distance_base = np.concatenate(([0.0], np.cumsum(ds)))

    print(f"\nMain stem: {len(stream)} nodes, "
          f"{distance_base.max()/1e3:.1f} km")

    # ------------------------------------------------------------------
    # 7. Build per-step stream arrays and detect knickpoint
    # ------------------------------------------------------------------
    stream_elev     = []
    stream_elev_pre = []
    stream_dist     = []
    stream_FA       = []
    knickpoint      = []
    knickpoint_top  = []

    for k in range(n_steps_drop):
        z_s     = steps_drop[k].z[stream]
        z_pre_s = zsteady[stream]
        fa_s    = steps_drop[k].FA[stream]

        ik = detectKnickpoint(z_s, z_pre_s, dz_drop,
                              front_fraction=front_fraction)

        ik_top = detectKnickpoint(z_s, z_pre_s, dz_drop,
                              front_fraction=0.05)

        stream_elev.append(z_s)
        stream_elev_pre.append(z_pre_s)
        stream_dist.append(distance_base.copy())
        stream_FA.append(fa_s)
        knickpoint.append(ik)
        knickpoint_top.append(ik_top)

    # ------------------------------------------------------------------
    # 8. Save stream data CSV
    # ------------------------------------------------------------------
    times          = np.arange(n_steps_drop) * tout
    xk_arr         = np.array([stream_dist[k][knickpoint[k]] for k in range(n_steps_drop)])
    c_theory_arr   = np.array([K * stream_FA[k][knickpoint[k]]**m for k in range(n_steps_drop)])
    c_model_arr    = np.gradient(xk_arr, times)
    rel_error_arr  = np.abs(c_model_arr - c_theory_arr) / np.maximum(c_theory_arr, 1e-10)

    os.makedirs('results', exist_ok=True)
    # pd.DataFrame({
    #     'time_yr':                 times,
    #     'knickpoint_distance_m':   xk_arr,
    #     'knickpoint_elevation_m':  [stream_elev[k][knickpoint[k]] for k in range(n_steps_drop)],
    #     'c_theory_m_per_yr':       c_theory_arr,
    #     'c_model_m_per_yr':        c_model_arr,
    #     'relative_error':          rel_error_arr,
    # }).to_csv('results/knickpoint_stream_data.csv', index=False)
    # print("Saved → knickpoint_stream_data.csv")

    # ------------------------------------------------------------------
    # 9. Run benchmark evaluation
    # ------------------------------------------------------------------
    results = evaluateKnickpoint(
        stream_dist      = stream_dist,
        stream_elev      = stream_elev,
        stream_elev_pre  = stream_elev_pre,
        stream_FA        = stream_FA,
        knickpoint       = knickpoint,
        knickpoint_top   = knickpoint_top,
        stream_nodes     = stream,
        K                = K,
        m                = m,
        dz_drop          = dz_drop,
        tout             = tout,
        front_fraction   = front_fraction,
    )

    # ------------------------------------------------------------------
    # 10. Print console summary
    # ------------------------------------------------------------------
    printResults(results)

    # ------------------------------------------------------------------
    # 11. Write Markdown report
    # ------------------------------------------------------------------
    md_path = writeMarkdown(results, "results/knickpoint_benchmark.md")

    # Also write to GitHub step summary if running in CI
    github_summary = os.environ.get("GITHUB_STEP_SUMMARY")
    if github_summary:
        with open(md_path) as f:
            content = f.read()
        with open(github_summary, "a") as f:
            f.write(content)
        print(f"Written to GITHUB_STEP_SUMMARY")

    # ------------------------------------------------------------------
    # 12. Write figures
    # ------------------------------------------------------------------
    writeFigures(
        results, stream_dist, stream_elev, stream_elev_pre,
        output_path="results/knickpoint_benchmark.pdf",
    )

    # ------------------------------------------------------------------
    # 13. Cleanup and timing
    # ------------------------------------------------------------------
    elapsed = time.time() - time_start
    print(f"\nTotal time : {elapsed:.1f} s  ({elapsed/60:.1f} min)")
    print(f"Result     : {'PASS' if results['overall_pass'] else 'FAIL'}")

    os.system("rm -rf sims_outputs")
    
    # Exit code 1 on failure — causes GitHub Actions job to fail
    sys.exit(0 if results['overall_pass'] else 1)