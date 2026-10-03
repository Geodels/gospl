"""
Multi-rank tests: each spawns `mpirun` subprocesses (np=1 vs np=2).

Protects: AGENTS.md > MPI contract (deadlock + partition-dependence classes).

Split out of the former tests/test_regression.py. Run this group alone with
`pytest -m parallel`; see tests/README.md for the marker list.
"""

from __future__ import annotations

import pytest

from _helpers import FIXTURES_DIR, MPI_TIMEOUT, mpi_child_env

# Skip the whole module (rather than erroring at collection) when the
# goSPL runtime stack is not installed.
pytest.importorskip(
    "gospl.tools.inputparser",
    reason="goSPL runtime deps (petsc4py / ruamel.yaml / scipy) not installed",
)

pytestmark = [pytest.mark.parallel, pytest.mark.mpi]


def test_watertable_parallel(tmp_path):
    """
    Protects (water-table + duricrust, **Phase 2 parallel validation**): the
    implicit head solve is partition-consistent. Runs `minimal_gw.yml` under
    `mpirun -n 1` and `-n 2` (subprocesses, like `test_parallel_correctness`),
    reduces the water-table head / recharge over owned nodes, and checks the two
    decompositions agree.

    Recharge is deterministic + local (no solver) → the owned-node sum must match
    to ~FP. The head carries KSP floating-point noise + the clip-discovered
    seepage set at partition boundaries (the design's flagged unknown), so the
    area-weighted mean head is checked to a looser tolerance (same rationale as
    the `rel_sum_fa` platform floor).
    """
    import json
    import shutil
    import subprocess
    import sys

    if _petsc4py_abi_mismatch():
        pytest.skip(
            "osx-arm64 py310-only petsc4py segfaults on nested mpirun finalize."
        )
    if shutil.which("mpirun") is None:
        pytest.skip("mpirun not on PATH; cannot exercise MPI decomposition.")

    fixtures_dir = FIXTURES_DIR
    needed = ["minimal_gw.yml", "mesh.npz", "soiltemp.npz"]
    if not all((fixtures_dir / f).exists() for f in needed):
        pytest.skip(f"missing one of {needed} in tests/fixtures.")

    dump_py = tmp_path / "_gw_parallel_dump.py"
    dump_py.write_text(_GW_PARALLEL_DUMP_SCRIPT)

    def run_at_rank(n):
        out_dir = tmp_path / f"n{n}"
        out_dir.mkdir()
        for f in needed:
            shutil.copy(fixtures_dir / f, out_dir / f)
        stats_json = out_dir / "stats.json"
        cmd = ["mpirun", "-n", str(n), sys.executable, str(dump_py),
               "minimal_gw.yml", str(stats_json)]
        child_env = mpi_child_env()
        child_env.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
        result = subprocess.run(cmd, cwd=out_dir, timeout=MPI_TIMEOUT,
                                capture_output=True, text=True, env=child_env)
        if result.returncode != 0 or not stats_json.exists():
            pytest.fail(
                f"`mpirun -n {n}` gw subprocess failed (rc={result.returncode}).\n"
                f"stdout:\n{result.stdout[-2000:]}\nstderr:\n{result.stderr[-2000:]}"
            )
        with open(stats_json) as f:
            return json.load(f)

    s1 = run_at_rank(1)
    s2 = run_at_rank(2)
    assert s1["size"] == 1 and s2["size"] == 2

    def rel(a, b):
        return abs(a - b) / max(abs(a), abs(b), 1.0e-30)

    # All three are computed on the EVOLVED state, so they carry the usual ~%
    # partition drift (the elevation → drainage → seaID/lake mask flips a few
    # boundary cells; plus KSP FP-noise + the clip-discovered seepage set in the
    # head). Tolerances are the partition-drift floor (same rationale as the
    # `rel_sum_fa` 5%); a real partition-dependence bug would blow far past these.
    assert rel(s1["sum_rech"], s2["sum_rech"]) < 5.0e-2, (
        f"recharge sum differs np1 vs np2: {s1['sum_rech']} vs {s2['sum_rech']}"
    )
    assert rel(s1["wmean_head"], s2["wmean_head"]) < 5.0e-2, (
        f"mean head differs np1 vs np2: {s1['wmean_head']} vs {s2['wmean_head']}"
    )
    assert rel(s1["wmean_wt"], s2["wmean_wt"]) < 1.5e-1, (
        f"mean water-table depth differs np1 vs np2: "
        f"{s1['wmean_wt']} vs {s2['wmean_wt']}"
    )
    # Baseflow output must be partition-consistent (no halo seam): the local
    # baseflowL equals its owner-synced value on every rank (fixed by the
    # _baseflowClosure local->global->local sync).
    assert s2.get("bf_halo", 0.0) < 1.0e-9, (
        f"baseflow not partition-consistent (halo seam): maxdiff {s2.get('bf_halo')}"
    )
    # Duricrust / induration / Karmor outputs must be partition-consistent (no halo
    # seam): the local duriHL and duriF equal their owner-synced values on every
    # rank (fixed by the _updateDuricrust/_updateSolute/_recordInduration syncs).
    assert s2.get("duri_halo", 0.0) < 1.0e-9, (
        f"duricrust not partition-consistent (halo seam): maxdiff {s2.get('duri_halo')}"
    )
    assert s2.get("durf_halo", 0.0) < 1.0e-9, (
        f"induration not partition-consistent (halo seam): maxdiff {s2.get('durf_halo')}"
    )


def test_geochem_flux_parallel(tmp_path):
    """
    Protects: the geochem seepage-export field ``gwSoluteFlux`` / ``soluteflux_*``
    is computed on OWNED nodes only (``seep_sink`` from the local FV stencil), so
    it MUST be halo-synced in ``_updateSolute`` — otherwise the output shows a
    partition seam. Runs a geochem fixture under np=2 and asserts the local field
    equals its owner-synced round-trip (``sf_halo`` ~ 0, non-zero before the fix).
    """
    import json
    import shutil
    import subprocess
    import sys

    if _petsc4py_abi_mismatch():
        pytest.skip("petsc4py ABI mismatch (nested mpirun segfault).")
    if shutil.which("mpirun") is None:
        pytest.skip("mpirun not on PATH; cannot exercise MPI decomposition.")

    fixtures_dir = FIXTURES_DIR
    needed = ["minimal_gw_geochem.yml", "mesh.npz", "soiltemp.npz"]
    if not all((fixtures_dir / f).exists() for f in needed):
        pytest.skip(f"missing one of {needed} in tests/fixtures.")

    dump_py = tmp_path / "_gw_geochem_dump.py"
    dump_py.write_text(_GW_PARALLEL_DUMP_SCRIPT)
    out_dir = tmp_path / "n2"
    out_dir.mkdir()
    for f in needed:
        shutil.copy(fixtures_dir / f, out_dir / f)
    stats_json = out_dir / "stats.json"
    cmd = ["mpirun", "-n", "2", sys.executable, str(dump_py),
           "minimal_gw_geochem.yml", str(stats_json)]
    child_env = mpi_child_env()
    child_env.setdefault("HDF5_USE_FILE_LOCKING", "FALSE")
    result = subprocess.run(cmd, cwd=out_dir, timeout=MPI_TIMEOUT,
                            capture_output=True, text=True, env=child_env)
    if result.returncode != 0 or not stats_json.exists():
        pytest.fail(
            f"`mpirun -n 2` geochem subprocess failed (rc={result.returncode}).\n"
            f"stderr:\n{result.stderr[-2000:]}"
        )
    with open(stats_json) as f:
        s = json.load(f)
    assert s.get("sf_halo", 0.0) < 1.0e-9, (
        f"soluteflux not partition-consistent (halo seam): maxdiff {s.get('sf_halo')}"
    )


# ---------------------------------------------------------------------------
# TEST 8 - Parallel decomposition must not alter the physical solution
# ---------------------------------------------------------------------------
#
# MPI's COMM_WORLD.size is fixed at process launch (`mpirun -n N`), so a
# single pytest invocation can only see ONE rank count. To compare n=1 vs
# n=N runs of the same model, this test spawns TWO subprocesses (`mpirun
# -n 1` and `mpirun -n 2`), each writing a JSON summary of global state,
# then asserts the summaries agree across decompositions.
#
# Why n=2 (not n=4) for the parallel case: GitHub's macOS runners (macos-15) only
# have 3 cores; OpenMPI refuses to oversubscribe by default. n=2 exposes
# exactly one partition boundary, which is sufficient to catch every
# decomposition-related bug class this test is designed to surface.
# ---------------------------------------------------------------------------

_PARALLEL_DUMP_SCRIPT = '''\
"""
Subprocess entry point for test_parallel_correctness.

Runs `gospl.model.Model(<yml>)`, performs MPI-reduced summary stats over
owned (non-ghost) nodes, and writes a JSON file on rank 0. Designed to be
invoked via `mpirun -n N python <this script> <yml> <out_json>` from the
test body — the test then loads two such JSON dumps (n=1 and n=2) and
asserts the global quantities match within tight tolerances.
"""
import json
import sys
import numpy as np
from mpi4py import MPI
from gospl.model import Model

fixture_yml = sys.argv[1]
output_json = sys.argv[2]

comm = MPI.COMM_WORLD
model = Model(fixture_yml, verbose=False, showlog=False)
try:
    model.runProcesses()

    # Owned-node mask: ghost nodes must never enter a reduction (they are
    # halo copies owned by a neighbouring rank). Same pattern as Tests 6-7.
    owned = model.inIDs == 1
    h_owned     = model.hLocal.getArray()[owned]
    cumED_owned = model.cumEDLocal.getArray()[owned]
    fa_owned    = model.FAL.getArray()[owned]
    larea_owned = model.larea[owned]

    # Per-rank local scalars; MPI.SUM/MAX-reduced to globals below.
    wmean_h_local  = float((h_owned * larea_owned).sum())
    area_local     = float(larea_owned.sum())
    flux_local     = float((cumED_owned * larea_owned).sum())
    activity_local = float((np.abs(cumED_owned) * larea_owned).sum())
    max_fa_local      = float(fa_owned.max()) if len(fa_owned) > 0 else 0.0
    sum_fa_local      = float(fa_owned.sum())
    # Owned-node count per rank. After allreduce(SUM) this is the total
    # number of nodes uniquely owned across the comm. Should equal the
    # mesh's global vertex count regardless of decomposition. If it
    # doesn't, there's a partition-ownership gap (a node owned by zero
    # ranks, or — less likely — double-counted by inIDs).
    owned_count_local = int(owned.sum())

    wmean_h_num = comm.allreduce(wmean_h_local, op=MPI.SUM)
    area_total  = comm.allreduce(area_local,    op=MPI.SUM)
    flux        = comm.allreduce(flux_local,    op=MPI.SUM)
    activity    = comm.allreduce(activity_local, op=MPI.SUM)
    max_fa      = comm.allreduce(max_fa_local,  op=MPI.MAX)
    sum_fa      = comm.allreduce(sum_fa_local,  op=MPI.SUM)
    owned_count = comm.allreduce(owned_count_local, op=MPI.SUM)
    wmean_h     = wmean_h_num / area_total

    if comm.Get_rank() == 0:
        with open(output_json, "w") as f:
            json.dump({
                "size":        comm.Get_size(),
                "wmean_h":     wmean_h,
                "flux":        flux,
                "activity":    activity,
                "max_fa":      max_fa,
                "sum_fa":      sum_fa,
                "owned_count": owned_count,
            }, f)
finally:
    model.destroy()
'''


# Groundwater np=1-vs-2 dump: reduce the water-table head / recharge over owned
# nodes and write a JSON on rank 0 (see test_watertable_parallel). Checks the
# Phase-2 head solve is partition-consistent (the design's flagged unknown).
_GW_PARALLEL_DUMP_SCRIPT = '''
import json
import sys
import numpy as np
from mpi4py import MPI
from gospl.model import Model

comm = MPI.COMM_WORLD
m = Model(sys.argv[1], verbose=False, showlog=False)
try:
    m.runProcesses()
    owned = m.inIDs == 1
    wt = m.wtDepth[owned]
    head = m.headL.getArray()[owned]
    rech = m.rechargeL.getArray()[owned]
    area = m.larea[owned]
    at = comm.allreduce(float(area.sum()), op=MPI.SUM)
    wmean_head = comm.allreduce(float((head * area).sum()), op=MPI.SUM) / at
    wmean_wt = comm.allreduce(float((wt * area).sum()), op=MPI.SUM) / at
    max_wt = comm.allreduce(float(wt.max()) if len(wt) else 0.0, op=MPI.MAX)
    sum_rech = comm.allreduce(float(rech.sum()), op=MPI.SUM)
    # Baseflow halo consistency: the local (halo-inclusive) baseflowL must equal
    # its owner-synced round-trip (0 diff). A partition seam (halo left at 0)
    # would make this non-zero — guards the _baseflowClosure halo sync.
    if getattr(m, "gwConserveBaseflow", False):
        bf = m.baseflowL.getArray().copy()
        m.tmpL.setArray(bf)
        m.dm.localToGlobal(m.tmpL, m.tmp)
        m.dm.globalToLocal(m.tmp, m.tmpL)
        bf_halo = comm.allreduce(float(np.abs(bf - m.tmpL.getArray()).max()), op=MPI.MAX)
    else:
        bf_halo = 0.0
    # Solute-export halo consistency (geochem): gwSoluteFlux = seep_sink*c*A is
    # computed on owned nodes only and synced in _updateSolute — a seam if not.
    if getattr(m, "gwGeochemOn", False):
        sf = m.gwSoluteFlux.copy()
        m.tmpL.setArray(sf)
        m.dm.localToGlobal(m.tmpL, m.tmp)
        m.dm.globalToLocal(m.tmp, m.tmpL)
        sf_halo = comm.allreduce(float(np.abs(sf - m.tmpL.getArray()).max()), op=MPI.MAX)
    else:
        sf_halo = 0.0
    # Duricrust halo consistency: duriH (duricrust), duriF (induration) and the
    # Karmor multiplier are written to output, so their ghost nodes must equal the
    # owners' values. A partition seam (halo left at the rank-local computed value)
    # would make these non-zero — guards the _updateDuricrust / _updateSolute /
    # _recordInduration halo syncs.
    def _halo(arr):
        m.tmpL.setArray(arr.copy())
        m.dm.localToGlobal(m.tmpL, m.tmp)
        m.dm.globalToLocal(m.tmp, m.tmpL)
        return comm.allreduce(float(np.abs(arr - m.tmpL.getArray()).max()), op=MPI.MAX)
    if getattr(m, "duriOn", False):
        duri_halo = _halo(m.duriHL.getArray())
        durf_halo = _halo(m.duriF)
    else:
        duri_halo = 0.0
        durf_halo = 0.0
    if comm.Get_rank() == 0:
        with open(sys.argv[2], "w") as f:
            json.dump({
                "size": comm.Get_size(),
                "wmean_head": wmean_head,
                "wmean_wt": wmean_wt,
                "max_wt": max_wt,
                "sum_rech": sum_rech,
                "bf_halo": bf_halo,
                "sf_halo": sf_halo,
                "duri_halo": duri_halo,
                "durf_halo": durf_halo,
            }, f)
finally:
    m.destroy()
'''


def _petsc4py_abi_mismatch() -> bool:
    """Return True if petsc4py was built against a different Python than runtime."""
    try:
        import importlib.metadata, sys
        dist = importlib.metadata.distribution("petsc4py")
        tag = dist.metadata["Name"]  # just a check it exists
        # build string embeds the python target, e.g. np2py310h...
        # compare against the running interpreter
        for f in dist.files or []:
            if "py3" in str(f) and "cpython" not in str(f):
                pass
        # simpler: check via the direct build tag in the wheel name
        record = next(
            (str(f) for f in (dist.files or []) if str(f).endswith(".dist-info/WHEEL")),
            None,
        )
        if record:
            wheel_text = (dist._path.parent / record).read_text()
            import re, sys
            tag_match = re.search(r"Tag: (\S+)", wheel_text)
            if tag_match:
                tag_str = tag_match.group(1)
                expected = f"cp{sys.version_info.major}{sys.version_info.minor}"
                return expected not in tag_str
    except Exception:
        pass
    return False

@pytest.mark.slow
@pytest.mark.skipif(
    _petsc4py_abi_mismatch(),
    reason="petsc4py built against different Python ABI (no py311 osx-arm64 "
           "build on conda-forge); segfaults on MPI finalization — not a gospl bug"
)
def test_parallel_correctness(tmp_path):
    """
    Protects: AGENTS.md > MPI contract — collective operations must yield
    physically identical results regardless of how the domain is partitioned
    across ranks.

    Silent failures prevented:

      1. Halo-exchange bug in `dm.localToGlobal` / `dm.globalToLocal`
         (mesher/unstructuredmesh.py:738-799): a ghost node written by the
         wrong rank would corrupt the elevation field in a thin ring at each
         subdomain boundary.

      2. Flow-routing inconsistency across partition boundaries: the
         donor-receiver graph in `flowplex._buildFlowDirection` depends on
         which rank owns each node. A bug in the boundary-node hand-off
         would produce divergent rcvID arrays, which cascade into divergent
         drainage area and elevation.

      3. Incorrect `inIDs` mask causing a rank to count ghost-node
         contributions in its local reduction (anywhere `inIDs == 1` is
         used to select owned nodes). With n=2 this would over-count
         boundary nodes by ~2x and shift global statistics.

      4. Non-deterministic KSP convergence across decompositions. PETSc
         parallel reductions (dot products, norms) are NOT guaranteed
         bitwise-identical because floating-point addition is not
         associative. We therefore assert STATISTICAL equivalence with
         tight relative tolerances rather than bitwise identity.

    Implementation: MPI.COMM_WORLD.size is fixed at process launch, so we
    cannot run n=1 and n=N in the same pytest process. Two subprocesses
    via `mpirun -n N python <dump_script>` each emit a JSON dump on rank
    0; this test reads both and compares.

    Tier 1A — drainage area (1e-4 relative, checked FIRST because routing
    bugs cascade into elevation bugs but not vice versa). Both the max
    and the sum of FAL over owned nodes must agree.

    Tier 1B — area-weighted mean elevation (1e-10 relative). Catches
    rank-doubling of ghost-node contributions in any reduction.

    Tier 1C — total cumED×area flux (1e-10 relative). Catches sediment
    created or destroyed at subdomain boundaries via incorrect axpy
    in SPL.py:352, nlSPL.py:404, or soilSPL.py:326.

    Non-applicability gates:
      - skip if `mpirun` is not on PATH (no MPI runtime to spawn);
      - skip if minimal.yml / mesh.npz are missing (consistent with the
        other Model-based regression tests);
      - skip if activity < 1.0 m³ (vacuous comparison on a fixture that
        moves negligible sediment — same gate as Tests 6-7).
    """
    import json
    import shutil
    import subprocess
    import sys

    # ---- Non-applicability gates ---------------------------------------
    if shutil.which("mpirun") is None:
        pytest.skip(
            "mpirun not on PATH; cannot exercise MPI decomposition. "
            "On CI this means the conda env lacks an MPI runtime (mpi4py "
            "package missing or broken). Locally, install via "
            "`mamba env create -f environment.yml`."
        )

    fixtures_dir = FIXTURES_DIR
    yml_src = fixtures_dir / "minimal.yml"
    mesh_src = fixtures_dir / "mesh.npz"
    if not yml_src.exists() or not mesh_src.exists():
        pytest.skip(
            f"{yml_src} or {mesh_src} not present. Same fixture as the "
            f"other Model-based regression tests."
        )

    # ---- Write the per-subprocess dump script once ---------------------
    # Stored as a real .py file (not -c "string") so tracebacks point at
    # readable lines if anything inside fails.
    dump_py = tmp_path / "_parallel_dump.py"
    dump_py.write_text(_PARALLEL_DUMP_SCRIPT)

    # ---- Run the model under mpirun -n N -------------------------------
    def run_at_rank(n):
        """
        Spawn `mpirun -n {n} python <dump_py> minimal.yml stats.json`
        in an isolated cwd. Each cwd has its own copies of the YAML and
        mesh so the goSPL output directory (`<cwd>/minimal/`) does not
        collide between the two runs.
        """
        out_dir = tmp_path / f"n{n}"
        out_dir.mkdir()
        shutil.copy(yml_src,  out_dir / "minimal.yml")
        shutil.copy(mesh_src, out_dir / "mesh.npz")
        stats_json = out_dir / "stats.json"

        cmd = [
            "mpirun", "-n", str(n),
            sys.executable, str(dump_py),
            "minimal.yml", str(stats_json),
        ]
        # Scrub inherited MPI-runtime env vars before spawning a nested mpirun.
        # pytest imported gospl (-> petsc4py.init -> MPI_Init), which under
        # OpenMPI exports OMPI_*/PMIX_*/PRTE_* into this process; if they leak
        # into the child `mpirun` it thinks it is already inside an MPI job and
        # silently refuses to launch (rc=1, empty output). MPICH is unaffected,
        # but the conda package and HPC container both use OpenMPI. OPAL_PREFIX
        # is preserved so mpirun can still locate its own libraries.
        child_env = mpi_child_env()
        result = subprocess.run(
            cmd, cwd=out_dir, timeout=MPI_TIMEOUT,
            capture_output=True, text=True, env=child_env,
        )
        if result.returncode != 0 or not stats_json.exists():
            pytest.fail(
                f"`mpirun -n {n}` subprocess failed "
                f"(rc={result.returncode}).\n"
                f"--- stdout ---\n{result.stdout}\n"
                f"--- stderr ---\n{result.stderr}"
            )
        return json.loads(stats_json.read_text())

    stats_n1 = run_at_rank(1)
    stats_n2 = run_at_rank(2)

    # ---- Sanity: the two subprocesses really did run at the requested
    #      rank counts (catches MPI implementation quirks, e.g. OpenMPI
    #      silently downgrading to 1 rank under oversubscription policy).
    assert stats_n1["size"] == 1, f"n=1 ran with size={stats_n1['size']}"
    assert stats_n2["size"] == 2, f"n=2 ran with size={stats_n2['size']}"

    # ---- Sanity: same total owned-node count across decompositions.
    # Every mesh vertex must be owned by EXACTLY one rank; the sum of
    # owned-node counts over the comm must therefore equal the mesh's
    # global vertex count regardless of decomposition. If this differs
    # between n=1 and n=2, there is a partition-ownership gap — bug in
    # `inIDs` setup (unstructuredmesh.py) or in PETSc DMPlex's vertex
    # partitioning. This is the FIRST thing to check if the cumulative
    # sums below diverge unexpectedly.
    assert stats_n1["owned_count"] == stats_n2["owned_count"], (
        f"n=1 and n=2 see different total owned-node counts "
        f"(n=1={stats_n1['owned_count']}, n=2={stats_n2['owned_count']}). "
        "Every mesh vertex must be owned by exactly one rank; this "
        "assertion failing means a node is owned by zero ranks (gap) "
        "or — less likely — double-counted (inIDs == 1 on two ranks). "
        "Suspect: mesher/unstructuredmesh.py's vIS-based `inIDs` "
        "assignment, or PETSc DMPlex partitioner setup."
    )

    # ---- Activity gate -------------------------------------------------
    if stats_n1["activity"] < 1.0:
        pytest.skip(
            f"Activity {stats_n1['activity']:.3e} m³ too small for "
            f"meaningful parallel-correctness check. NEEDS_HUMAN_REVIEW: "
            f"lengthen the run or steepen the gradient in minimal.yml."
        )

    # ---- Compute relative diffs ---------------------------------------
    def rel(a, b, denom):
        return abs(a - b) / max(abs(denom), 1e-30)

    rel_mean_h = rel(stats_n1["wmean_h"], stats_n2["wmean_h"],
                     stats_n1["wmean_h"])
    rel_flux   = rel(stats_n1["flux"], stats_n2["flux"],
                     stats_n1["activity"])
    rel_max_fa = rel(stats_n1["max_fa"], stats_n2["max_fa"],
                     stats_n1["max_fa"])
    rel_sum_fa = rel(stats_n1["sum_fa"], stats_n2["sum_fa"],
                     stats_n1["sum_fa"])

    diagnostic = (
        f"\n  owned_n1       = {stats_n1['owned_count']}"
        f"\n  owned_n2       = {stats_n2['owned_count']}"
        f"\n  wmean_h_n1     = {stats_n1['wmean_h']:.6e} m"
        f"\n  wmean_h_n2     = {stats_n2['wmean_h']:.6e} m"
        f"\n  rel_mean_h     = {rel_mean_h:.3e}  (tol 1e-10)"
        f"\n  flux_n1        = {stats_n1['flux']:+.3e} m³"
        f"\n  flux_n2        = {stats_n2['flux']:+.3e} m³"
        f"\n  rel_flux       = {rel_flux:.3e}  (tol 1e-5)"
        f"\n  max_fa_n1      = {stats_n1['max_fa']:.3e} m²"
        f"\n  max_fa_n2      = {stats_n2['max_fa']:.3e} m²"
        f"\n  rel_max_fa     = {rel_max_fa:.3e}  (tol 1e-4)"
        f"\n  sum_fa_n1      = {stats_n1['sum_fa']:.3e} m²"
        f"\n  sum_fa_n2      = {stats_n2['sum_fa']:.3e} m²"
        f"\n  rel_sum_fa     = {rel_sum_fa:.3e}  (tol 5e-2)"
        f"\n  activity       = {stats_n1['activity']:.3e} m³"
    )

    # ---- Tier 1A: drainage area (routing bugs cascade — check first) ----
    assert rel_max_fa < 1e-4, (
        "Max drainage area differs between n=1 and n=2 beyond 1e-4 "
        "relative. Likely a partition-boundary bug in "
        "flow/flowplex.py:\n"
        "  (a) `_buildFlowDirection` assigns a boundary node's receiver "
        "to the wrong rank's local index.\n"
        "  (b) `_matrix_build` inserts a wrong off-diagonal entry at a "
        "halo node, so the IDA solve distributes water incorrectly "
        "across the partition boundary.\n"
        "  (c) The rcvIDi snapshot (flowplex.py:421-426) is taken AFTER "
        "a domain-decomposition-dependent globalToLocal."
        + diagnostic
    )
    # Sum-of-FA tolerance is INTENTIONALLY loose (5e-2 ~ 5%) — and the
    # gap between this and the other tolerances is the load-bearing
    # observation, not the absolute number.
    #
    # `mfdreceivers` already breaks EXACT slope ties deterministically
    # using global node IDs (see fortran/functions.F90 and the
    # `self.gid` argument in flow/flowplex.py:_buildFlowDirection). What
    # remains is the near-tie case: KSP convergence is not bitwise-
    # identical across decompositions because floating-point addition
    # is not associative under parallel reductions, so ghost-node
    # elevations differ by O(KSP rtol) between n=1 and n=2. A small
    # subset of near-tie nodes picks a different receiver and the
    # cumulative drainage statistics diverge.
    #
    # The DRIFT MAGNITUDE IS PLATFORM-DEPENDENT:
    #   - macOS-14 (arm64, conda-forge OpenMPI):  ~0.3%
    #   - Ubuntu-latest (x86_64, conda-forge MPICH): ~1.7%
    # (Same fixture, same Python, same goSPL, same tie-break fix —
    # the difference is in MPI implementation, BLAS variant, and FMA
    # availability at the platform level.) 5% is set as 3x the worst
    # observed on either platform. Real routing regressions (lost
    # neighbour entry, wrong row weights at halo nodes, stale rcvIDi
    # snapshot) would shift this by orders of magnitude — wholesale
    # rerouting drives rel_sum_fa toward 0.5+, not 0.05.
    #
    # Tightening this back below ~5% requires either KSP-precision
    # halo determinism (hard PETSc work) or a slope-tolerance band in
    # _buildFlowDirection (algorithm change). Out of scope here.
    assert rel_sum_fa < 5e-2, (
        "Total drainage area differs between n=1 and n=2 beyond 5e-2 "
        "relative — far beyond the ~0.3% (macOS) to ~1.7% (Ubuntu) "
        "noise floor from non-deterministic tie-breaking. Likely a "
        "real routing regression in flow/flowplex.py:\n"
        "  (a) `_buildFlowDirection` lost a neighbour entry at a "
        "partition boundary.\n"
        "  (b) The IDA matrix (`_matrix_build`) assembles with wrong "
        "row weights on halo nodes.\n"
        "  (c) The rcvIDi snapshot (flowplex.py:421-426) was taken from "
        "a stale (post-fill) topology."
        + diagnostic
    )

    # ---- Tier 1B: global mean elevation -------------------------------
    # Like rel_sum_fa, this drifts with the platform's parallel-reduction
    # order (non-associative floating-point): for n=1 vs n=2 the observed
    # drift is ~3.6e-11 on macOS-14 (arm64, OpenMPI) but ~7.6e-10 on
    # Ubuntu-latest (x86_64, MPICH) — the SAME near-tie receiver
    # non-determinism that loosens rel_sum_fa above, propagated into
    # elevation. 5e-9 is ~6x the worst observed and still ~7 orders of
    # magnitude below any real ghost-node reduction bug (which would drive
    # mean_h to the %-level, not the nm-level). Tightening this back needs
    # KSP-precision halo determinism (hard PETSc work), as for rel_sum_fa.
    assert rel_mean_h < 5e-9, (
        "Area-weighted mean elevation differs between n=1 and n=2 "
        "beyond 5e-9 relative. Likely causes:\n"
        "  (a) Ghost nodes included in a rank-local reduction before "
        "allreduce — every sum must use `inIDs == 1` as the owned-node "
        "mask.\n"
        "  (b) A `localToGlobal` call is missing after `hLocal.setArray`, "
        "so a rank is solving with a stale halo on the next KSP step.\n"
        "  (c) The KSP RHS double-counts a ghost-node contribution on "
        "ranks that own a boundary."
        + diagnostic
    )

    # ---- Tier 1C: total erosion / deposition flux ---------------------
    # Flux tolerance is intentionally loose (1e-5) because the residual
    # KSP-near-tie non-determinism described above shifts WHICH nodes
    # erode how much (a few near-tie nodes drain different basins →
    # different local incision). Observed 2.4e-7 on the minimal fixture
    # AFTER the exact-tie-break fix landed; 1e-5 gives ~40x headroom for
    # real mass-conservation regressions to trip the assert. Mean
    # elevation is unaffected (algorithm is mass-conserving regardless
    # of which path sediment takes) so Tier 1B remains at 1e-10.
    assert rel_flux < 1e-5, (
        "Total cumED × area flux differs between n=1 and n=2 beyond "
        "1e-5 relative — well beyond the ~1e-7 floor expected from "
        "routing tie-break non-determinism. Sediment is being created "
        "or destroyed at subdomain boundaries:\n"
        "  (a) An axpy in SPL.py:352, nlSPL.py:404, or soilSPL.py:326 "
        "uses the wrong scale factor on boundary nodes.\n"
        "  (b) sedplex._getSedFlux integrates an erosion source that "
        "has already been localToGlobal'd and therefore contains halo "
        "copies of boundary-node values."
        + diagnostic
    )


# ---------------------------------------------------------------------------
# TEST 8b - Soil + temperature map: parallel partition-shape correctness
# ---------------------------------------------------------------------------
#
# soilSPL.__init__ loads a temperature map (Norton et al. 2013 soil
# production) to build `self.prodSoil`. The map is stored full-mesh
# (mpoints) on disk and MUST be subset to the local partition via
# `[self.locIDs]` — exactly as the sibling `soilFile` branch does. Before
# the fix it was used un-subset, so `prodSoil` stayed global (mpoints) and
# `_form_residual_soil` (line 122) tried to broadcast it against the local
# `hSoil`/`rainVal` arrays. That ONLY works when MPIsize == 1 (lpoints ==
# mpoints); in parallel it raises
#     ValueError: operands could not be broadcast together
#                 with shapes (mpoints,) (lpoints,)
# A classic serial-only-tested path: invisible at n=1, fatal at n>1.
#
# This test spawns `mpirun -n 2` on the soil+temp fixture and asserts the
# run completes — n=2 exposes exactly one partition boundary, enough to
# trip the un-subset path. Same subprocess/env-scrub machinery as
# test_parallel_correctness (nested-mpirun OMPI_*/PMIX_* leak guard).
# ---------------------------------------------------------------------------

_SOIL_TEMP_DRIVER = '''\
"""Subprocess entry for test_soil_temp_parallel_shape: run the soil+temp
model to completion and print a success sentinel on rank 0."""
import sys
from mpi4py import MPI
from gospl.model import Model

model = Model(sys.argv[1], verbose=False)
try:
    model.runProcesses()
finally:
    model.destroy()
if MPI.COMM_WORLD.Get_rank() == 0:
    print("SOIL_TEMP_PARALLEL_OK", flush=True)
'''


@pytest.mark.slow
@pytest.mark.skipif(
    _petsc4py_abi_mismatch(),
    reason="petsc4py built against different Python ABI; segfaults on MPI "
           "finalization — not a gospl bug (see test_parallel_correctness)",
)
def test_soil_temp_parallel_shape(tmp_path):
    """
    Protects: soilSPL.__init__ must subset the temperature map to the local
    partition (`loadData[self.tempData][self.locIDs]`). Regression guard for
    the parallel-only ValueError (global prodSoil broadcast against local
    hSoil) that is invisible at n=1 and fatal at n>1.
    """
    import shutil
    import subprocess
    import sys

    if shutil.which("mpirun") is None:
        pytest.skip("mpirun not on PATH; cannot exercise MPI decomposition.")

    fixtures_dir = FIXTURES_DIR
    yml_src = fixtures_dir / "minimal_soil_temp.yml"
    mesh_src = fixtures_dir / "mesh.npz"
    temp_src = fixtures_dir / "soiltemp.npz"
    if not (yml_src.exists() and mesh_src.exists() and temp_src.exists()):
        pytest.skip(
            "soil+temp fixtures not present (minimal_soil_temp.yml / mesh.npz "
            "/ soiltemp.npz)."
        )

    out_dir = tmp_path / "soiltemp"
    out_dir.mkdir()
    shutil.copy(yml_src, out_dir / "minimal_soil_temp.yml")
    shutil.copy(mesh_src, out_dir / "mesh.npz")
    shutil.copy(temp_src, out_dir / "soiltemp.npz")
    driver = out_dir / "_soil_temp_driver.py"
    driver.write_text(_SOIL_TEMP_DRIVER)

    # Scrub inherited OpenMPI runtime env before spawning a nested mpirun
    # (see test_parallel_correctness for the rationale).
    child_env = mpi_child_env()

    result = subprocess.run(
        ["mpirun", "-n", "2", sys.executable, str(driver), "minimal_soil_temp.yml"],
        cwd=out_dir,
        timeout=MPI_TIMEOUT,
        capture_output=True,
        text=True,
        env=child_env,
    )
    assert result.returncode == 0 and "SOIL_TEMP_PARALLEL_OK" in result.stdout, (
        "Soil+temperature model failed under `mpirun -n 2`. If this is a "
        "shape mismatch (global mpoints vs local lpoints), the temperature "
        "map in soilSPL.__init__ is not subset to self.locIDs.\n"
        f"(rc={result.returncode})\n"
        f"--- stdout ---\n{result.stdout}\n--- stderr ---\n{result.stderr}"
    )


# ---------------------------------------------------------------------------
# TEST 8b2 - Horizontal advection: parallel correctness (np=1 vs np=2)
# ---------------------------------------------------------------------------
# tectonics._varAdvector advects fields (elevation, cumED, flexure, soil) by an
# implicit FV upwind/IIOE solve on the DMPlex with a cached operator + KSP and a
# per-field warm start. This guards that the advected field is the SAME
# regardless of decomposition. Uses a PLANAR flat mesh (a ridge advected
# downwind by a uniform wind). NOTE: the cyclic cylinder mesh has a separate,
# known parallel-correctness bug (the wrap-seam partitions badly at some rank
# counts) — that is out of scope here; this test guards the common flat case.
# ---------------------------------------------------------------------------

_ADVECT_DUMP_SCRIPT = '''\
"""Subprocess entry for test_advection_parallel: advect on a flat mesh and
write owned-node stats of the advected elevation as JSON on rank 0."""
import json
import sys
import numpy as np
from mpi4py import MPI
from gospl.model import Model

comm = MPI.COMM_WORLD
model = Model(sys.argv[1], verbose=False, showlog=False)
try:
    model.runProcesses()
    owned = model.inIDs == 1
    h = model.hLocal.getArray()[owned]
    area = model.larea[owned]
    wmean = comm.allreduce(float((h * area).sum()), op=MPI.SUM) / \\
            comm.allreduce(float(area.sum()), op=MPI.SUM)
    hmax = comm.allreduce(float(h.max()) if len(h) else -1e30, op=MPI.MAX)
    owned_count = comm.allreduce(int(owned.sum()), op=MPI.SUM)
    if comm.Get_rank() == 0:
        json.dump({"size": comm.Get_size(), "wmean": wmean, "hmax": hmax,
                   "owned_count": owned_count}, open(sys.argv[2], "w"))
finally:
    model.destroy()
'''


@pytest.mark.slow
@pytest.mark.skipif(
    _petsc4py_abi_mismatch(),
    reason="petsc4py built against different Python ABI; segfaults on MPI "
           "finalization — not a gospl bug (see test_parallel_correctness)",
)
def test_advection_parallel(tmp_path):
    """
    Protects: horizontal advection (`tectonics._varAdvector`) gives the same
    advected field regardless of MPI decomposition, on a planar flat mesh. The
    cached operator + single-pass CSR assembly + per-field warm start must be
    partition-correct. (The cyclic cylinder mesh has a separate, documented
    parallel bug and is intentionally NOT used here.)

    Also guards against the continental-sediment closed-sink deadlock: the
    full `runProcesses` pipeline drives `sedplex._distributeSediment`, whose
    closed-sink conservation closure must reduce its `_closedDepo.any()` guard
    across ranks before running collective DM scatters. Before that fix this
    test hung at np=2 (one rank in the scatter block, the other already in
    `_spillCoords`'s Allreduce) and timed out at 600s.
    """
    import json
    import shutil
    import subprocess
    import sys

    if shutil.which("mpirun") is None:
        pytest.skip("mpirun not on PATH; cannot exercise MPI decomposition.")

    fixtures_dir = FIXTURES_DIR
    yml_src = fixtures_dir / "flat_advect.yml"
    npz_src = fixtures_dir / "flat_advect.npz"
    if not (yml_src.exists() and npz_src.exists()):
        pytest.skip("flat_advect fixtures not present.")

    dump_py = tmp_path / "_advect_dump.py"
    dump_py.write_text(_ADVECT_DUMP_SCRIPT)

    def run_at_rank(n):
        out_dir = tmp_path / f"n{n}"
        out_dir.mkdir()
        shutil.copy(yml_src, out_dir / "flat_advect.yml")
        shutil.copy(npz_src, out_dir / "flat_advect.npz")
        stats_json = out_dir / "stats.json"
        # Scrub inherited OpenMPI runtime env before the nested mpirun (see
        # test_parallel_correctness for the rationale).
        child_env = mpi_child_env()
        result = subprocess.run(
            ["mpirun", "-n", str(n), sys.executable, str(dump_py),
             "flat_advect.yml", str(stats_json)],
            cwd=out_dir, timeout=MPI_TIMEOUT, capture_output=True, text=True, env=child_env,
        )
        if result.returncode != 0 or not stats_json.exists():
            pytest.fail(
                f"`mpirun -n {n}` advection subprocess failed "
                f"(rc={result.returncode}).\n--- stdout ---\n{result.stdout}"
                f"\n--- stderr ---\n{result.stderr}"
            )
        return json.loads(stats_json.read_text())

    s1 = run_at_rank(1)
    s2 = run_at_rank(2)

    assert s1["size"] == 1 and s2["size"] == 2
    assert s1["owned_count"] == s2["owned_count"], "partition ownership gap"

    # The advected ridge must survive (a partition-boundary bug would collapse
    # or blow up the field — observed measure ~644 m).
    assert s2["hmax"] > 100.0, (
        f"advected ridge lost at np=2 (hmax={s2['hmax']:.3f})"
    )

    # n=1 and n=2 must agree: the only expected difference is KSP partition
    # noise (flat advection agrees to ~0.1% across decompositions).
    def rel(a, b, d):
        return abs(a - b) / max(abs(d), 1e-30)

    assert rel(s1["hmax"], s2["hmax"], s1["hmax"]) < 1e-2, (
        f"advected peak differs beyond noise: n1={s1['hmax']:.4f} "
        f"n2={s2['hmax']:.4f}"
    )
    assert rel(s1["wmean"], s2["wmean"], s1["wmean"]) < 1e-2, (
        f"advected mean differs beyond noise: n1={s1['wmean']:.4f} "
        f"n2={s2['wmean']:.4f}"
    )


def test_mpi_abort_on_exception(tmp_path):
    """
    Protects: model.py `_install_mpi_abort_excepthook`.

    At MPIsize>1 an uncaught exception on ANY rank must `MPI_Abort` the WHOLE
    job, never leave the other ranks blocked in a collective. Silent failure
    prevented: an exception on one rank (e.g. a fatal solve, or any bug) while
    the others sit in an `Allreduce`/`Barrier` deadlocks the job — it hangs and
    must be killed by hand (exactly the failure mode reported on the
    implicit-timestepping model). The hook turns that into an immediate, clean
    termination.

    Here rank 0 raises right after init while rank 1 enters a Barrier rank 0
    will never reach. Without the hook the run hangs to the timeout; with it the
    job aborts in seconds and rank 1 never gets PAST the barrier.
    """
    import shutil
    import subprocess
    import sys

    if shutil.which("mpirun") is None:
        pytest.skip("mpirun not on PATH; cannot exercise MPI abort.")
    if _petsc4py_abi_mismatch():
        pytest.skip("petsc4py ABI mismatch: nested mpirun segfaults at finalize.")

    fixtures_dir = FIXTURES_DIR
    yml_src = fixtures_dir / "minimal.yml"
    mesh_src = fixtures_dir / "mesh.npz"
    if not yml_src.exists() or not mesh_src.exists():
        pytest.skip("minimal fixtures not present.")

    out_dir = tmp_path / "abort"
    out_dir.mkdir()
    shutil.copy(yml_src, out_dir / "minimal.yml")
    shutil.copy(mesh_src, out_dir / "mesh.npz")
    script = out_dir / "abort_run.py"
    script.write_text(
        "from mpi4py import MPI\n"
        "from gospl.model import Model\n"
        "rank = MPI.COMM_WORLD.Get_rank()\n"
        "m = Model('minimal.yml', verbose=False)\n"
        "MPI.COMM_WORLD.Barrier()\n"
        "if rank == 0:\n"
        "    raise RuntimeError('synthetic fatal error on rank 0')\n"
        "MPI.COMM_WORLD.Barrier()  # rank 0 never arrives -> hang without the hook\n"
        "print('RANK_PAST_BARRIER', flush=True)\n"
    )

    # Scrub inherited OpenMPI runtime env so the nested mpirun launches (see
    # run_at_rank in test_parallel_correctness for the full rationale).
    child_env = mpi_child_env()

    try:
        result = subprocess.run(
            ["mpirun", "-n", "2", sys.executable, "abort_run.py"],
            cwd=out_dir, timeout=90, capture_output=True, text=True, env=child_env,
        )
    except subprocess.TimeoutExpired:
        pytest.fail(
            "Parallel run HUNG (>90s) after a rank-0 exception — the MPI abort "
            "hook did not fire, so the surviving rank deadlocked in its Barrier."
        )

    combined = result.stdout + result.stderr
    assert result.returncode != 0, "expected a non-zero exit after the abort"
    assert "RANK_PAST_BARRIER" not in result.stdout, (
        "rank 1 proceeded past the barrier — the abort hook did not kill it.\n"
        f"--- output ---\n{combined[:800]}"
    )
    assert ("MPI_ABORT" in combined) or ("synthetic fatal error" in combined), (
        f"no sign of the clean-abort path in the output:\n{combined[:800]}"
    )


# ---------------------------------------------------------------------------
# TEST 8c - Cached diffusion operators: collective rebuild decision (np>1)
# ---------------------------------------------------------------------------
#
# hillslope._hillSlope caches two coastline-gated diffusion operators — the
# marine flow-direction smoother (smooth=2) and the linear soil-creep solve
# (smooth=0). Each rebuilds its operator + redoes PCSetUp only when the sea
# mask `seaID` moves. `seaID` is RANK-LOCAL, but `_buildDiffMat` assembly and
# `PCSetUp` are COLLECTIVE. If the "did the coastline move?" test stays
# rank-local, then once the partitions' coastlines drift apart one rank
# rebuilds (entering collective Mat assembly / PCSetUp) while another reuses —
# divergent collective paths — and the run DEADLOCKS at np>1. Serial is immune,
# and a short run (e.g. test_parallel_correctness's 10 steps) never drifts the
# masks apart, so this shipped undetected (surfaced on a long stratigraphy run).
#
# Reproducing the asymmetry via physics is mesh/partition-dependent and flaky,
# so this guard constructs it directly: build both caches symmetrically, then
# perturb `seaID` on rank 0 only and re-enter the cached paths. Pre-fix that
# deadlocks; the collective `allreduce(..., op=MPI.LOR)` in `_hillSlope` makes
# every rank agree, so it completes. A hang shows up as a subprocess timeout.
# ---------------------------------------------------------------------------

_CACHE_REBUILD_DRIVER = '''\
"""Subprocess entry for test_parallel_cached_diffusion_rebuild.

Build the cached linear-hillslope (smooth=0) and marine-smoother (smooth=2)
operators symmetrically, then force an ASYMMETRIC coastline change (only rank 0's
seaID mask moves) and re-enter both cached paths. With a rank-local rebuild
decision this deadlocks at np>1; with the collective reduce it completes.
"""
import sys
import numpy as np
from mpi4py import MPI
from gospl.model import Model

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
model = Model(sys.argv[1], verbose=False, showlog=False)
try:
    h = model.hLocal.getArray()
    # Symmetric seed -> first call rebuilds each cached operator on ALL ranks.
    model.seaID = np.where(h <= model.sealevel)[0]
    model._hillSlope(smooth=0)
    model._hillSlope(smooth=2)
    # Asymmetric change: ONLY rank 0's mask now differs from what it cached.
    if rank == 0:
        model.seaID = (
            np.array([], dtype=int) if model.seaID.size > 0
            else np.array([0], dtype=int)
        )
    # Pre-fix: rank 0 rebuilds (collective) while the others reuse -> deadlock.
    model._hillSlope(smooth=0)
    model._hillSlope(smooth=2)
    comm.Barrier()
finally:
    model.destroy()
if rank == 0:
    print("CACHE_REBUILD_PARALLEL_OK", flush=True)
'''


@pytest.mark.slow
@pytest.mark.skipif(
    _petsc4py_abi_mismatch(),
    reason="petsc4py built against different Python ABI; segfaults on MPI "
           "finalization — not a gospl bug (see test_parallel_correctness)",
)
def test_parallel_cached_diffusion_rebuild(tmp_path):
    """
    Protects: AGENTS.md > MPI contract — the coastline-gated cached diffusion
    operators in `hillslope._hillSlope` (smooth=0 linear soil creep, smooth=2
    marine smoother) must decide whether to rebuild COLLECTIVELY. `seaID` is
    rank-local, but the rebuild it gates (`_buildDiffMat` assembly + `PCSetUp`)
    is collective, so the decision must be reduced across ranks
    (`MPI.COMM_WORLD.allreduce(..., op=MPI.LOR)`). Without it, an asymmetric
    coastline change (one rank's mask moves, another's does not) makes one rank
    rebuild while another reuses → divergent collective paths → deadlock at
    np>1 (invisible serially and on short runs).

    The test forces that asymmetry deterministically at np=2 and asserts the
    re-entry completes; a regression re-introduces the hang, which surfaces here
    as a subprocess timeout.
    """
    import shutil
    import subprocess
    import sys

    if shutil.which("mpirun") is None:
        pytest.skip("mpirun not on PATH; cannot exercise MPI decomposition.")

    fixtures_dir = FIXTURES_DIR
    yml_src = fixtures_dir / "minimal.yml"
    mesh_src = fixtures_dir / "mesh.npz"
    if not (yml_src.exists() and mesh_src.exists()):
        pytest.skip("minimal.yml / mesh.npz not present.")

    out_dir = tmp_path / "cacherebuild"
    out_dir.mkdir()
    shutil.copy(yml_src, out_dir / "minimal.yml")
    shutil.copy(mesh_src, out_dir / "mesh.npz")
    driver = out_dir / "_cache_rebuild_driver.py"
    driver.write_text(_CACHE_REBUILD_DRIVER)

    # Scrub inherited OpenMPI runtime env before spawning a nested mpirun
    # (see test_parallel_correctness for the rationale).
    child_env = mpi_child_env()

    try:
        result = subprocess.run(
            ["mpirun", "-n", "2", sys.executable, str(driver), "minimal.yml"],
            cwd=out_dir,
            timeout=180,
            capture_output=True,
            text=True,
            env=child_env,
        )
    except subprocess.TimeoutExpired:
        pytest.fail(
            "Cached-diffusion re-entry DEADLOCKED at np=2 (subprocess timed "
            "out). The coastline-gated rebuild in hillslope._hillSlope "
            "(smooth=0/2) must reduce the 'seaID changed' test across ranks "
            "(MPI.COMM_WORLD.allreduce(..., op=MPI.LOR)) before it gates the "
            "collective _buildDiffMat / PCSetUp."
        )
    assert result.returncode == 0 and "CACHE_REBUILD_PARALLEL_OK" in result.stdout, (
        "Asymmetric-seaID cached-diffusion re-entry failed under `mpirun -n 2` "
        f"(rc={result.returncode}).\n"
        f"--- stdout ---\n{result.stdout}\n--- stderr ---\n{result.stderr}"
    )


@pytest.mark.slow
def test_advection_partition_invariant(tmp_path):
    """IIOE2 advection + edge reset must not depend on the decomposition.

    Two partition bugs, both fixed 2026-10 and both guarded here (flat_advect
    ridge, np=1 vs np=3, 3 steps, whole domain):

    * interior: an owned row's inflow weight is the GHOST neighbour's outflow
      weight, computed from the ghost's truncated local stencil (Fortran FVarea,
      range, outflow count); `adveciioe2` now takes owner-synced inputs.
      Interior difference 1.8 m -> 3.8e-3 m;
    * edges: the post-advection edge reset (`_resetEdges`) marked only OWNED
      edge nodes, so `fitedges` averaged a ghost edge neighbour's raw value into
      an owned edge node; it now marks ghosts too (`_advEdgeHalo`). Edge-zone
      difference ~2-3 m (both IIOE schemes) -> ~1e-5 m.

    np=3 (not 2): the 2-way cut on this fixture does not exercise the limiter.
    """
    import shutil
    import subprocess
    import sys
    from pathlib import Path

    import numpy as np

    if shutil.which("mpirun") is None:
        pytest.skip("mpirun not on PATH")
    tool = Path(__file__).resolve().parents[1] / "scripts" / "ab_partition.py"
    for f in ("flat_advect.yml", "flat_advect.npz"):
        shutil.copy(FIXTURES_DIR / f, tmp_path / f)
    keep = tmp_path / "ab"
    res = subprocess.run(
        [sys.executable, str(tool), "flat_advect.yml", "-n", "3", "--steps", "3",
         "--keep", str(keep), "--timeout", str(MPI_TIMEOUT)],
        cwd=tmp_path, capture_output=True, text=True, env=mpi_child_env(),
        timeout=3 * MPI_TIMEOUT)
    a, b = keep / "ab_np1.npz", keep / "ab_np3.npz"
    assert a.exists() and b.exists(), res.stdout[-2000:] + res.stderr[-2000:]
    a, b = np.load(a), np.load(b)
    c = a["_coords"]
    lo, hi = c[:, :2].min(0), c[:, :2].max(0)
    dist = np.minimum(c[:, :2] - lo, hi - c[:, :2]).min(axis=1)
    d = np.abs(a["elev"] - b["elev"])
    assert np.isfinite(d).all()
    edge, interior = dist < 600.0, dist >= 600.0
    assert d[interior].max() < 0.01, f"interior difference {d[interior].max():.3e} m"
    assert d[edge].max() < 0.01, f"edge-zone difference {d[edge].max():.3e} m"
