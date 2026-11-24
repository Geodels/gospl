"""
Lithology list (``strata: lithologies:``, docs/DESIGN_CARBONATES.md 2.7, phase
C0): more lithologies than the dual-lithology fine fraction, each with its own
pile, compaction, erodibility, diffusivity and grain-size behaviour.

The checks split an existing fraction into two identical lithologies: the run
must then be the dual run (totals and the fine pile to round-off), and each
lithology must be conserved on its own. The dual-lithology run itself is
unchanged (its own tests, and the bitwise baseline of the C0 refactor).
"""
from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys

import numpy as np
import pytest

from _helpers import FIXTURES_DIR, MPI_TIMEOUT, mpi_child_env

pytest.importorskip(
    "gospl._fortran", reason="goSPL Fortran extension not built"
)

pytestmark = [pytest.mark.strata]

GRAV = (
    "gravity:\n"
    "    runout: 3.0e5\n"
    "    runout_fine: 1.2e6\n"
    "    fraction: 0.5\n\n"
)

FINE2 = ("    lithologies:\n      - {name: fine2, size: fine, phi0: 0.63, z0: 1960.0, "
         "k_factor: 1.5, diff_factor: 2.0, bedrock_frac: 0.2}\n")
COARSE2 = ("    lithologies:\n      - {name: coarse2, size: coarse, phi0: 0.49, "
           "z0: 3700.0, k_factor: 1.0, diff_factor: 1.0, bedrock_frac: 0.3}\n")


def _case(tmp_path, name, split=None, extra=""):
    """minimal_dual.yml in tmp_path, optionally with a lithology split
    (``fine``: fine + identical fine2, the fine bedrock share halved;
    ``coarse``: an identical coarse-like coarse2 with half the coarse
    bedrock share) and extra blocks."""
    src = FIXTURES_DIR / "minimal_dual.yml"
    if not src.exists() or not (FIXTURES_DIR / "mesh.npz").exists():
        pytest.skip("minimal_dual.yml / mesh.npz not present")
    if not (tmp_path / "mesh.npz").exists():
        shutil.copy(FIXTURES_DIR / "mesh.npz", tmp_path / "mesh.npz")
    text = re.sub(r"dir: '[^']*'", f"dir: '{name}_out'", src.read_text())
    if split == "fine":
        text = text.replace("    dual: True\n", "    dual: True\n" + FINE2)
        text = text.replace("bedrock_coarse_frac: 0.6", "bedrock_coarse_frac: 0.8")
    elif split == "coarse":
        text = text.replace("    dual: True\n", "    dual: True\n" + COARSE2)
    if extra:
        text = text.replace("output:", extra + "output:", 1)
    (tmp_path / f"{name}.yml").write_text(text)
    return f"{name}.yml"


def _run(tmp_path, yml):
    """Run a case serially; return the final state needed by the checks."""
    from gospl.model import Model

    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        m = Model(yml, verbose=False, showlog=False)
        try:
            m.runProcesses()
            top = m.stratStep + 1
            out = {
                "h": m.hLocal.getArray().copy(),
                "H": m.stratH[:, :top].copy(),
                "fine": m._fineGroup([h[:, :top] for h in m.stratHL]).copy(),
                "piles": [h[:, :top].copy() for h in m.stratHL],
                "ero": m._litEroded.copy(),
                "dep": m._litDeposited.copy(),
                "litNb": m.litNb,
            }
        finally:
            m.destroy()
    finally:
        os.chdir(cwd)
    return out


def _rel(a, b):
    return float(np.max(np.abs(a - b)) / max(np.max(np.abs(a)), 1.0e-30))


def test_lithology_table_parse(tmp_path):
    """The dual fine fraction is lithology 0; ``strata: lithologies:``
    appends entries with their properties; bad entries are refused."""
    from gospl.model import Model

    yml = _case(tmp_path, "parse", split="coarse")
    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        m = Model(yml, verbose=False, showlog=False)
        try:
            assert m.litNb == 2 and m.litName == ["fine", "coarse2"]
            assert m.litFine == [True, False]
            assert m.litK[0] == 1.5 and m.litD[0] == 2.0      # dual keys
            assert abs(m.litBed[0] - 0.4) < 1e-12 and m.litBed[1] == 0.3
            assert len(m.stratHL) == 2 and m.stratHf is m.stratHL[0]
        finally:
            m.destroy()
        for bad, msg in (("size: medium", "size"), ("name: fine,", "already used"),
                         ("bedrock_frac: 0.3", "sum to more than 1")):
            text = (tmp_path / yml).read_text()
            if msg == "size":
                text = text.replace("size: coarse", bad)
            elif msg == "already used":
                text = text.replace("name: coarse2,", bad)
            else:
                text = text.replace("bedrock_frac: 0.3", "bedrock_frac: 0.7")
            (tmp_path / "bad.yml").write_text(text)
            with pytest.raises(ValueError, match=msg):
                Model("bad.yml", verbose=False, showlog=False)
    finally:
        os.chdir(cwd)


@pytest.mark.slow
@pytest.mark.parametrize("split", ["fine", "coarse"])
def test_lithology_split_matches_dual(tmp_path, split):
    """
    Protects: the routed composition of several lithologies (river routing,
    pit cascade, marine lockstep and balancing, gravity currents).

    A fraction split into two identical lithologies gives the dual run: same
    elevation and layer thickness, same fine pile (sum of the fine-like
    piles), and each lithology conserved (eroded = deposited, closed fixture
    with gravity currents on).
    """
    ref = _run(tmp_path, _case(tmp_path, "ref", extra=GRAV))
    spl = _run(tmp_path, _case(tmp_path, "spl", split=split, extra=GRAV))
    assert spl["litNb"] == 2
    assert _rel(ref["h"], spl["h"]) < 1.0e-9
    assert _rel(ref["H"], spl["H"]) < 1.0e-9
    assert _rel(ref["fine"], spl["fine"]) < 1.0e-9
    for l in range(2):
        err = abs(spl["ero"][l] - spl["dep"][l]) / max(spl["ero"][l], 1.0)
        assert err < 1.0e-3, f"lithology {l}: eroded {spl['ero'][l]:.4e}, " \
                             f"deposited {spl['dep'][l]:.4e}"
    # the two halves of the split behave identically
    if split == "fine":
        assert _rel(spl["piles"][0][:, 1:], spl["piles"][1][:, 1:]) < 1.0e-9


_PILE_SCRIPT = '''\
"""Child for test_lithology_partition_invariant (mpirun -n N)."""
import sys
import numpy as np
from mpi4py import MPI
from gospl.model import Model

m = Model(sys.argv[1], verbose=False, showlog=False)
try:
    m.runProcesses()
    own = m.inIDs == 1
    top = m.stratStep + 1
    rows = np.column_stack([m.locIDs[own].astype(float)]
                           + [h[own, :top] for h in m.stratHL])
    allr = MPI.COMM_WORLD.gather(rows, root=0)
    if MPI.COMM_WORLD.rank == 0:
        a = np.vstack(allr)
        np.save(sys.argv[2], a[np.argsort(a[:, 0])])
finally:
    m.destroy()
'''


@pytest.mark.mpi
@pytest.mark.slow
def test_lithology_partition_invariant(tmp_path):
    """Several lithologies (coarse split, gravity currents) give the same
    piles at np=1 and np=3, to the solver round-off."""
    if shutil.which("mpirun") is None:
        pytest.skip("mpirun not on PATH")
    yml = _case(tmp_path, "part", split="coarse", extra=GRAV)
    child = tmp_path / "_piles.py"
    child.write_text(_PILE_SCRIPT)
    out = {}
    for n in (1, 3):
        f = tmp_path / f"piles{n}.npy"
        r = subprocess.run(["mpirun", "-n", str(n), sys.executable, str(child), yml, str(f)],
                           cwd=tmp_path, capture_output=True, text=True,
                           env=mpi_child_env(), timeout=MPI_TIMEOUT)
        assert r.returncode == 0, r.stderr[-2000:]
        out[n] = np.load(f)
    a, b = out[1], out[3]
    assert np.array_equal(a[:, 0], b[:, 0])
    # skip the bedrock sentinel column of each pile (1e5 m scale)
    nlay = (a.shape[1] - 1) // 2
    for l in range(2):
        cols = slice(1 + l * nlay + 1, 1 + (l + 1) * nlay)
        d = np.max(np.abs(a[:, cols] - b[:, cols]))
        assert d < 1.0e-5, f"lithology {l}: np=1 vs np=3 differ by {d:.3e} m"


def test_marine_balancing_degenerate_priors(tmp_path):
    """
    Protects: the marine balancing of the lithologies (``_marineLitShare``)
    on degenerate routed mixes. A node whose routed mix holds only a
    lithology with no input, and priors of ~1e-310 (subnormal), made the
    former in-place rescaling overflow; the NaN shares collapsed a whole
    carbonate run to a flat surface (2026-11-24). The shares must stay finite
    in [0, 1] and every lithology must deposit its share of the group input.
    """
    from gospl.model import Model

    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        m = Model(_case(tmp_path, "deg", split="fine"), verbose=False, showlog=False)
        try:
            n = m.lpoints
            mdep = np.full(n, 0.5)
            m.depoFineFrac[:] = 0.6
            prior = np.zeros((n, m.litNb))
            prior[:, 0] = 1.0
            prior[::3, 0] = 0.0
            prior[::3, 1] = 1.0              # only the lithology with no input
            prior[1::5, 0] = 1.0e-310        # subnormal priors
            m._marDiffLit = prior
            casIn = np.zeros((n, m.litNb))
            casIn[:, 0] = 1.0                # fine2 has no input at all
            sed = np.full(n, 2.0)
            m._marineLitShare(mdep, casIn, sed, casIn[:, 0].copy())
            sh = np.asarray(m.litShare)
            assert np.all(np.isfinite(sh)) and sh.min() >= 0.0 and sh.max() <= 1.0 + 1e-12
            own = m.inIDs == 1
            vol = mdep * m.larea * m.depoFineFrac
            dep = [(vol * sh[l])[own].sum() for l in range(m.litNb)]
            assert dep[1] <= 1.0e-6 * dep[0]       # no input, (almost) nothing deposited
            assert np.isclose(dep[0] + dep[1], vol[own].sum(), rtol=1e-9)
        finally:
            m.destroy()
    finally:
        os.chdir(cwd)
