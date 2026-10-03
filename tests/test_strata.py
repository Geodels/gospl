"""
Stratigraphy and dual (coarse/fine) lithology.

Protects: AGENTS.md > Dual lithology.

Split out of the former tests/test_regression.py. Run this group alone with
`pytest -m strata`; see tests/README.md for the marker list.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

from _helpers import FIXTURES_DIR, _strata_parser

# Skip the whole module (rather than erroring at collection) when the
# goSPL runtime stack is not installed.
pytest.importorskip(
    "gospl.tools.inputparser",
    reason="goSPL runtime deps (petsc4py / ruamel.yaml / scipy) not installed",
)

pytestmark = [pytest.mark.strata]


def test_dual_lithology_opt_in():
    """
    Protects: DESIGN_DUAL_LITHOLOGY.md Phase 0 — dual lithology is an
    opt-in parsed in `_extraStrata` (continuation of `_readCompaction`).

    Silent failure prevented: a refactor dropping the `_extraStrata` call
    from `_readCompaction`, or flipping the default, would silently change
    the sediment model for every existing input file.

    Invariants:
      1. No `strata` block  → stratLith False; coarse curve == compaction
         curve (the dual-off path must stay bitwise-identical).
      2. `strata: dual: True` with stratigraphy on → stratLith True and the
         per-lithology parameters are parsed.
      3. `strata: dual: True` with stratigraphy OFF (stratNb == 0) → flag
         forced back to False (dual requires stratigraphy).
    """
    # ---- Case 1: no strata block → single-fraction, defaults mirror compaction
    parser = _strata_parser(stratNb=5)
    parser._extraStrata()
    assert parser.stratLith is False
    assert parser.phi0c == parser.phi0s and parser.z0c == parser.z0s, (
        "Dual-off coarse porosity curve must default to the single-fraction "
        "compaction curve so behaviour is unchanged."
    )

    # ---- Case 2: dual on with stratigraphy enabled → parsed
    parser = _strata_parser(stratNb=5)
    parser.input = {
        "strata": {
            "dual": True,
            "coarse": {"phi0": 0.45, "z0": 3000.0},
            "fine": {"phi0": 0.65, "z0": 1500.0, "k_factor": 1.5},
            "bedrock_coarse_frac": 0.7,
            "fine_efficiency": 0.3,
            "pitInletBias": {"coarse": 0.8, "fine": 0.1},
            "fine_diff_factor": 2.0,
            "bedrock_sentinel": True,
        }
    }
    parser._extraStrata()
    assert parser.stratLith is True
    assert parser.phi0c == 0.45 and parser.z0c == 3000.0
    assert parser.phi0f == 0.65 and parser.z0f == 1500.0
    assert parser.fine_k_factor == 1.5
    assert parser.bedrock_coarse_frac == 0.7
    assert parser.fine_efficiency == 0.3
    assert parser.pit_inlet_bias_coarse == 0.8
    assert parser.pit_inlet_bias_fine == 0.1
    assert parser.fine_diff_factor == 2.0
    assert parser.bedrock_sentinel is True

    # ---- Case 3: dual requested but stratigraphy off → forced False
    parser = _strata_parser(stratNb=0)
    parser.input = {"strata": {"dual": True}}
    parser._extraStrata()
    assert parser.stratLith is False, (
        "Dual lithology must require stratigraphy (stratNb > 0); with strat "
        "disabled the flag has to fall back to single-fraction."
    )


def test_dual_lithology_layer_allocation():
    """
    Protects: DESIGN_DUAL_LITHOLOGY.md Phase 1 — the fine-fraction layer
    fields (stratHf, phiF) are allocated by readStratLayers only when
    dual lithology is enabled, and stay None otherwise.

    Silent failure prevented: allocating the fields unconditionally would
    change the single-fraction memory footprint and (via _outputStrat)
    write extra HDF5 datasets into every existing single-fraction run.

    Exercises the no-file branch of readStratLayers directly via __new__,
    which only needs lpoints/stratNb/phi0s/stratLith/memclear/strataFile.
    """
    from gospl.sed import stratplex

    from gospl.tools.constants import BEDROCK_SENTINEL

    def _strata_mesh(stratLith):
        m = stratplex.STRAMesh.__new__(stratplex.STRAMesh)
        m.strataFile = None
        m.lpoints = 8
        m.stratNb = 3
        m.phi0s = 0.49
        m.phi0f = 0.63
        m.bedrock_coarse_frac = 0.5
        m.memclear = False
        m.stratLith = stratLith
        m.stratHf = None
        m.phiF = None
        return m

    # Single-fraction: fine fields stay None.
    m = _strata_mesh(stratLith=False)
    m.readStratLayers()
    assert m.stratHf is None and m.phiF is None, (
        "Single-fraction run must not allocate stratHf/phiF."
    )

    # Dual: fine fields allocated with stratH's shape. The bedrock sentinel
    # layer (layer 0) carries the bedrock fine fraction; layers above are 0.
    m = _strata_mesh(stratLith=True)
    m.readStratLayers()
    assert m.stratHf is not None and m.phiF is not None
    assert m.stratHf.shape == m.stratH.shape == (m.lpoints, m.stratNb)
    assert m.phiF.shape == (m.lpoints, m.stratNb)
    assert np.allclose(m.stratHf[:, 0], BEDROCK_SENTINEL * (1.0 - 0.5)), (
        "Bedrock-sentinel layer must carry the (1 - bedrock_coarse_frac) fine split."
    )
    assert (m.stratHf[:, 1:] == 0.0).all(), "Layers above bedrock start coarse-empty."


def test_dual_lithology_initial_strata_composition(tmp_path):
    """
    Protects: a user can supply a per-layer coarse/fine composition for the
    INITIAL stratigraphy via the npstrata file (strataHf = fine bulk thickness,
    phiF = fine porosity). readStratLayers must load it per layer, clamp
    0 <= strataHf <= strataH, and default phiF to phi0f when absent.
    """
    stratplex = pytest.importorskip("gospl.sed.stratplex")
    n, nl = 3, 2
    f = tmp_path / "init.npz"
    np.savez(
        str(f),
        strataH=np.full((n, nl), 10.0),
        strataZ=np.zeros((n, nl)),
        phiS=np.full((n, nl), 0.49),
        # per-layer fine thickness; row 2 layer 0 is 15 > 10 -> must clamp to 10.
        strataHf=np.array([[2.0, 5.0], [0.0, 10.0], [15.0, 3.0]]),
        # phiF deliberately omitted -> should default to phi0f.
    )
    m = stratplex.STRAMesh.__new__(stratplex.STRAMesh)
    m.strataFile = str(f)
    m.lpoints = n
    m.mpoints = n
    m.stratNb = 2                     # extra capacity beyond the initial layers
    m.locIDs = np.arange(n)
    m.phi0s, m.phi0f = 0.49, 0.63
    m.memclear = False
    m.stratLith = True
    m.stratHf = None
    m.phiF = None

    m.readStratLayers()

    il = m.initLay
    assert il == nl
    # Per-layer composition loaded.
    assert np.isclose(m.stratHf[0, 0], 2.0) and np.isclose(m.stratHf[0, 1], 5.0)
    # Clamp: strataHf (15) exceeded strataH (10) -> clamped to 10.
    assert np.isclose(m.stratHf[2, 0], 10.0)
    # Partition stays physical everywhere.
    assert (m.stratHf[:, :il] <= m.stratH[:, :il] + 1e-9).all()
    assert (m.stratHf[:, :il] >= 0.0).all()
    # phiF absent -> defaulted to phi0f on the initial layers.
    assert np.allclose(m.phiF[:, :il], 0.63)


def test_strata_bedrock_sentinel(tmp_path):
    """
    Protects: strata.bedrock_sentinel — a dedicated infinite-bedrock sentinel
    layer is inserted BENEATH the file-provided initial layers. The file's
    layers shift up by one, layer 0 becomes the frozen 1e6 m reservoir with the
    bedrock_coarse_frac composition, and bedrockLay flips to 1 (so compaction
    freezes it). Off by default -> the legacy file path (no sentinel) is
    unchanged.
    """
    stratplex = pytest.importorskip("gospl.sed.stratplex")
    from gospl.tools.constants import BEDROCK_SENTINEL

    n, nl = 3, 2
    f = tmp_path / "init.npz"
    np.savez(
        str(f),
        strataH=np.full((n, nl), 10.0),
        strataZ=np.zeros((n, nl)),
        phiS=np.full((n, nl), 0.49),
        strataHf=np.full((n, nl), 4.0),     # each file layer 40% fine
    )

    def _mesh(sentinel):
        m = stratplex.STRAMesh.__new__(stratplex.STRAMesh)
        m.strataFile = str(f)
        m.lpoints = n
        m.mpoints = n
        m.stratNb = 2
        m.locIDs = np.arange(n)
        m.phi0s, m.phi0f = 0.49, 0.63
        m.bedrock_coarse_frac = 0.7         # -> 0.3 fine in the sentinel
        m.memclear = False
        m.stratLith = True
        m.stratHf = None
        m.phiF = None
        m.bedrock_sentinel = sentinel
        m.readStratLayers()
        return m

    # ---- sentinel ON: extra frozen bedrock layer at index 0 ----
    m = _mesh(True)
    assert m.initLay == nl + 1 and m.bedrockLay == 1
    # Layer 0 is the 1e6 m reservoir with the bedrock_coarse_frac split.
    assert np.allclose(m.stratH[:, 0], BEDROCK_SENTINEL)
    assert np.allclose(m.stratHf[:, 0], BEDROCK_SENTINEL * (1.0 - 0.7))
    assert np.allclose(m.phiS[:, 0], 0.49) and np.allclose(m.phiF[:, 0], 0.63)
    # File layers shifted to indices 1..nl, composition preserved.
    assert np.allclose(m.stratH[:, 1 : nl + 1], 10.0)
    assert np.allclose(m.stratHf[:, 1 : nl + 1], 4.0)

    # ---- sentinel OFF (default path): unchanged, no sentinel ----
    m0 = _mesh(False)
    assert m0.initLay == nl and m0.bedrockLay == 0
    assert np.allclose(m0.stratH[:, 0], 10.0)        # deepest file layer is layer 0
    assert np.allclose(m0.stratHf[:, 0], 4.0)


def test_strata_file_validation(tmp_path):
    """
    Protects: _checkStrataFile complains about a malformed npstrata file -- a
    missing required field or a layer array whose shape does not match
    'strataH' raises ValueError (fail fast) rather than a cryptic KeyError /
    broadcast error deep in the loader.
    """
    stratplex = pytest.importorskip("gospl.sed.stratplex")
    n, nl = 4, 2

    def _mesh(npz):
        m = stratplex.STRAMesh.__new__(stratplex.STRAMesh)
        m.strataFile = str(npz)
        m.lpoints = n
        m.mpoints = n
        m.locIDs = np.arange(n)
        m.phi0s, m.phi0f = 0.49, 0.63
        m.stratNb = 2
        m.memclear = False
        m.stratLith = False
        m.stratHf = None
        m.phiF = None
        return m

    # Missing required field 'phiS' -> ValueError naming it.
    f1 = tmp_path / "miss.npz"
    np.savez(str(f1), strataH=np.full((n, nl), 5.0), strataZ=np.zeros((n, nl)))
    with pytest.raises(ValueError, match="phiS"):
        _mesh(f1).readStratLayers()

    # Shape mismatch: phiS has the wrong number of layers.
    f2 = tmp_path / "shape.npz"
    np.savez(
        str(f2),
        strataH=np.full((n, nl), 5.0),
        strataZ=np.zeros((n, nl)),
        phiS=np.full((n, nl + 1), 0.49),
    )
    with pytest.raises(ValueError, match="phiS"):
        _mesh(f2).readStratLayers()

    # Wrong mesh dimension on strataH -> ValueError mentioning mesh_points.
    f3 = tmp_path / "mesh.npz"
    np.savez(
        str(f3),
        strataH=np.full((n + 1, nl), 5.0),
        strataZ=np.zeros((n + 1, nl)),
        phiS=np.full((n + 1, nl), 0.49),
    )
    with pytest.raises(ValueError, match="mesh_points"):
        _mesh(f3).readStratLayers()

    # Valid file -> loads cleanly (no exception), proving the gate is not
    # over-eager.
    f4 = tmp_path / "ok.npz"
    np.savez(
        str(f4),
        strataH=np.full((n, nl), 5.0),
        strataZ=np.zeros((n, nl)),
        phiS=np.full((n, nl), 0.49),
    )
    m = _mesh(f4)
    m.readStratLayers()
    assert m.initLay == nl


def test_dual_lithology_erosion_split():
    """
    Protects: DESIGN_DUAL_LITHOLOGY.md Phase 2 — erodeStrat splits the
    eroded solid into thCoarse/thFine by the consumed layers' composition,
    and stays mass-consistent. K-blend / composition helpers also tested.

    Builds a tiny single-node column by hand and drives erodeStrat through
    a stubbed PETSc boundary (tmpL/globalToLocal), so the split arithmetic
    is exercised without a full model.

    Invariants:
      1. All-coarse column -> thFine == 0 and thCoarse matches what the
         single-fraction branch would produce (parity).
      2. Mixed column -> the deposited (uncompacted) split reconstructs the
         eroded solid: thCoarse*(1-phi0c) + thFine*(1-phi0f) == solid removed.
      3. _surfaceComposition / _surfaceLithoK reflect the exposed layer and
         the fine_k_factor; both are neutral (1.0) when dual is off.
    """
    from gospl.sed import stratplex

    class _StubVec:
        def __init__(self, arr):
            self._a = arr
        def getArray(self):
            return self._a
        def setArray(self, a):
            self._a = np.asarray(a, dtype=np.float64)

    class _StubDM:
        def globalToLocal(self, g, l):
            l.setArray(g.getArray().copy())

    def _mesh(stratLith, ero, stratH, stratHf=None, phiS=None, phiF=None):
        m = stratplex.STRAMesh.__new__(stratplex.STRAMesh)
        lp, nb = stratH.shape
        m.lpoints, m.stratNb, m.stratStep = lp, nb, nb - 1
        m.dt = 1.0
        m.memclear = False
        m.phi0s = m.phi0c = 0.49
        m.phi0f = 0.63
        m.bedrock_coarse_frac = 0.5
        m.fine_k_factor = 1.0
        m.inIDs = np.ones(lp, dtype=int)
        m.larea = np.ones(lp)
        m._fineEroded = 0.0
        m.stratLith = stratLith
        m.stratH = stratH.astype(np.float64).copy()
        m.phiS = (phiS if phiS is not None else np.full_like(stratH, 0.49)).copy()
        m.stratK = np.ones_like(m.stratH)
        m.stratHf = stratHf.copy() if stratHf is not None else None
        m.phiF = phiF.copy() if phiF is not None else None
        # PETSc boundary stub: tmp carries -ero (erosion is negative).
        m.tmp = _StubVec(np.array(ero, dtype=np.float64))
        m.tmpL = _StubVec(np.zeros(lp))
        m.dm = _StubDM()
        return m

    # ---- Case 1: all-coarse parity (single vs dual must agree on thCoarse)
    stratH = np.array([[2.0, 3.0]])          # two layers, total 5 m
    single = _mesh(False, ero=[-4.0], stratH=stratH)
    single.erodeStrat()
    dual = _mesh(
        True, ero=[-4.0], stratH=stratH,
        stratHf=np.zeros_like(stratH), phiF=np.full_like(stratH, 0.63),
    )
    dual.erodeStrat()
    assert np.allclose(dual.thFine, 0.0), "All-coarse column must erode no fine."
    assert np.allclose(dual.thCoarse, single.thCoarse), (
        "All-coarse dual erosion must match the single-fraction thCoarse."
    )

    # ---- Case 2: mixed column -> per-fraction solid reconstruction
    stratH = np.array([[2.0, 3.0]])
    stratHf = np.array([[1.0, 1.2]])         # fine bulk per layer
    phiS = np.full_like(stratH, 0.49)
    phiF = np.full_like(stratH, 0.63)
    dual = _mesh(True, ero=[-4.0], stratH=stratH, stratHf=stratHf,
                 phiS=phiS, phiF=phiF)
    dual.erodeStrat()
    # Solid removed by erosion of 4 m: fully erode top layer (3 m) + 1 m of
    # the lower (well-mixed) layer.
    # top layer (idx1): coarse (3-1.2)*(1-.49) + fine 1.2*(1-.63)
    # lower (idx0): remove 1 m of 2 m -> half: coarse (2-1)/2*(1-.49)*1 ... compute generically:
    coarse_solid = (3 - 1.2) * (1 - 0.49) + ((2 - 1.0) * 0.5) * (1 - 0.49)
    fine_solid = 1.2 * (1 - 0.63) + (1.0 * 0.5) * (1 - 0.63)
    got = dual.thCoarse[0] * (1 - dual.phi0c) + dual.thFine[0] * (1 - dual.phi0f)
    assert np.isclose(got, coarse_solid + fine_solid, rtol=1e-9), (
        f"Per-fraction solid reconstruction failed: got {got}, "
        f"expected {coarse_solid + fine_solid}"
    )
    assert dual.thFine[0] > 0 and dual.thCoarse[0] > 0

    # ---- Case 3: composition + K-blend helpers
    m = _mesh(True, ero=[0.0], stratH=np.array([[2.0, 2.0]]),
              stratHf=np.array([[0.0, 0.5]]), phiF=np.full((1, 2), 0.63))
    fc = m._surfaceComposition()
    assert np.isclose(fc[0], 1.0 - 0.5 / 2.0), "Surface coarse fraction wrong."
    m.fine_k_factor = 2.0
    litK = m._surfaceLithoK()
    assert np.isclose(litK[0], fc[0] + (1 - fc[0]) * 2.0)
    m.stratLith = False
    assert np.allclose(m._surfaceLithoK(), 1.0), "K-blend must be neutral when off."


def test_dual_lithology_deposit_and_compaction():
    """
    Protects: DESIGN_DUAL_LITHOLOGY.md Phase 4 — per-fraction deposition
    (deposeStrat) and per-fraction compaction (_depthPorosityDual).

    Deposition: the fresh layer accumulates a fine fraction equal to the
    step's global eroded composition, with each lithology's surface porosity.

    Compaction (the headline physics): each fraction compacts on its own
    porosity-depth curve, conserving its solid phase while fines lose more
    bulk thickness than coarse at the same burial depth.
    """
    from gospl.sed import stratplex
    from mpi4py import MPI

    class _StubVec:
        def __init__(self, arr):
            self._a = np.asarray(arr, dtype=np.float64)
        def getArray(self):
            return self._a
        def setArray(self, a):
            self._a = np.asarray(a, dtype=np.float64)

    class _StubDM:
        def globalToLocal(self, g, l):
            l.setArray(g.getArray().copy())

    # ---- Deposition: 4 m deposit, per-node fine fraction (fineFrac) = 0.25 ----
    m = stratplex.STRAMesh.__new__(stratplex.STRAMesh)
    m.lpoints, m.stratNb, m.stratStep = 1, 2, 1
    m.stratLith = True
    m.memclear = False
    m.phi0c, m.phi0f = 0.49, 0.63
    m.depoFineFrac = np.array([0.25])   # per-node deposit fine fraction (Phase 3a/3b)
    m.inIDs = np.ones(1, dtype=int)
    m.larea = np.ones(1)
    m._fineDeposited = 0.0
    m.stratH = np.zeros((1, 2))
    m.stratHf = np.zeros((1, 2))
    m.phiS = np.zeros((1, 2))
    m.phiF = np.zeros((1, 2))
    m.stratK = np.ones((1, 2))
    m.tmp = _StubVec([4.0])
    m.tmpL = _StubVec([0.0])
    m.dm = _StubDM()
    m.deposeStrat()
    assert np.isclose(m.stratH[0, 1], 4.0)
    assert np.isclose(m.stratHf[0, 1], 4.0 * 0.25), (
        "Fine deposit must equal depo * per-node fineFrac."
    )
    assert np.isclose(m.phiF[0, 1], 0.63) and np.isclose(m.phiS[0, 1], 0.49)

    # ---- Compaction: one mixed layer buried 2 km ----
    m = stratplex.STRAMesh.__new__(stratplex.STRAMesh)
    m.stratStep = 0
    m.stratLith = True
    m.memclear = False
    m.bedrockLay = 0
    m.phi0c, m.z0c = 0.49, 3700.0
    m.phi0f, m.z0f = 0.63, 1960.0
    Hc0, Hf0 = 6.0, 4.0
    m.stratH = np.array([[Hc0 + Hf0]])
    m.stratHf = np.array([[Hf0]])
    m.phiS = np.array([[0.49]])
    m.phiF = np.array([[0.63]])
    depth = np.array([[-2000.0]])
    newH = m._depthPorosity(depth)

    phiS_new = 0.49 * np.exp(-2000.0 / 3700.0)
    phiF_new = 0.63 * np.exp(-2000.0 / 1960.0)
    # Per-fraction solid is conserved through compaction.
    Hc_new = newH[0, 0] - m.stratHf[0, 0]
    Hf_new = m.stratHf[0, 0]
    assert np.isclose(Hc_new * (1 - phiS_new), Hc0 * (1 - 0.49), rtol=1e-9), (
        "Coarse solid must be conserved through compaction."
    )
    assert np.isclose(Hf_new * (1 - phiF_new), Hf0 * (1 - 0.63), rtol=1e-9), (
        "Fine solid must be conserved through compaction."
    )
    assert newH[0, 0] < Hc0 + Hf0, "Compaction must reduce total thickness."
    # Fines lose proportionally more bulk than coarse at the same depth.
    assert (Hf_new / Hf0) < (Hc_new / Hc0), (
        "Fines must compact more than coarse for these curves."
    )


def test_dual_lithology_advection_fine_pile():
    """
    Protects: DESIGN_DUAL_LITHOLOGY.md Phase 5 — stratalRecord advects the
    fine pile (stratHf, phiF) alongside the total/coarse pile via a second
    strataonesed call (NOT stratathreesed, whose extra fields are 0-1
    fractions and renormalised — wrong for the bulk-thickness representation).

    Uses an identity advection (each node maps to itself with weight 1) so
    the records must come back unchanged, confirming the fine fields are
    actually routed through the interpolation and written back.
    """
    stratplex = pytest.importorskip("gospl.sed.stratplex")

    class _MockVec:
        def __init__(self, n):
            self._a = np.zeros(n)
        def setArray(self, a):
            self._a = np.asarray(a, dtype=np.float64).copy()
        def getArray(self):
            return self._a

    class _MockDM:
        def globalToLocal(self, src, dst):
            dst.setArray(src.getArray())  # identity halo exchange

    n = 4
    m = stratplex.STRAMesh.__new__(stratplex.STRAMesh)
    m.lpoints = n
    m.stratStep = 2          # advect layers 0..1
    m.stratLith = True
    m.stratH = np.array([[5.0, 7.0, 0.0]] * n)
    m.stratHf = np.array([[2.0, 3.0, 0.0]] * n)
    m.stratZ = np.array([[-10.0, -3.0, 0.0]] * n)
    m.phiS = np.array([[0.45, 0.48, 0.0]] * n)
    m.phiF = np.array([[0.60, 0.62, 0.0]] * n)
    m.tmp = _MockVec(n)
    m.tmpL = _MockVec(n)
    m.dm = _MockDM()

    # Identity interpolation: 3 neighbours all = self, weights summing to 1.
    indices = np.repeat(np.arange(n)[:, None], 3, axis=1)
    weights = np.full((n, 3), 1.0 / 3.0)
    onIDs = np.array([], dtype=int)

    H0, Hf0 = m.stratH.copy(), m.stratHf.copy()
    phiS0, phiF0 = m.phiS.copy(), m.phiF.copy()
    m.stratalRecord(indices, weights, onIDs)

    # Advected layers (0..stratStep-1) must be preserved by identity mapping.
    s = slice(0, m.stratStep)
    assert np.allclose(m.stratHf[:, s], Hf0[:, s]), (
        "Fine pile thickness must round-trip through identity advection."
    )
    assert np.allclose(m.phiF[:, s], phiF0[:, s]), (
        "Fine porosity must round-trip through identity advection."
    )
    # Coarse/total pile still correct (regression on the original behaviour).
    assert np.allclose(m.stratH[:, s], H0[:, s])
    assert np.allclose(m.phiS[:, s], phiS0[:, s])


def test_dual_lithology_pit_fine_bias():
    """
    Protects: DESIGN_DUAL_LITHOLOGY.md Phase 3b — _pitFineFraction biases the
    pit/lake deposit composition so fine concentrates toward the depocenter
    (deep) and coarse toward the inlet/margins (shallow), while conserving the
    per-pit incoming fine volume.

    Single pit, 4 nodes at increasing bathymetric depth, uniform deposit and
    area. The pit retained fine volume 1.2 of total 4 (ff_pit = 0.3, the
    coarse-settles-first retained fraction from the cascade). The resulting
    per-node fine fraction must (a) increase monotonically with depth and (b)
    conserve the retained fine volume: Σ(delta·larea·ffrac) == _pitRetFine.
    """
    sedplex = pytest.importorskip("gospl.sed.sedplex")
    n = 4
    m = sedplex.SEDMesh.__new__(sedplex.SEDMesh)
    m.lpoints = n
    m.stratLith = True
    m.pitParams = np.zeros((1, 3))          # one pit
    m.pitIDs = np.zeros(n, dtype=int)        # all nodes in pit 0
    m.inIDs = np.ones(n, dtype=int)
    m.larea = np.ones(n)
    m.lFill = np.full(n, 10.0)               # rim at 10 m
    hl = np.array([9.0, 7.0, 4.0, 8.0])      # depth = 1, 3, 6, 2
    delta = np.ones(n)                       # uniform deposit thickness
    depo = np.array([4.0])                   # retained volume (Σ delta*larea)
    m._pitRetFine = np.array([1.2])          # retained fine → ff_pit = 1.2/4 = 0.3
    m.depoFineFrac = np.zeros(n)
    # Default lacustrine bias: coarse delta vs fine depocenter (seg = 0.5).
    m.pit_inlet_bias_coarse = 0.5
    m.pit_inlet_bias_fine = 0.0

    m._pitFineFraction(hl, delta, depo)

    depth = m.lFill - hl
    ff = m.depoFineFrac
    # (a) fine fraction increases with depth (depocenter is fine-rich).
    order = np.argsort(depth)
    assert np.all(np.diff(ff[order]) > 0), (
        f"Fine fraction must increase with depth; got {ff} for depth {depth}."
    )
    # (b) per-pit fine volume conserved (ff_pit = 0.3, Σ delta*larea = 4).
    fine_vol = float(np.sum(delta * m.larea * ff))
    assert np.isclose(fine_vol, 0.3 * 4.0, rtol=1e-9), (
        f"Pit fine volume not conserved: {fine_vol} != 1.2"
    )

    # (c) the pitInletBias contrast controls the segregation STRENGTH: a larger
    # coarse-vs-fine contrast steepens the depth gradient, while equal biases
    # remove it (uniform composition) — all conserving the same fine volume.
    m.depoFineFrac = np.zeros(n)
    m.pit_inlet_bias_coarse = 1.0           # seg = 1.0 (full depth-proportional)
    m.pit_inlet_bias_fine = 0.0
    m._pitFineFraction(hl, delta, depo)
    ff_strong = m.depoFineFrac.copy()
    span_default = ff[order][-1] - ff[order][0]
    span_strong = ff_strong[order][-1] - ff_strong[order][0]
    assert span_strong > span_default, (
        "Stronger pitInletBias contrast must steepen the depth gradient."
    )
    assert np.isclose(float(np.sum(delta * m.larea * ff_strong)), 1.2, rtol=1e-9)

    m.depoFineFrac = np.zeros(n)
    m.pit_inlet_bias_coarse = 0.3           # seg = 0 → no compositional split
    m.pit_inlet_bias_fine = 0.3
    m._pitFineFraction(hl, delta, depo)
    assert np.allclose(m.depoFineFrac, 0.3, atol=1e-12), (
        "Equal coarse/fine inlet bias must give a uniform pit composition."
    )


def test_dual_lithology_marine_fine_bias():
    """
    Protects: DESIGN_DUAL_LITHOLOGY.md Phase 3c — _marineFineFraction biases
    the marine deposit composition so fine concentrates in deep / distal water
    and coarse stays proximal (shallow), conserving the marine fine volume.
    Subaqueous analogue of the pit-fine bias.

    Four marine nodes at increasing water depth, uniform deposit and area,
    uniform arriving composition (ff_mar = 0.3). The per-node fine fraction
    must increase with depth and conserve fine volume.
    """
    seaplex = pytest.importorskip("gospl.sed.seaplex")
    n = 4

    class _StubVec:
        def __init__(self, arr):
            self._a = np.asarray(arr, dtype=np.float64)
        def getArray(self):
            return self._a
        def setArray(self, a):
            self._a = np.asarray(a, dtype=np.float64)

    class _StubDM:
        def globalToLocal(self, g, l):
            l.setArray(g.getArray().copy())

    m = seaplex.SEAMesh.__new__(seaplex.SEAMesh)
    m.lpoints = n
    m.sealevel = 0.0
    m.inIDs = np.ones(n, dtype=int)
    m.larea = np.ones(n)
    m.depoFineFrac = np.zeros(n)
    mdep = np.ones(n)                        # uniform marine deposit
    m.tmp = _StubVec(mdep)
    m.tmpL = _StubVec(np.zeros(n))
    m.dm = _StubDM()
    hl = np.array([-1.0, -3.0, -6.0, -8.0])  # depth = 1, 3, 6, 8
    sedFlux = np.ones(n)                     # uniform marine input (total)
    fineFlux = np.full(n, 0.3)               # post-cascade fine → ff_mar = 0.3

    m._marineFineFraction(hl, sedFlux, fineFlux)

    depth = m.sealevel - hl
    ff = m.depoFineFrac
    order = np.argsort(depth)
    assert np.all(np.diff(ff[order]) > 0), (
        f"Marine fine fraction must increase with depth; got {ff}."
    )
    fine_vol = float(np.sum(mdep * m.larea * ff))
    assert np.isclose(fine_vol, 0.3 * n, rtol=1e-9), (
        f"Marine fine volume not conserved: {fine_vol} != {0.3 * n}"
    )


@pytest.mark.slow
def test_dual_model_runs_and_invariants(minimal_dual_model):
    """
    Integration (DESIGN_DUAL_LITHOLOGY.md Phase 6): a full dual-lithology model
    runs end-to-end and preserves the per-fraction invariants. Exercises the
    whole dual path together — erodeStrat split, _getSedFlux fine routing
    (Phase 3a), deposeStrat per-node fineFrac, per-fraction compaction, and
    fine-pile advection. This is the first live coverage of the dual sediment
    transport path (both fixtures used by other tests have stratNb == 0).
    """
    model = minimal_dual_model
    assert model.stratLith is True and model.stratNb > 0
    assert model.stratHf is not None and model.phiF is not None
    # Bedrock sentinel carries the configured fine split before any run.
    assert np.isclose(
        model.stratHf[0, 0], 1.0e6 * (1.0 - model.bedrock_coarse_frac)
    )

    model.runProcesses()

    top = model.stratStep + 1
    H = model.stratH[:, :top]
    Hf = model.stratHf[:, :top]
    # Fine pile stays physical: non-negative and never exceeding the total.
    assert (Hf >= -1.0e-9).all(), "Negative fine thickness in the strata pile."
    assert (Hf <= H + 1.0e-6).all(), "Fine thickness exceeds layer total."
    # Porosity and the routed fine fraction stay in range.
    assert (model.phiF[:, :top] >= -1.0e-12).all()
    assert (model.phiF[:, :top] <= 1.0 + 1.0e-12).all()
    assert (model.fineFrac >= 0.0).all() and (model.fineFrac <= 1.0).all()
    # The dual path actually moved fine material (bedrock contributes it).
    assert float(Hf.sum()) > 0.0


@pytest.mark.slow
def test_dual_lithology_surf_fine_frac_output(minimal_dual_model):
    """
    Protects: the dual-lithology surface composition field. A dual-lithology
    run must write ``surfFineFrac`` (the surface exposed mud share, 0-1) into
    the mesh HDF5, and reference it from the per-step XMF, so ``gospl.xdmf`` can
    be coloured by the in-place sand/mud composition in ParaView (the in-place
    complement to the ``sedLoadF`` flux).
    """
    import glob as _glob

    h5py = pytest.importorskip("h5py")

    m = minimal_dual_model
    assert m.stratLith
    m.runProcesses()
    out = m.outputDir

    # Surface fine fraction written into the mesh output, physical range.
    mh5 = sorted(_glob.glob(os.path.join(out, "h5", "%s.*.p0.h5" % m.file)))[-1]
    with h5py.File(mh5, "r") as hf:
        assert "surfFineFrac" in hf, "surfFineFrac missing from mesh output"
        sf = np.array(hf["surfFineFrac"])
        assert (sf >= -1e-6).all() and (sf <= 1.0 + 1e-6).all()

    # The per-step XMF advertises it so ParaView exposes the field.
    xmf = sorted(_glob.glob(os.path.join(out, "xmf", "%s*.xmf" % m.file)))[-1]
    assert "surfFineFrac" in open(xmf).read()


@pytest.mark.slow
def test_dual_all_coarse_matches_single_fraction(
    minimal_dual_coarse_model, minimal_strat_model
):
    """
    Parity guard (DESIGN_DUAL_LITHOLOGY.md Phase 6): dual lithology configured
    all-coarse (bedrock_coarse_frac=1.0, no erodibility/diffusivity contrast)
    must reproduce the single-fraction stratigraphy run exactly. Confirms the
    dual code paths are a faithful superset of the single-fraction path — any
    accidental divergence (extra deposition, wrong compaction branch, etc.)
    trips this.
    """
    dual = minimal_dual_coarse_model
    single = minimal_strat_model
    dual.runProcesses()
    single.runProcesses()

    top = min(dual.stratStep, single.stratStep) + 1
    # No fine produced anywhere in the all-coarse configuration.
    assert float(dual.stratHf[:, :top].sum()) == 0.0, (
        "All-coarse dual run produced fine sediment."
    )
    # Elevation and stratal thickness must match single-fraction bitwise
    # (the smoke check measured exactly 0 difference; atol guards float noise).
    assert np.allclose(
        dual.hLocal.getArray(), single.hLocal.getArray(), rtol=0.0, atol=1.0e-9
    ), "Elevation diverged from the single-fraction stratigraphy run."
    assert np.allclose(
        dual.stratH[:, :top], single.stratH[:, :top], rtol=0.0, atol=1.0e-9
    ), "Stratal thickness diverged from the single-fraction stratigraphy run."


@pytest.mark.slow
def test_dual_mass_conservation(minimal_dual_model):
    """
    Protects: DESIGN_DUAL_LITHOLOGY.md Phase 6 — sediment conservation through
    the dual-lithology pipeline on a closed sphere, and a valid per-fraction
    partition maintained end-to-end.

    The standard test_mass_conservation SKIPS when stratNb > 0 (compaction
    moves h without cumED), so dual/stratigraphy-mode total conservation is
    otherwise untested. cumED only changes via PAIRED sediment erosion/
    deposition (never by compaction), so on a closed sphere Σ(dcumED·larea)
    must still vanish relative to the redistributed volume — even with the
    extra fine-flux solve and the per-pit / marine composition reallocation
    that dual lithology adds. A real sediment leak (e.g. fine routed but not
    deposited) would push this to O(0.1+); measured ~1.6e-5 here.

    Per-fraction: the coarse/fine partition must stay valid through erosion,
    transport, deposition, compaction, advection and diffusion — fine bulk
    non-negative and never exceeding the layer total; the fine solid phase
    non-negative; and the pile must actually carry fine (not trivially zero).
    """
    model = minimal_dual_model

    # Closed-domain gate (mirrors test_mass_conservation).
    reasons = []
    if getattr(model, "flatModel", True):
        reasons.append("flatModel=True (boundary outflux)")
    if getattr(model, "tecdata", None) is not None:
        reasons.append("tectonics active")
    if getattr(model, "flexOn", False):
        reasons.append("flexure active")
    if reasons:
        pytest.skip(
            "Dual mass conservation requires a closed sphere with only "
            "sediment-conserving kernels. This fixture has: " + "; ".join(reasons)
        )
    assert model.stratLith and model.stratNb > 0, "fixture must enable dual strat"

    larea = model.larea
    owned = model.inIDs == 1
    cumED_before = model.cumEDLocal.getArray().copy()

    model.runProcesses()

    dED = model.cumEDLocal.getArray() - cumED_before
    from mpi4py import MPI
    dV = MPI.COMM_WORLD.allreduce(float(np.sum((dED * larea)[owned])), op=MPI.SUM)
    activity = MPI.COMM_WORLD.allreduce(
        float(np.sum((np.abs(dED) * larea)[owned])), op=MPI.SUM
    )

    # ---- Total sediment conserved in dual mode (closed sphere) ----
    assert activity > 0.0, "no sediment was redistributed; test is vacuous"
    rel = abs(dV) / activity
    assert rel < 5.0e-4, (
        f"Dual-mode sediment not conserved: |ΣdcumED|/activity = {rel:.2e} "
        f"(> 5e-4). A fraction routed-but-not-deposited would leak here."
    )

    # ---- Valid per-fraction partition end-to-end ----
    top = model.stratStep + 1
    H = model.stratH[:, :top]
    Hf = model.stratHf[:, :top]
    assert (Hf >= -1.0e-9).all(), "Negative fine bulk thickness."
    assert (Hf <= H + 1.0e-6).all(), "Fine bulk exceeds layer total."
    fine_solid = (Hf * (1.0 - model.phiF[:, :top]))[owned]
    assert (fine_solid >= -1.0e-9).all(), "Negative fine solid phase."
    fine_total = MPI.COMM_WORLD.allreduce(float(np.sum(fine_solid)), op=MPI.SUM)
    assert fine_total > 0.0, "Dual pile carries no fine — sub-system is dead."


@pytest.mark.slow
def test_dual_fine_conservation(minimal_dual_model):
    """
    Protects: per-fraction (FINE) mass conservation on a closed sphere.

    The total-cumED test_dual_mass_conservation checks the *total* sediment
    budget; it would pass even if fine were silently created/destroyed while
    coarse compensated. This test closes that gap directly: the fine solid
    REMOVED from the strata pile (erodeStrat, accumulated in _fineEroded) must
    equal the fine DEPOSITED back into it (deposeStrat continental + marine,
    _fineDeposited), to within the floor/transit budget.

    This is the guard required before any future fine-routing change (e.g.
    fine-enriched overspill): a mismatched fine mirror in the pit/marine
    cascade — exactly the bug that reverted the first overspill attempt —
    trips this where the total budget does not. Measured imbalance ~1e-4 on
    minimal_dual; asserted < 5e-3.
    """
    model = minimal_dual_model
    if getattr(model, "flatModel", True):
        pytest.skip("fine conservation needs a closed sphere (flatModel=False)")
    assert model.stratLith and model.stratNb > 0, "fixture must enable dual strat"

    model.runProcesses()

    from mpi4py import MPI
    ero = MPI.COMM_WORLD.allreduce(model._fineEroded, op=MPI.SUM)
    dep = MPI.COMM_WORLD.allreduce(model._fineDeposited, op=MPI.SUM)
    assert ero > 0.0, "no fine was eroded; test is vacuous"
    rel = abs(ero - dep) / ero
    assert rel < 5.0e-3, (
        f"Fine not conserved: eroded={ero:.4e} deposited={dep:.4e} "
        f"rel imbalance={rel:.2e} (> 5e-3). A fine-only leak the total-cumED "
        f"test cannot see."
    )


@pytest.mark.slow
def test_dual_sedloadf_output(minimal_dual_model, minimal_strat_model):
    """
    Protects: the fine sediment load `sedLoadF` is written to the HDF5 output
    when dual lithology is enabled, and absent for single-fraction runs.
    `sedLoad` is the TOTAL flux; `sedLoadF` the fine sub-flux, so
    0 <= sedLoadF <= sedLoad at every node.
    """
    h5py = pytest.importorskip("h5py")
    import glob

    def _latest_h5(model):
        files = sorted(
            glob.glob(os.path.join(str(model.outputDir), "h5", "gospl.*.p*.h5"))
        )
        return files[-1] if files else None

    # ---- dual run writes sedLoadF, 0 <= fine <= total ----
    dual = minimal_dual_model
    dual.runProcesses()
    f = _latest_h5(dual)
    assert f is not None, "no gospl HDF5 output was written"
    with h5py.File(f, "r") as hf:
        assert "sedLoadF" in hf, "sedLoadF missing from dual-lithology output"
        sl = np.array(hf["sedLoad"])[:, 0]
        slf = np.array(hf["sedLoadF"])[:, 0]
    assert (slf >= 0.0).all(), "negative fine sediment load"
    assert (slf <= sl + 1.0e-6).all(), "fine load exceeds total load"

    # ---- single-fraction stratigraphy run must NOT write sedLoadF ----
    single = minimal_strat_model
    single.runProcesses()
    fs = _latest_h5(single)
    with h5py.File(fs, "r") as hf:
        assert "sedLoad" in hf and "sedLoadF" not in hf, (
            "single-fraction output must not contain sedLoadF"
        )


# ---------------------------------------------------------------------------
# TEST 9 - Stratigraphy: deposition + compaction physics
# ---------------------------------------------------------------------------
#
# Pure-numpy test (no Model, no PETSc DMPlex) — instantiates the STRAMesh
# class via `__new__` to bypass the heavy init, sets minimal state, and
# exercises `deposeStrat` + `_depthPorosity` directly with mocked PETSc
# Vec / DM objects. Lives in the fast tier alongside TESTs 1-2; runs in
# well under a second on any platform.
# ---------------------------------------------------------------------------


def test_stratigraphy_deposition_and_compaction():
    """
    Protects: stratplex.deposeStrat + stratplex._depthPorosity — the
    deposition and compaction physics that update underlying stratal
    thickness and porosity when surface deposition events accumulate
    over time.

    Constructs a fictive 3-layer stratigraphy (5 nodes × 3 layers,
    each 50 m thick at phi=0.5) and runs two phases:

      Phase 1 — deposeStrat: add a 10 m deposition to the top layer.
      Asserts:
        - stratH[top] grows by exactly the deposition amount;
        - phiS[top] is set to phi0s (fresh deposit's surface porosity);
        - stratK[top] is reset to 1.0 (default erodibility multiplier);
        - lower layers are untouched.

      Phase 2 — _depthPorosity: apply Athy's-law compaction
      (phi = phi0 * exp(depth/z0)) at synthetic layer mid-depths.
      Asserts:
        - porosity decreases monotonically with depth (deeper compacts
          more);
        - layer thickness decreases (compaction shrinks layers);
        - solid-phase volume is conserved per layer to float-noise
          precision: h_new * (1 - phi_new) == h_old * (1 - phi_old).

    Silent failures prevented:
      1. deposeStrat lays down sediment but forgets to set phiS to phi0s
         on the new layer → downstream code reads NaN porosity from the
         forward-fill in `_fillZeroPorosity`.
      2. _depthPorosity inverts the depth sign convention → porosity
         INCREASES with depth, the column unrealistically inflates.
      3. Compaction's solid-phase bookkeeping (stratplex.py:268-278)
         loses or creates mass: `newH = solidPhase / (1 - phi_new)` must
         exactly preserve `h * (1 - phi)`. Any algebraic refactor that
         changes the order of operations is checked here.
      4. The bedrock-sentinel freeze (stratplex.py:283-286) is silently
         removed: bedrockLay > 0 should keep layer indices < bedrockLay
         from changing in either thickness or porosity. Tested in a
         second assertion block below.
    """
    # Late import — STRAMesh's module triggers `from gospl._fortran import
    # strataonesed`, which fails at collection if the Fortran extension
    # is not built. `importorskip` lets the test skip cleanly in that
    # environment instead of failing the whole file.
    stratplex = pytest.importorskip("gospl.sed.stratplex")
    STRAMesh = stratplex.STRAMesh

    # ---- Mocks for the PETSc-side interface ----------------------------
    # deposeStrat calls `self.dm.globalToLocal(self.tmp, self.tmpL)` then
    # reads `self.tmpL.getArray()`. We pre-populate tmpL with a fake
    # object whose getArray() returns the desired deposition vector;
    # the globalToLocal call becomes a no-op (tmpL is already correct).
    class _MockLocalVec:
        def __init__(self, arr):
            self._arr = np.asarray(arr, dtype=np.float64)
        def getArray(self):
            return self._arr

    class _MockDM:
        def globalToLocal(self, src, dst):
            pass  # tmpL pre-populated by the test

    # ---- Build the fictive STRAMesh state ------------------------------
    # __new__ bypasses __init__; STRAMesh.__init__ only sets all four
    # state arrays to None and returns, so we replace that with explicit
    # synthetic state. No PETSc / mesh allocations involved.
    s = STRAMesh.__new__(STRAMesh)
    s.lpoints   = 5
    s.stratNb   = 3
    s.stratStep = 2        # 3 layers indexed 0..2; layer 2 is the top
    s.phi0s     = 0.5      # surface porosity
    s.z0s       = 100.0    # e-folding compaction depth (m)
    s.bedrockLay = 0       # no infinite-bedrock sentinel layer
    s.memclear  = False
    s.stratLith = False    # single-fraction path (dual-lithology disabled)

    s.stratH = np.full((5, 3), 50.0, dtype=np.float64)
    s.phiS   = np.full((5, 3),  0.5, dtype=np.float64)
    s.stratK = np.full((5, 3),  0.7, dtype=np.float64)

    s.tmp  = None
    s.tmpL = _MockLocalVec(np.full(5, 10.0))   # 10 m deposition at every node
    s.dm   = _MockDM()

    # ==== Phase 1: deposition ===========================================
    H_before  = s.stratH.copy()
    phi_before = s.phiS.copy()
    K_before  = s.stratK.copy()

    s.deposeStrat()

    # Top layer thickness grew by exactly the deposition amount
    np.testing.assert_array_almost_equal(
        s.stratH[:, s.stratStep] - H_before[:, s.stratStep],
        np.full(5, 10.0),
        err_msg=(
            "deposeStrat did not add the deposition to the top layer "
            "thickness. Check stratplex.py:157 — `self.stratH[:, "
            "self.stratStep] += depo`."
        ),
    )

    # Top layer porosity reset to phi0s
    np.testing.assert_array_equal(
        s.phiS[:, s.stratStep],
        np.full(5, s.phi0s),
        err_msg=(
            "deposeStrat did not reset the top layer's porosity to "
            "phi0s after deposition (stratplex.py:159). Downstream "
            "code that forward-fills zero porosity will inherit the "
            "wrong value from below."
        ),
    )

    # Top layer K multiplier reset to 1.0
    np.testing.assert_array_equal(
        s.stratK[:, s.stratStep],
        np.ones(5),
        err_msg=(
            "deposeStrat did not reset the top layer's erodibility "
            "multiplier to 1.0 (stratplex.py:163). Freshly deposited "
            "sediment should always carry the default K."
        ),
    )

    # Lower layers (indices 0, 1) untouched
    np.testing.assert_array_equal(s.stratH[:, :s.stratStep], H_before[:, :s.stratStep])
    np.testing.assert_array_equal(s.phiS[:, :s.stratStep],   phi_before[:, :s.stratStep])
    np.testing.assert_array_equal(s.stratK[:, :s.stratStep], K_before[:, :s.stratStep])

    # ==== Phase 2: compaction ===========================================
    # Layer geometry after Phase 1 (per node):
    #   layer 0: H=50,  phi=0.5,  bottom of column
    #   layer 1: H=50,  phi=0.5
    #   layer 2: H=60,  phi=0.5,  top of column (50 + 10 m deposit)
    # Mid-point depths below the post-deposition surface (z=0):
    #   layer 2 mid-depth = -30   (60/2 below surface)
    #   layer 1 mid-depth = -85   (60 + 50/2)
    #   layer 0 mid-depth = -135  (60 + 50 + 50/2)
    # _depthPorosity expects depth as (lpoints, stratStep+1) with the
    # column ordered [bottom, ..., top] — same layout as stratH.
    depth = np.tile(np.array([-135.0, -85.0, -30.0]), (5, 1))

    H_pre   = s.stratH.copy()
    phi_pre = s.phiS.copy()

    # Snapshot solid-phase volume per layer BEFORE compaction.
    # The conservation law: h * (1 - phi) is invariant across compaction
    # (compaction only changes void volume, never solid mass).
    solid_before = H_pre * (1.0 - phi_pre)

    newH = s._depthPorosity(depth)

    # 1. Porosity strictly decreases with depth
    assert (s.phiS[:, 0] < s.phiS[:, 1]).all(), (
        "Bottom layer porosity should be less than middle layer after "
        f"compaction. Got phi[bottom]={s.phiS[:, 0]}, "
        f"phi[middle]={s.phiS[:, 1]}. Check the sign of `depth/z0s` in "
        "stratplex.py:265."
    )
    assert (s.phiS[:, 1] < s.phiS[:, 2]).all(), (
        "Middle layer porosity should be less than top layer after "
        f"compaction. Got phi[middle]={s.phiS[:, 1]}, "
        f"phi[top]={s.phiS[:, 2]}."
    )

    # 2. Layer thickness strictly decreases (compaction shrinks layers)
    assert (newH[:, 0] < H_pre[:, 0]).all(), (
        f"Bottom layer should compact. Got newH={newH[:, 0]}, "
        f"H_pre={H_pre[:, 0]}."
    )
    assert (newH[:, 1] < H_pre[:, 1]).all(), (
        f"Middle layer should compact. Got newH={newH[:, 1]}, "
        f"H_pre={H_pre[:, 1]}."
    )
    assert (newH[:, 2] < H_pre[:, 2]).all(), (
        f"Top layer should compact (its phi went from 0.5 to "
        f"{s.phiS[0, 2]:.4f}). Got newH={newH[:, 2]}, "
        f"H_pre={H_pre[:, 2]}."
    )

    # 3. Solid-phase volume conservation per layer
    solid_after = newH * (1.0 - s.phiS)
    np.testing.assert_allclose(
        solid_after, solid_before, rtol=1e-12, atol=1e-12,
        err_msg=(
            "Solid-phase volume not conserved across compaction. "
            "_depthPorosity must change phi and H in lockstep so "
            "h*(1-phi) is invariant — check the construction "
            "`newH = solidPhase / tot` at stratplex.py:278."
        ),
    )

    # ==== Phase 3: bedrock sentinel freeze ==============================
    # Re-run the compaction with bedrockLay = 1 (layer 0 is the infinite-
    # bedrock sentinel). The bedrock layer's thickness and porosity must
    # NOT change even though `depth` would otherwise drive compaction
    # there. Catches regressions of the freeze logic at stratplex.py:283.
    s2 = STRAMesh.__new__(STRAMesh)
    s2.lpoints   = 5
    s2.stratNb   = 3
    s2.stratStep = 2
    s2.phi0s     = 0.5
    s2.z0s       = 100.0
    s2.bedrockLay = 1      # layer 0 is bedrock
    s2.memclear  = False
    s2.stratLith = False   # single-fraction path (dual-lithology disabled)

    # Layer 0 holds the BEDROCK_SENTINEL thickness (1e6) — never compact it.
    s2.stratH = np.array([[1.0e6, 50.0, 50.0]] * 5, dtype=np.float64)
    s2.phiS   = np.array([[0.0,    0.5,  0.5]] * 5, dtype=np.float64)
    s2.stratK = np.full((5, 3), 1.0, dtype=np.float64)

    depth2 = np.tile(np.array([-1.0e6 - 25.0, -75.0, -25.0]), (5, 1))
    H_bedrock_before   = s2.stratH[:, 0].copy()
    phi_bedrock_before = s2.phiS[:, 0].copy()

    newH2 = s2._depthPorosity(depth2)

    np.testing.assert_array_equal(
        newH2[:, 0], H_bedrock_before,
        err_msg=(
            "Bedrock sentinel layer (index < bedrockLay) was compacted. "
            "The freeze at stratplex.py:285 should preserve `stratH` "
            "exactly on bedrock indices."
        ),
    )
    np.testing.assert_array_equal(
        s2.phiS[:, 0], phi_bedrock_before,
        err_msg=(
            "Bedrock sentinel porosity was modified. The freeze at "
            "stratplex.py:286 should preserve `phiS` exactly on "
            "bedrock indices."
        ),
    )


@pytest.mark.slow
@pytest.mark.parametrize("fixture", ["minimal_strat", "minimal_dual", "minimal_prov"])
def test_transport_limited_strata_conserve_volume(fixture, tmp_path):
    """
    Protects: STRAMesh._netRoutedSource. With stratigraphy on, the sediment
    routed downstream is built from erodeStrat's eroded thickness. With
    transport-limited SPL (`spl: G > 0`) the eroder ALSO deposits part of that
    sediment in place (deposeStrat), and the gross erosion used to be routed,
    so the in-place deposit was laid down a second time downstream: a volume
    GAIN (goSPL-examples stratigraphic_record/input-stratiG, G = 3: 38% more
    sediment into the sea than eroded). Every strata fixture has G = 0, which
    is why nothing caught it.

    On the closed sphere (no outlet, so net volume must be ~0) with G = 1:
    +9.0e-4 of the activity before the fix (+9.4e-4 dual, +9.0e-4 prov),
    -5e-5 after (the documented DEPOSIT_FLOOR / pit-residue floor).
    """
    import json
    import os
    import re
    import shutil

    from gospl.model import Model

    for f in FIXTURES_DIR.iterdir():
        if f.is_file():
            shutil.copy(f, tmp_path / f.name)
    yml = (tmp_path / f"{fixture}.yml").read_text()
    yml, n = re.subn(r"(\n\s*G:\s*)0\.", r"\g<1>1.", yml)
    assert n == 1, "fixture no longer sets spl: G: 0."
    (tmp_path / "g.yml").write_text(yml)

    cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        m = Model("g.yml", verbose=False, showlog=False,
                  summary=str(tmp_path / "s.jsonl"))
        try:
            m.runProcesses()
        finally:
            m.destroy()
    finally:
        os.chdir(cwd)

    recs = [json.loads(line) for line in open(tmp_path / "s.jsonl")]
    assert recs[0]["closed"], "minimal fixtures are a closed sphere"
    steps = recs[1:-1]
    ero = sum(r["volume"]["eroded"] for r in steps)
    dep = sum(r["volume"]["deposited"] for r in steps)
    assert dep > 0 and ero < 0
    net_rel = (ero + dep) / (dep - ero)
    assert abs(net_rel) < 2.0e-4, (
        f"{fixture} with G=1: net volume {net_rel:+.2e} of activity "
        "(in-place SPL deposit routed downstream too?)")
