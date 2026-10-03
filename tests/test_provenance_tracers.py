"""
In-model provenance tracers (stratP / vSedP).

Protects: AGENTS.md > Analysis tools (In-model provenance tracers).

Split out of the former tests/test_regression.py. Run this group alone with
`pytest -m provenance`; see tests/README.md for the marker list.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

from _helpers import _strata_parser

# Skip the whole module (rather than erroring at collection) when the
# goSPL runtime stack is not installed.
pytest.importorskip(
    "gospl.tools.inputparser",
    reason="goSPL runtime deps (petsc4py / ruamel.yaml / scipy) not installed",
)

pytestmark = [pytest.mark.provenance]


def test_provenance_opt_in():
    """
    Protects: DESIGN_PROVENANCE.md §6 Phase 0 — in-model provenance tracers are
    opt-in via the `provenance:` block (parsed in `_extraProvenance`, a
    continuation of `_extraStrata`), and require stratigraphy.
    """
    # No block -> off.
    p = _strata_parser(stratNb=5)
    p._extraStrata()
    assert p.provOn is False and p.provNb == 0

    # On with stratigraphy -> parsed.
    p = _strata_parser(stratNb=5)
    p.input = {"provenance": {"classes": 3, "uniform": 0}}
    p._extraStrata()
    assert p.provOn is True and p.provNb == 3
    assert p._provSourceUniform == 0 and p._provSourceMap is None

    # On but stratigraphy off -> forced off.
    p = _strata_parser(stratNb=0)
    p.input = {"provenance": {"classes": 3}}
    p._extraStrata()
    assert p.provOn is False


@pytest.mark.slow
def test_provenance_seeding(minimal_prov_model):
    """
    Protects: Phase 0 provenance state — stratP is allocated (lpoints, stratNb,
    n_classes) and seeded so every initial layer carries the node's bedrock
    source class (Σ over classes == stratH). Passive in Phase 0, so the run is
    unaffected; the conservation through erosion/deposition is a later phase.
    """
    m = minimal_prov_model
    assert m.provOn and m.provNb == 2
    assert m.stratP.shape == (m.lpoints, m.stratNb, 2)
    assert (m.source_class == 1).all()                 # uniform class 1
    # Seeded: all initial thickness in class 1, none in class 0.
    assert np.allclose(m.stratP[:, :, 1], m.stratH)
    assert np.allclose(m.stratP[:, :, 0], 0.0)
    # Passive tracer in Phase 0 — the model still runs end to end.
    m.runProcesses()
    H = m.hLocal.getArray()
    assert np.isfinite(H).all()


@pytest.mark.slow
def test_provenance_erosion_split(minimal_prov_model):
    """
    Protects: Phase B1 — erodeStrat splits the eroded sediment by provenance.
    Eroded bedrock is attributed to the node's source class, the per-class
    eroded flux (provEro) sums to the total uncompacted erosion, and stratP
    stays consistent (Σ over classes == stratH).
    """
    m = minimal_prov_model              # uniform bedrock source class 1, 2 classes
    n = m.lpoints

    # Impose a uniform erosion and run only the stratigraphic erosion step.
    m.tmpL.setArray(np.full(n, -5.0))
    m.dm.localToGlobal(m.tmpL, m.tmp)
    m.erodeStrat()

    # Eroded sediment is all source class 1 (no class-0 bedrock anywhere).
    assert np.allclose(m.provEro[:, 0], 0.0)
    assert (m.provEro[:, 1] >= -1.0e-12).all()
    assert float(m.provEro[:, 1].max()) > 0.0, "no provenance eroded"
    # Per-class eroded flux sums to the total uncompacted erosion (thCoarse,
    # single-fraction here) — so the routed sub-fluxes will sum to the total.
    assert np.allclose(m.provEro.sum(axis=1), m.thCoarse)
    # stratP partitions stratH exactly after erosion.
    top = m.stratStep + 1
    assert np.allclose(m.stratP[:, :top, :].sum(axis=2), m.stratH[:, :top])


@pytest.mark.slow
def test_provenance_conservation(minimal_prov_model):
    """
    Protects: Phases B2+B3 — provenance is carried conservatively through the
    full erosion -> transport -> deposition -> stratigraphy loop. With a single
    source class, every stratigraphic layer must remain 100% that class after a
    run (a leak would put thickness in another class or break Σ == stratH).
    """
    m = minimal_prov_model              # uniform bedrock source class 1
    m.runProcesses()

    top = m.stratStep + 1
    H = m.stratH[:, :top]
    P = m.stratP[:, :top, :]
    # No spurious creation of the absent source class 0.
    assert float(np.abs(P[:, :, 0]).max()) == 0.0
    # All recorded sediment is attributed to the single source (class 1) and the
    # provenance exactly partitions the layer thickness.
    relH = max(float(H.sum()), 1.0)
    assert float(np.abs(P[:, :, 1] - H).sum()) / relH < 1.0e-6
    assert float(np.abs(P.sum(axis=2) - H).sum()) / relH < 1.0e-6


def test_provenance_pit_fraction(minimal_prov_multi_model):
    """
    Protects: Phase B2b-pit — _pitProvFraction sets a continental pit/lake
    deposit's per-node source composition to the pit's cascade-retained mix.
    The invariant is that the per-pit retained provenance (Σ over classes ==
    retained volume) yields per-node fractions that sum to 1 and reproduce the
    retained mix uniformly across the pit, so depoProvFrac stays summed-to-1 and
    stratP keeps partitioning stratH exactly.
    """
    m = minimal_prov_multi_model
    n = m.lpoints

    # Synthesise a single pit covering the first third of the local nodes with a
    # known retained mix (60% class 0, 40% class 1), then drive the method.
    m.pitIDs = np.full(n, -1, dtype=np.int64)
    in_pit = np.zeros(n, dtype=bool)
    in_pit[: max(n // 3, 1)] = True
    m.pitIDs[in_pit] = 0

    # Match pitParams length to the (one) pit so _pitProvFraction can index it.
    depo = np.array([1000.0], dtype=np.float64)         # retained volume (m^3)
    m._pitRetProv = np.array([[600.0, 400.0]], dtype=np.float64)
    m.depoProvFrac = np.zeros((n, m.provNb), dtype=np.float64)
    # Stub pitParams so len(self.pitParams) == 1.
    m.pitParams = np.zeros((1, 3), dtype=np.float64)

    m._pitProvFraction(depo)

    # Every in-pit node carries the retained mix; rows sum to 1.
    assert np.allclose(m.depoProvFrac[in_pit, 0], 0.6)
    assert np.allclose(m.depoProvFrac[in_pit, 1], 0.4)
    assert np.allclose(m.depoProvFrac[in_pit].sum(axis=1), 1.0)
    # Non-pit nodes are untouched (left at zero here).
    assert np.allclose(m.depoProvFrac[~in_pit], 0.0)


@pytest.mark.slow
def test_provenance_multisource(minimal_prov_multi_model):
    """
    Protects: multi-source provenance (2 classes from a map) — both classes are
    tracked, the per-layer composition still partitions the layer thickness
    exactly (Σ over classes == stratH), and the attribution is spatially
    sensible (each source dominates the strata over its own region).
    """
    m = minimal_prov_multi_model
    assert m.provOn and m.provNb == 2
    m.runProcesses()

    top = m.stratStep + 1
    H = m.stratH[:, :top]
    P = m.stratP[:, :top, :]
    # Conservation: provenance exactly partitions the stratigraphy.
    assert float(np.abs(P.sum(axis=2) - H).sum()) / max(float(H.sum()), 1.0) < 1.0e-6
    # Both source classes are present in the record.
    assert float(P[:, :, 0].max()) > 0.0 and float(P[:, :, 1].max()) > 0.0
    # Spatially sensible: each source's share dominates over its own region.
    src = m.source_class
    tot = P.sum(axis=2).sum(axis=1)
    has = tot > 0
    frac0 = np.divide(P[:, :, 0].sum(axis=1), tot, out=np.zeros_like(tot), where=has)
    reg0 = has & (src == 0)
    reg1 = has & (src == 1)
    if reg0.any() and reg1.any():
        assert frac0[reg0].mean() > frac0[reg1].mean()
    # B2b (marine): the recorded composition matches the eroded supply mix — the
    # marine sink (the dominant deposition here) carries the basin-delivered
    # provenance, so the deposited per-class ratio tracks the eroded ratio.
    ero = m._provEroded
    dep = m._provDeposited
    if ero.sum() > 0 and dep.sum() > 0:
        assert abs(dep[1] / dep.sum() - ero[1] / ero.sum()) < 1.0e-2

    # B2b marine diffusion lockstep: the per-class deposit thickness is spread
    # by the SAME marine diffusion as the total (_diffuseProvTracers), so the
    # diffused far field carries a routed composition instead of the domain
    # average. Confirm the path ran and the post-diffusion composition partitions
    # the deposit (Σ_c == 1, never > 1) on every marine deposit node.
    md = getattr(m, "_marDiffProv", None)
    assert md is not None and md.shape == (m.lpoints, m.provNb)
    m.dm.globalToLocal(m.tmp, m.tmpL)
    denom = md.sum(axis=1)
    routed = denom > 0.0
    if routed.any():
        comp = md[routed] / denom[routed, None]
        assert comp.max() <= 1.0 + 1e-9
        assert np.allclose(comp.sum(axis=1), 1.0)


def test_provenance_hillslope_routing(minimal_prov_multi_model):
    """
    Protects: _hillslopeProvFraction threads provenance through the hillslope
    diffusion — a creep deposit gets the flux-weighted eroded composition of its
    higher (donor) neighbours, not the bedrock fallback. Mesh-independent
    invariant: if ALL eroded material is one class, every fed deposition node's
    routed composition must be exactly that class (pure), and valid (≤1) /
    never > 1 elsewhere.
    """
    m = minimal_prov_multi_model
    n, C = m.lpoints, m.provNb

    # Tilted pre-diffusion surface (increases with x): every interior node has a
    # higher (donor) neighbour towards +x.
    x = m.lcoords[:, 0]
    zOld = (x - float(x.min())).astype(np.float64)
    m.hOld.setArray(zOld)
    # Post-diffusion surface: deposit (dz > 0) on the lower-x half.
    znew = zOld.copy()
    lo = x < float(np.median(x))
    znew[lo] += 0.5
    m.hLocal.setArray(znew)

    # All eroded material is class 0.
    m.provEro = np.zeros((n, C), dtype=np.float64)
    m.provEro[:, 0] = 1.0
    m.depoProvFrac = np.zeros((n, C), dtype=np.float64)

    m._hillslopeProvFraction()

    fs = m.depoProvFrac.sum(axis=1)
    routed = fs > 1.0e-6
    assert routed.any(), "no creep deposit received a routed composition"
    # Every routed node is pure class 0 (all donors are class 0).
    assert np.allclose(m.depoProvFrac[routed, 0], 1.0)
    assert np.allclose(m.depoProvFrac[routed, 1:], 0.0)
    # Never over-filled.
    assert m.depoProvFrac.max() <= 1.0 + 1.0e-9


@pytest.mark.slow
def test_provenance_output_io(minimal_prov_model):
    """
    Protects: Phase B4 (I/O) — the per-layer provenance composition stratP is
    written to the stratal HDF5 (lpoints, layers, classes), consistent with the
    recorded stratH (Σ over classes == stratH).
    """
    import glob

    model = minimal_prov_model
    model.runProcesses()
    files = sorted(
        glob.glob(os.path.join(str(model.outputDir), "h5", "stratal.*.p*.h5"))
    )
    if not files:
        pytest.skip("no stratal output written")
    h5py = pytest.importorskip("h5py")
    with h5py.File(files[-1], "r") as hf:
        assert "stratP" in hf, "provenance stratP not in stratal output"
        P = np.array(hf["stratP"])
        H = np.array(hf["stratH"])
        assert P.shape[2] == model.provNb
        relH = max(float(H.sum()), 1.0)
        assert float(np.abs(P.sum(axis=2) - H).sum()) / relH < 1.0e-6
        assert float(np.abs(P[:, :, 0]).max()) == 0.0      # single source -> class 1


def test_provenance_deposit_no_holes():
    """
    Protects: deposeStrat — a depositing node with no arriving composition
    (Σ_c depoProvFrac ≈ 0, e.g. off-channel hillslope creep with no river
    through-flux) must NOT produce a layer with thickness but zero provenance
    (Σ_c stratP < stratH), which surfaces in post-processing as a cell with no
    source class (`dominant == -1`). Such locally-derived deposits fall back to
    the in-situ bedrock `source_class`; routed nodes are renormalised to sum to
    1 exactly. Every depositing node ends with Σ_c stratP == stratH.
    """
    stratplex = pytest.importorskip("gospl.sed.stratplex")
    n, C = 3, 2

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

    m = stratplex.STRAMesh.__new__(stratplex.STRAMesh)
    m.lpoints = n
    m.stratStep = 1
    m.stratLith = False
    m.provOn = True
    m.provNb = C
    m.memclear = False
    m.phi0s = 0.49
    m.larea = np.ones(n)
    m.inIDs = np.ones(n, dtype=int)
    m._provDeposited = np.zeros(C)
    m.source_class = np.array([1, 0, 0])          # node 0 bedrock = class 1
    m.dm = _StubDM()
    depo = np.array([10.0, 5.0, 0.0])             # nodes 0,1 deposit; 2 none
    m.tmp = _StubVec(depo)
    m.tmpL = _StubVec(np.zeros(n))
    m.stratH = np.zeros((n, 2))
    m.phiS = np.zeros((n, 2))
    m.stratK = np.zeros((n, 2))
    m.stratP = np.zeros((n, 2, C))
    # node 0: hole (no composition); node 1: routed mix; node 2: no deposit
    m.depoProvFrac = np.array([[0.0, 0.0], [0.3, 0.7], [0.0, 0.0]])

    m.deposeStrat()

    P = m.stratP[:, 1, :]
    H = m.stratH[:, 1]
    # Hole node attributed to its bedrock source_class (1), fully.
    assert np.allclose(P[0], [0.0, 10.0])
    # Routed node keeps its mix.
    assert np.allclose(P[1], [1.5, 3.5])
    # Every depositing node: Σ_c stratP == stratH (no holes, no overshoot).
    dep = H > 0
    assert np.allclose(P[dep].sum(axis=1), H[dep])
    # Post-processing dominant source is defined (never -1) where deposited.
    assert (P[dep].sum(axis=1) > 0).all()


def test_provenance_compaction_rescale():
    """
    Protects: provenance B4 — getCompaction shrinks stratH (pore-water expelled)
    but compaction is composition-neutral, so stratP MUST be rescaled by the same
    per-layer ratio. Without the rescale stratP keeps the pre-compaction (larger)
    thickness while stratH shrinks, so the post-processed fraction stratP[c]/stratH
    exceeds 1 (the over-represented bedrock source). The minimal-run conservation
    tests miss this: burial depth ~0 there makes compaction a near no-op.

    Three thick, deeply-buried layers with a non-uniform 2-class mix. After
    compaction: stratH shrinks, Σ_c stratP == stratH stays exact, every per-class
    fraction stays <= 1, and the composition is unchanged.
    """
    stratplex = pytest.importorskip("gospl.sed.stratplex")
    n, L, C = 2, 3, 2

    class _StubVec:
        def __init__(self, arr):
            self._a = np.asarray(arr, dtype=np.float64)
        def getArray(self):
            return self._a
        def setArray(self, a):
            self._a = np.asarray(a, dtype=np.float64)

    class _StubDM:
        def localToGlobal(self, l, g):
            g.setArray(l.getArray().copy())

    m = stratplex.STRAMesh.__new__(stratplex.STRAMesh)
    m.lpoints = n
    m.stratStep = L - 1
    m.stratLith = False
    m.provOn = True
    m.provNb = C
    m.bedrockLay = 0
    m.memclear = False
    m.verbose = False
    m.phi0s = 0.5
    m.z0s = 2000.0
    m.dm = _StubDM()
    m.hLocal = _StubVec(np.zeros(n))          # surface at z = 0
    m.hGlobal = _StubVec(np.zeros(n))

    # Thick layers (2000 m each) => deep burial => strong compaction.
    m.stratH = np.full((n, L), 2000.0)
    m.phiS = np.full((n, L), m.phi0s)         # freshly-deposited porosity
    # Non-uniform composition, exactly partitioning each layer's thickness.
    frac0 = np.array([[0.2, 0.5, 0.9], [0.7, 0.1, 0.4]])   # (n, L) class-0 share
    P = np.zeros((n, L, C))
    P[:, :, 0] = m.stratH * frac0
    P[:, :, 1] = m.stratH * (1.0 - frac0)
    m.stratP = P.copy()

    H_before = m.stratH.copy()
    m.getCompaction()

    # Compaction actually shrank the pile (otherwise the test is vacuous).
    assert (m.stratH < H_before - 1.0).any(), "compaction did not bite"

    # Σ over classes == compacted stratH (exact).
    psum = m.stratP.sum(axis=2)
    assert np.allclose(psum, m.stratH, rtol=1e-12, atol=1e-9)

    # No per-class fraction exceeds 1 (the reported symptom).
    frac = np.divide(m.stratP, m.stratH[:, :, None],
                     out=np.zeros_like(m.stratP), where=m.stratH[:, :, None] > 0)
    assert frac.max() <= 1.0 + 1e-9, f"fraction > 1 after compaction: {frac.max()}"

    # Composition preserved (class-0 share unchanged by pore-water loss).
    assert np.allclose(frac[:, :, 0], frac0, rtol=1e-9, atol=1e-9)
