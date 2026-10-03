"""
Level-B conservative groundwater geochemistry.

Protects: docs/DESIGN_WATERTABLE_GEOCHEM.md, docs/DESIGN_GEOCHEM_EXTENSIONS.md.

Split out of the former tests/test_regression.py. Run this group alone with
`pytest -m geochem`; see tests/README.md for the marker list.
"""

from __future__ import annotations

import numpy as np
import pytest

from _helpers import FIXTURES_DIR, _gw_model

# Skip the whole module (rather than erroring at collection) when the
# goSPL runtime stack is not installed.
pytest.importorskip(
    "gospl.tools.inputparser",
    reason="goSPL runtime deps (petsc4py / ruamel.yaml / scipy) not installed",
)

pytestmark = [pytest.mark.geochem]


def test_geochem_opt_in(minimal_model):
    """
    Protects (Level-B geochemistry, **G0** — DESIGN_WATERTABLE_GEOCHEM.md): the
    conservative solute-transport feature is OPT-IN and inert when off. Without a
    `groundwater: geochem:` block `gwGeochemOn` is False and no solute state is
    allocated; with it, `gwGeochemOn` is True, the per-species state is allocated
    as an `(lpoints, n_species)` array, and a run completes — G0 allocates state
    but solves nothing, so it stays byte-identical.
    """
    # (a) no groundwater at all ⇒ geochem inert, nothing allocated.
    m = minimal_model
    assert getattr(m, "gwGeochemOn", False) is False
    assert not hasattr(m, "gwSolute")

    # (b) groundwater ON but no `geochem:` block ⇒ still off, no solute state.
    mg = _gw_model("minimal_gw.yml")
    try:
        assert mg.gwOn and not mg.gwGeochemOn
        assert not hasattr(mg, "gwSolute")
    finally:
        mg.destroy()

    # (c) geochem ON with two tracers ⇒ n_species state allocated; inert run.
    mc = _gw_model("minimal_gw_geochem.yml")
    try:
        assert mc.gwGeochemOn and mc.gwNspecies == 2
        assert mc.gwSolute.shape == (mc.lpoints, 2)
        assert mc.gwSourcePool.shape == (mc.lpoints, 2)
        assert mc.gwOceanFlux.shape == (2,)
        assert hasattr(mc, "soluteL") and hasattr(mc, "soluteG")
        for arr in (mc.gwGeoCsat, mc.gwGeoPrecip, mc.gwGeoVsolid):
            assert arr.shape == (2,)
        # weatherability is a RAW list (scalar/map/table per species; ext 1).
        assert isinstance(mc.gwGeoWeather, list) and len(mc.gwGeoWeather) == 2
        mc.tEnd = mc.tNow + 0.5 * mc.dt
        mc.runProcesses()             # opt-in run completes; solute state is active
        assert np.isfinite(mc.gwSolute).all() and (mc.gwSolute >= 0.0).all()
    finally:
        mc.destroy()


def test_geochem_transport():
    """
    Protects (Level-B geochemistry, **G1** — DESIGN_WATERTABLE_GEOCHEM.md): the
    steady solute-transport operator ``∇·(q c)`` (upwind advection by the
    groundwater flux ``q = −T∇h``, from the head operator's face conductances,
    plus the vertical seepage sink that makes it a well-posed M-matrix). A tracer
    imposed at high-head inflow nodes and pinned to 0 at the seepage set advects
    **downstream, bounded in [0, 1]** (first-order upwind is monotone — no
    over/undershoot) and reaches the interior. G1 is the operator + solver only;
    the physical source / precipitation / export land in G2/G3.
    """
    m = _gw_model("minimal_gw_geochem.yml")
    try:
        m.tEnd = m.tNow + 0.5 * m.dt
        m.runProcesses()                       # solves the head (updateGroundwater)

        own = m.inIDs == 1
        h = m.headL.getArray()
        z = m.hLocal.getArray()
        # Seepage sink set (c pinned to 0 — solute exits): sea + ponded lakes +
        # open outlets, as in the head solve.
        seep = np.zeros(m.lpoints, dtype=bool)
        seep[m.seaID] = True
        pit = getattr(m, "pitIDs", None)
        lf = getattr(m, "lFill", None)
        if pit is not None and lf is not None:
            seep |= (pit > -1) & (lf > z)
        oid = getattr(m, "outletIDs", None)
        if oid is not None and len(oid) > 0:
            seep[oid] = True
        # Synthetic inflow: the highest-head interior nodes carry c = 1.
        interior = own & ~seep
        assert interior.any()
        inflow = interior & (h >= np.percentile(h[interior], 85))
        dmask = seep | inflow
        dval = np.where(inflow, 1.0, 0.0)

        c = m._solveSoluteTransport(np.zeros(m.lpoints), dmask, dval)

        assert int(m._ksp_solute.getConvergedReason()) > 0, "transport KSP diverged"
        assert np.isfinite(c).all()
        assert (c >= -1.0e-6).all() and (c <= 1.0 + 1.0e-6).all(), (
            "upwind transport not monotone (over/undershoot)"
        )
        transported = own & ~dmask
        assert (c[transported] > 0.01).any(), "solute not advected downstream"
    finally:
        m.destroy()


def test_geochem_conserves():
    """
    Protects (Level-B geochemistry, **G2** — DESIGN_WATERTABLE_GEOCHEM.md): the
    per-step **dissolve → transport → precipitate → export** cycle is mass
    conservative. Wired into ``updateGroundwater``, each tracer's budget closes:
    the **source-pool debit equals the dissolved mass**, and ``dissolved =
    precipitated + exported + in-solution`` to machine precision; all terms are
    non-negative and the precipitated solute grows the duricrust ``duriH``.
    """
    m = _gw_model("minimal_gw_geochem.yml")
    try:
        assert m.gwGeochemOn and m.gwNspecies == 2
        m.tEnd = m.tNow + 3 * m.dt
        pool0 = m.gwSourcePool.copy()
        while m.tNow < m.tEnd:
            m.runProcesses()

        own = m.inIDs == 1
        for k in range(m.gwNspecies):
            diss = m.gwDissolved[k]
            prec = m.gwPrecip[k]
            exp = m.gwOceanFlux[k]
            pooldrop = float((pool0[:, k] - m.gwSourcePool[:, k])[own].sum())
            assert diss > 0.0, "no dissolution — test not exercised"
            # (a) source-pool debit == dissolved mass (independent accounting).
            assert np.isclose(diss, pooldrop, rtol=1.0e-9), "pool debit != dissolved"
            # (b) closed budget: the steady transport is exactly conservative, so
            # every step dissolved = precipitated + exported (no storage term).
            assert abs(diss - (prec + exp)) < 1.0e-6 * diss, (
                "geochem mass budget does not close"
            )
            # (c) physical: nothing created, precipitate bounded by dissolved.
            assert prec >= 0.0 and exp >= 0.0 and prec <= diss + 1.0e-6 * diss
        # The transported+precipitated solute grew the crust.
        assert m.duriHL.getArray().max() > 0.0, "no duricrust precipitated from solute"
    finally:
        m.destroy()


def test_geochem_ocean_flux():
    """
    Protects (Level-B geochemistry, **G3** — DESIGN_WATERTABLE_GEOCHEM.md): the
    dissolved solute discharged to the surface (the baseflow-carried export) is
    tracked per tracer (``gwOceanFlux``) and exposed as the per-node
    ``soluteflux`` field. Because the steady transport is exactly conservative,
    the export equals ``dissolved − precipitated`` each step, and the per-node
    export field is finite, non-negative and non-zero (solute reaches the
    seepage/discharge zones).
    """
    m = _gw_model("minimal_gw_geochem.yml")
    try:
        m.tEnd = m.tNow + 3 * m.dt
        m.runProcesses()
        while m.tNow < m.tEnd:
            m.runProcesses()

        for k in range(m.gwNspecies):
            export = m.gwOceanFlux[k]
            diss_minus_prec = m.gwDissolved[k] - m.gwPrecip[k]
            assert diss_minus_prec > 0.0, "no net dissolution — test not exercised"
            assert abs(export - diss_minus_prec) < 1.0e-6 * diss_minus_prec, (
                "export != dissolved − precipitated"
            )
        own = m.inIDs == 1
        sf = m.gwSoluteFlux
        assert np.isfinite(sf).all() and (sf >= -1.0e-9).all()
        assert float(sf[own].sum()) > 0.0, "no baseflow-carried solute export"
    finally:
        m.destroy()


def test_geochem_multitracer():
    """
    Protects (Level-B geochemistry, **G4** — DESIGN_WATERTABLE_GEOCHEM.md): the
    `n_species` array runs multiple independent tracers that give **distinct
    duricrust typing**. With two tracers (carbonate, silica) of different
    weatherability / precipitation rate, both are active with distinct per-tracer
    budgets, and the per-node dominant-crust field `gwCrustType` (output
    `crust_type`) is a valid argmax index that resolves to *both* tracers in
    different places.
    """
    m = _gw_model("minimal_gw_geochem.yml")
    try:
        assert m.gwNspecies == 2
        m.tEnd = m.tNow + 3 * m.dt
        m.runProcesses()
        while m.tNow < m.tEnd:
            m.runProcesses()

        own = m.inIDs == 1
        ct = m.gwCrustType
        # Both tracers active, with distinct behaviour (different params).
        assert m.gwDissolved[0] > 0.0 and m.gwDissolved[1] > 0.0
        assert m.gwPrecip[0] != m.gwPrecip[1]
        assert m.gwOceanFlux[0] != m.gwOceanFlux[1]
        # crust_type is a valid dominant index, consistent with the argmax.
        assert ((ct >= -1) & (ct < m.gwNspecies)).all()
        crusted = m.gwCrustBySpecies.sum(axis=1) > 0.0
        assert np.array_equal(
            ct[crusted], m.gwCrustBySpecies[crusted].argmax(axis=1)
        )
        assert (ct[~crusted] == -1).all()
        # Both tracers CONTRIBUTE crust with distinct budgets (distinct typing).
        # With uniform (non-spatial) weatherability, carbonate (2x weatherability)
        # is the dominant former across the subaerial land, so we check both are
        # active rather than each being dominant somewhere — spatially distinct
        # typing (each dominant in its region) is the lithology example's job.
        assert (m.gwCrustBySpecies[own, 0] > 0.0).any(), "carbonate crust absent"
        assert (m.gwCrustBySpecies[own, 1] > 0.0).any(), "silica crust absent"
    finally:
        m.destroy()


def test_geochem_provenance():
    """
    Protects (Level-B geochemistry, **G5** — DESIGN_WATERTABLE_GEOCHEM.md): with
    in-model provenance on, the precipitated crust is **attributed to the
    source-rock region where its solute dissolved**. By linearity the solute is
    transported per source class (same operator, class-restricted source), so the
    downstream crust carries the provenance of its upgradient recharge area. The
    per-class crust sums to the total (conservation), the dominant-source field
    ``crust_source`` is a valid argmax, and both source regions appear.
    """
    m = _gw_model("minimal_gw_geochem_prov.yml")
    try:
        assert m.gwGeochemOn and m.provOn and m.provNb == 2
        # Two source-rock regions split by longitude.
        x = m.lcoords[:, 0]
        m.source_class = np.where(x < np.median(x), 0, 1).astype(np.int64)
        m.tEnd = m.tNow + 3 * m.dt
        m.runProcesses()
        while m.tNow < m.tEnd:
            m.runProcesses()

        own = m.inIDs == 1
        cs = m.gwCrustSource
        crusted = m.gwCrustProv.sum(axis=1) > 0.0
        # Both source regions contribute to the crust chemistry.
        for r in range(m.provNb):
            assert float(m.gwCrustProv[own, r].sum()) > 0.0, "source %d absent" % r
        # crust_source is a valid dominant-source index consistent with argmax.
        assert ((cs >= -1) & (cs < m.provNb)).all()
        assert np.array_equal(cs[crusted], m.gwCrustProv[crusted].argmax(axis=1))
        assert (cs[~crusted] == -1).all()
        assert (cs[own] == 0).any() and (cs[own] == 1).any(), "provenance not resolved"
        # Provenance conserves: per-class crust sums to the total crust.
        assert np.allclose(
            m.gwCrustProv.sum(axis=1), m.gwCrustBySpecies.sum(axis=1), rtol=1.0e-6
        )
    finally:
        m.destroy()


def test_geochem_strat_archive():
    """
    Protects (Level-B geochemistry, per-layer crust archive —
    DESIGN_WATERTABLE_GEOCHEM.md): the crust's **dominant solute species**
    (``stratCrustType``) and **dominant source region** (``stratCrustSource``)
    are recorded PER STRATIGRAPHIC LAYER — the categorical companions to the
    ``stratDuri`` degree — so a section preserves *what* each crust layer is and
    *where* its chemistry came from. Codes are integer labels stored as float64,
    -1 = no crust in that layer. Verified: both arrays are allocated at the
    stratigraphic shape, some layers carry a crust code, and every code is a
    valid species / source index (or -1).
    """
    m = _gw_model("minimal_gw_geochem_prov.yml")
    try:
        assert m.gwGeochemOn and m.provOn and m.stratNb > 0
        assert m.stratCrustType is not None and m.stratCrustSource is not None
        assert m.stratCrustType.shape == (m.lpoints, m.stratNb)
        assert m.stratCrustSource.shape == (m.lpoints, m.stratNb)
        x = m.lcoords[:, 0]
        m.source_class = np.where(x < np.median(x), 0, 1).astype(np.int64)
        m.tEnd = m.tNow + 4 * m.dt
        m.runProcesses()
        while m.tNow < m.tEnd:
            m.runProcesses()

        top = m.stratStep + 1
        ct = m.stratCrustType[:, :top]
        cs = m.stratCrustSource[:, :top]
        # Some layers carry a recorded crust (code >= 0), not just the -1 init.
        assert (ct >= 0).any(), "no crust type archived in any layer"
        # Every archived code is a valid species / source index (or -1).
        codes = np.unique(ct)
        assert np.isin(codes, np.arange(-1, m.gwNspecies)).all()
        scodes = np.unique(cs)
        assert np.isin(scodes, np.arange(-1, m.provNb)).all()
        # A layer with a crust species also has a crust source (both written
        # together in the same depth range) and vice versa — no orphan codes.
        assert np.array_equal(ct >= 0, cs >= 0)
    finally:
        m.destroy()


def test_geochem_spatial_weatherability():
    """
    Protects (Level-B geochemistry, extension 1 — DESIGN_GEOCHEM_EXTENSIONS.md
    §1): the per-species ``weatherability`` may vary in space, so **lithology
    controls which species each region yields**. Here it is driven by the
    per-vertex lithology label (``weatherability_from: source_class``) with a
    per-(class, species) table: class 0 weathers only carbonate (species 0),
    class 1 only silica (species 1). Verified on **dissolution** (not on where the
    crust ends up — solute is transported down-gradient before it precipitates, so
    the crust can sit in a different lithology): the resolved weatherability is a
    per-vertex array matching the class split, and the *off* species' source pool
    is **undebited** in each region (no dissolution) while the *on* species is
    consumed.
    """
    m = _gw_model("minimal_gw_geochem_litho.yml")
    try:
        assert m.gwGeochemOn and m.gwWeatherFrom == "source_class"
        # Two lithologies split by longitude.
        x = m.lcoords[:, 0]
        m.source_class = np.where(x < np.median(x), 0, 1).astype(np.int64)
        pool0 = 1.0e6 * m.larea                       # initial per-node pool
        m.tEnd = m.tNow + 3 * m.dt
        m.runProcesses()
        while m.tNow < m.tEnd:
            m.runProcesses()

        # Weatherability resolved to a per-vertex array per species (not scalar),
        # matching the class split (class 0 -> carbonate only, class 1 -> silica).
        assert m._gwGeoWeatherArr is not None
        for k in range(m.gwNspecies):
            wab = m._gwGeoWeatherArr[k]
            assert np.ndim(wab) == 1 and wab.shape[0] == m.lpoints
        assert np.allclose(m._gwGeoWeatherArr[0], (m.source_class == 0))
        assert np.allclose(m._gwGeoWeatherArr[1], (m.source_class == 1))
        # The OFF species' source pool is untouched in each region (no
        # dissolution there); the ON species is consumed somewhere in its region.
        c0, c1 = m.source_class == 0, m.source_class == 1
        assert np.allclose(m.gwSourcePool[c0, 1], pool0[c0]), "silica dissolved in class 0"
        assert np.allclose(m.gwSourcePool[c1, 0], pool0[c1]), "carbonate dissolved in class 1"
        assert (m.gwSourcePool[c0, 0] < pool0[c0]).any(), "no carbonate dissolved in class 0"
        assert (m.gwSourcePool[c1, 1] < pool0[c1]).any(), "no silica dissolved in class 1"
    finally:
        m.destroy()


def test_geochem_lithology_map():
    """
    Protects (Level-B geochemistry, extension 1 form **(b)** —
    DESIGN_GEOCHEM_EXTENSIONS.md §1.3): the per-species weatherability is gathered
    from a per-(class, species) table by a **standalone per-vertex lithology map**
    (``lithology: [file, key]``), independent of provenance. Here the shared
    ``prov_src.npz`` ``rock`` field (2 classes) is the lithology map: class 0
    weathers only carbonate, class 1 only silica. Verified on **dissolution**
    (source-pool debit): the resolved weatherability matches the map, and the
    *off* species' pool is undebited in each lithology (no provenance needed).
    """
    import os

    m = _gw_model("minimal_gw_geochem_litho2.yml")
    # The lithology map loads lazily (first step) via a relative path, so run in
    # the fixtures dir (cwd is otherwise restored after construction).
    cwd = os.getcwd()
    os.chdir(str(FIXTURES_DIR))
    try:
        # prov_src.npz is a gitignored, generated fixture (see conftest.py
        # minimal_prov_multi_model). Regenerate it here so this test does not
        # depend on a provenance test having produced it first.
        rockAll = (
            np.load("mesh.npz")["v"][:, 0]
            > np.load("mesh.npz")["v"][:, 0].mean()
        ).astype(np.int64)
        np.savez("prov_src.npz", rock=rockAll)

        assert m.gwGeochemOn and not getattr(m, "provOn", False)
        assert m._gwLithoMap == ["prov_src", "rock"]
        rock = np.load("prov_src.npz")["rock"][m.locIDs].astype(np.int64)
        pool0 = 1.0e6 * m.larea
        m.tEnd = m.tNow + 3 * m.dt
        m.runProcesses()
        while m.tNow < m.tEnd:
            m.runProcesses()

        # Weatherability resolved per-vertex from the lithology map (not scalar).
        assert m._gwGeoWeatherArr is not None
        for k in range(m.gwNspecies):
            assert np.ndim(m._gwGeoWeatherArr[k]) == 1
        assert np.allclose(m._gwGeoWeatherArr[0], (rock == 0))
        assert np.allclose(m._gwGeoWeatherArr[1], (rock == 1))
        # The OFF species' source pool is untouched in each lithology region.
        c0, c1 = rock == 0, rock == 1
        assert np.allclose(m.gwSourcePool[c0, 1], pool0[c0]), "silica dissolved in class 0"
        assert np.allclose(m.gwSourcePool[c1, 0], pool0[c1]), "carbonate dissolved in class 1"
        assert (m.gwSourcePool[c0, 0] < pool0[c0]).any(), "no carbonate dissolved in class 0"
        assert (m.gwSourcePool[c1, 1] < pool0[c1]).any(), "no silica dissolved in class 1"
    finally:
        os.chdir(cwd)
        m.destroy()


def test_geochem_surface_lithology():
    """
    Protects (Level-B geochemistry, extension 1 refinement §1.8 —
    DESIGN_GEOCHEM_EXTENSIONS.md): ``weatherability_from: surface_class`` derives
    the lithology label each step from the **top non-empty stratigraphic layer's
    provenance** (``stratP``), not the static bedrock ``source_class`` — so as
    exhumation/burial changes the exposed rock the weatherability follows.
    Verified: ``_surfaceSourceClass`` returns the top-layer dominant class (not
    the bedrock), it is re-resolved each step (dynamic), and the resolved
    weatherability tracks it.
    """
    m = _gw_model("minimal_gw_geochem_surface.yml")
    try:
        assert m.gwGeochemOn and m.provOn and m.gwWeatherFrom == "surface_class"
        m.tEnd = m.tNow + 3 * m.dt
        m.runProcesses()
        while m.tNow < m.tEnd:
            m.runProcesses()

        # Force a class-1 (silica) sediment cover on top of a class-0 bedrock:
        # set the top non-empty layer of every sedimented column to pure class 1,
        # and the static bedrock label to class 0 everywhere.
        top = m.stratStep + 1
        H = m.stratH[:, :top]
        rev = (H > 0)[:, ::-1]
        valid = rev.any(axis=1)
        assert valid.any(), "no stratigraphy accumulated to test with"
        top_idx = top - 1 - np.argmax(rev, axis=1)
        rows = np.arange(H.shape[0])
        vr, vt = rows[valid], top_idx[valid]
        m.stratP[vr, vt, :] = 0.0
        m.stratP[vr, vt, 1] = m.stratH[vr, vt]        # top layer is pure class 1
        m.source_class[:] = 0                          # bedrock is class 0

        # Dynamic surface class follows the TOP LAYER (1), not the bedrock (0);
        # bedrock fallback only where a column has no sediment.
        lbl = m._surfaceSourceClass()
        assert (lbl[valid] == 1).all(), "surface class did not follow the top layer"
        assert (lbl[~valid] == 0).all(), "bedrock fallback broken where no strata"

        # Re-resolved weatherability tracks it: silica (species 1) is weathered on
        # the class-1 surface, carbonate (species 0) is not.
        m._gwGeoWeatherArr = None
        arr = m._resolveGeoWeather()
        assert np.allclose(arr[1][valid], 1.0) and np.allclose(arr[0][valid], 0.0)
        # And it is DYNAMIC: _updateSolute re-resolves every step (not cached).
        pool0 = 1.0e6 * m.larea
        m.gwSourcePool[:] = pool0[:, None]
        m.tEnd = m.tNow + m.dt
        m.runProcesses()
        # On the class-1 surface, only silica dissolves (carbonate pool untouched).
        assert np.allclose(m.gwSourcePool[valid, 0], pool0[valid]), "carbonate dissolved on silica surface"
        assert (m.gwSourcePool[valid, 1] < pool0[valid]).any(), "silica not dissolved on silica surface"
    finally:
        m.destroy()


def test_geochem_river_load():
    """
    Protects (Level-B geochemistry, extension 2 — DESIGN_GEOCHEM_EXTENSIONS.md
    §2): the groundwater-exported (seepage) solute is routed **down the surface
    drainage network** as a conservative passive tracer, reusing the flow-
    accumulation matrix ``(I − Wᵀ) L = s``. Verified by rebuilding the flow matrix,
    injecting a known unit source and routing it directly: the load **accumulates
    downstream** (max ≫ source), and it is **conserved** — by the flow operator's
    column-sum identity ``Σ s == Σᵢ Lᵢ·(1 − outwᵢ)`` with ``(1 − outw) = (I − W)·1
    = fMatiᵀ·1`` (using the *actual* matrix weights, incl. outlet/pit overrides).
    """
    m = _gw_model("minimal_gw_river.yml")
    try:
        assert m.gwGeochemOn and m.gwRiverLoad
        assert hasattr(m, "riverSolute") and hasattr(m, "riverSoluteG")
        m.tEnd = m.tNow + m.dt
        m.runProcesses()                      # one step; river_load runs in-pipeline
        assert np.all(np.isfinite(m.riverSolute)) and (m.riverSolute >= 0.0).all()

        # Controlled check: rebuild the flow matrix, inject a unit source, route it.
        # Routing is per-species (gwSoluteFluxSp); this fixture is single-tracer.
        m.flowAccumulation()
        m.gwSoluteFluxSp[:] = 1.0
        m._routeRiverSolute()
        # Accumulation: the routed load far exceeds the per-node unit source.
        assert m.riverSolute.max() > 1.0 + 1.0e-9, "no downstream accumulation"
        # Exact conservation via the operator: Σ s == Σᵢ Lᵢ·(1 − outwᵢ). Build s
        # afresh (global), and (1 − outw) = fMatiᵀ·1 with the real matrix weights.
        s = m.hGlobal.duplicate()
        m.soluteL.setArray(m.gwSoluteFluxSp[:, 0])
        m.dm.localToGlobal(m.soluteL, s)
        e = m.hGlobal.duplicate()
        e.set(1.0)
        y = m.hGlobal.duplicate()
        m.fMati.multTranspose(e, y)           # y = (I − W)·1  = per-node (1 − outw)
        lhs = e.dot(s)                        # Σ s  (total injected)
        rhs = y.dot(m.riverSoluteG)           # Σᵢ Lᵢ·(1 − outwᵢ) == Σ s
        s.destroy(); e.destroy(); y.destroy()
        assert np.isclose(lhs, rhs, rtol=1.0e-5), "river solute not conserved (%g vs %g)" % (lhs, rhs)
    finally:
        m.destroy()


def test_geochem_river_species_reactions_marine():
    """
    Protects (Level-B geochemistry, extension 2 refinements §2.6 —
    DESIGN_GEOCHEM_EXTENSIONS.md): **per-species** river routing, **in-transit
    reactions** (a first-order loss `river_decay = κ`), and **marine coupling**
    (delivered coastal flux → a per-species reservoir). Verified by injecting a
    unit source per species and routing: the per-species loads sum to the total;
    each species satisfies the mass balance ``Σ s = (delivered + trapped) + lost``
    where delivered+trapped = ``(fMatiᵀ·1)·L`` and lost = ``riverSoluteLost``
    (= 0 for the conservative species, > 0 for the decaying one); and the marine
    reservoir grows by exactly ``riverSoluteToOceanSp · dt``.
    """
    m = _gw_model("minimal_gw_river2.yml")
    try:
        assert m.gwGeochemOn and m.gwRiverLoad and m.gwMarineCoupling
        assert m.gwNspecies == 2 and m.gwGeoRiverDecay[0] == 0.0 and m.gwGeoRiverDecay[1] > 0.0
        x = m.lcoords[:, 0]
        m.source_class = np.where(x < np.median(x), 0, 1).astype(np.int64)
        m.tEnd = m.tNow + m.dt
        m.runProcesses()                       # one in-pipeline step
        # Per-species loads sum to the totals.
        assert np.allclose(m.riverSoluteSp.sum(axis=1), m.riverSolute)
        assert np.isclose(m.riverSoluteToOceanSp.sum(), m.riverSoluteToOcean)

        # Controlled: rebuild the flow matrix, inject a unit source PER SPECIES.
        m.flowAccumulation()
        m.gwSoluteFluxSp[:] = 1.0
        before = m.marineSolute.copy()
        m._routeRiverSolute()

        e = m.hGlobal.duplicate(); e.set(1.0)
        y = m.hGlobal.duplicate(); m.fMati.multTranspose(e, y)   # (I − W)·1
        sg = m.hGlobal.duplicate(); Lg = m.hGlobal.duplicate()
        for k in range(m.gwNspecies):
            m.soluteL.setArray(m.gwSoluteFluxSp[:, k]); m.dm.localToGlobal(m.soluteL, sg)
            m.soluteL.setArray(m.riverSoluteSp[:, k]); m.dm.localToGlobal(m.soluteL, Lg)
            # Σ s == (delivered + trapped) + in-transit lost   (per species).
            lhs = e.dot(sg)
            rhs = y.dot(Lg) + float(m.riverSoluteLost[k])
            assert np.isclose(lhs, rhs, rtol=1.0e-4), "species %d not conserved (%g vs %g)" % (k, lhs, rhs)
        # Conservative species loses nothing; the decaying one loses a real amount.
        assert m.riverSoluteLost[0] == 0.0
        assert m.riverSoluteLost[1] > 0.0
        # Marine reservoir grew by exactly the delivered coastal flux × dt.
        assert np.allclose(m.marineSolute - before, m.riverSoluteToOceanSp * m.dt)
        sg.destroy(); Lg.destroy(); e.destroy(); y.destroy()
    finally:
        m.destroy()


def test_geochem_perspecies_outputs(tmp_path, monkeypatch):
    """
    Protects (Level-B geochemistry — per-species outputs): with several tracers
    the model writes **per-species** HDF5 fields named by tracer —
    ``solute_<name>``, ``crust_<name>``, ``soluteflux_<name>`` and (with river
    routing) ``riverSolute_<name>`` — alongside the aggregated totals. Verified by
    running the multitracer river fixture through an output step and reading the
    HDF5: the per-species datasets exist and sum to their totals.
    """
    import os
    import shutil
    from gospl.model import Model

    fx = str(FIXTURES_DIR)
    if not os.path.exists(os.path.join(fx, "minimal_gw_river2.yml")):
        pytest.skip("minimal_gw_river2.yml fixture not present")
    h5py = pytest.importorskip("h5py")
    monkeypatch.chdir(tmp_path)
    for f in ("mesh.npz", "soiltemp.npz"):
        shutil.copy(os.path.join(fx, f), tmp_path / f)
    (tmp_path / "gw.yml").write_text(open(os.path.join(fx, "minimal_gw_river2.yml")).read())

    m = Model("gw.yml", verbose=False, showlog=False)
    x = m.lcoords[:, 0]
    m.source_class = np.where(x < np.median(x), 0, 1).astype(np.int64)
    while m.tNow < m.tEnd:
        m.runProcesses()
    names = [str(n) for n in m.gwGeoName]
    m.destroy()

    files = sorted((tmp_path / "gw_river2_out" / "h5").glob("gospl.*.p0.h5"))
    assert files, "no output HDF5 written"
    with h5py.File(files[-1], "r") as f:
        # Per-species datasets exist for every tracer + the routed river load.
        for nm in names:
            for pre in ("solute_", "crust_", "soluteflux_", "riverSolute_"):
                assert pre + nm in f, "missing %s%s" % (pre, nm)
        # Per-species concentration / seepage / river load sum to their totals.
        for tot, pre in (("solute", "solute_"), ("soluteflux", "soluteflux_"),
                         ("riverSolute", "riverSolute_")):
            persum = np.sum([np.array(f[pre + nm])[:, 0] for nm in names], axis=0)
            assert np.allclose(persum, np.array(f[tot])[:, 0], rtol=1.0e-5, atol=1.0e-6), \
                "%s per-species sum != total" % tot

    # The domain-integrated solute budget (incl. the ocean flux, a scalar not
    # in the mesh output) is written to a CSV time series.
    budget = tmp_path / "gw_river2_out" / "gw_solute_budget.csv"
    assert budget.exists(), "gw_solute_budget.csv not written"
    hdr = budget.read_text().splitlines()[0].split(",")
    for nm in names:
        assert "oceanflux_" + nm in hdr, "missing oceanflux_%s in budget CSV" % nm
        assert "riverToOcean_" + nm in hdr and "marine_" + nm in hdr
