# goSPL `dev` vs `master` (v2026.6.13) — change summary

All differences between `dev` and `master` are the new **dual-lithology
(coarse/fine) sediment** capability, delivered across 7 PRs (#415–#421,
27 commits, ~2,100 lines across 23 files). The feature is **opt-in** via
`strata: dual: True`; with it off, behaviour is **byte-identical** to
v2026.6.13 (enforced by a bitwise parity test).

Mapped to the roadmap items:

| Roadmap item | Status |
|---|---|
| Expand LEM output parameters | ✅ done |
| Sediment thickness — porosity | ✅ done |
| Sediment thickness — lithology | ✅ done |
| Sediment thickness — flexure | ⚠️ existing flexure leveraged, **no new flexure code** (see §3) |

Full design and rationale: `docs/DESIGN_DUAL_LITHOLOGY.md`.

---

## 1. Expand LEM output parameters

The coarse/fine split is now written wherever sediment is output:

| Field | File / attribute | Meaning |
|---|---|---|
| `sedLoadF` | mesh HDF5 / XDMF (`SLf`) | fine sediment volume flux (m³/yr). `sedLoad` stays the **total**, so **coarse = `sedLoad − sedLoadF`** |
| `stratHf` | stratal HDF5 | fine-fraction bulk thickness per layer (coarse = `stratH − stratHf`) |
| `phiF` | stratal HDF5 | per-layer fine porosity (beside the coarse `phiS`) |

All three are restart-round-tripped; gated on `stratLith` so single-fraction
outputs are unchanged. Documented in `HOW_TO_ADD_OUTPUT.md`.

## 2. Sediment thickness — porosity

Porosity/compaction is now **lithology-dependent**, the first-order control on
*preserved* (compacted) thickness:

- Two independent depth-porosity curves — coarse `{phi0c, z0c}` (defaults to the
  existing `compaction.phis/z0s`) and fine `{phi0f, z0f}`; fines start more
  porous and compact more.
- `_depthPorosityDual` compacts each fraction on its own Athy-law curve and sums
  the recompacted bulk thicknesses, conserving each fraction's solid volume.
- **Net effect:** the same solid load yields a different preserved thickness
  depending on the local sand/mud mix.

## 3. Sediment thickness — flexure & lithology

### Lithology — fully implemented

**Sediment is fully partitioned into coarse and fine, and that composition is
carried and recorded everywhere** — eroded, transported, deposited, compacted
and stored in the stratigraphic column, with per-fraction mass conserved. The
chain: per-layer input composition (`npstrata` `strataHf`/`phiF`), erosion split
by source composition, composition-weighted erodibility (`fine_k_factor`) and
diffusivity (`fine_diff_factor`), a separately routed fine flux (`vSedF`),
fine-enriched overspill (coarse trapped in filled depressions, fine carried to
the distal basin), depocenter/distal-biased lake and marine deposition,
per-fraction compaction, and fine-pile advection. Total **and** per-fraction
mass are conserved on a closed sphere.

The **glacial till** produced by the SIA ice-sheet model (now on `dev`) also
feeds this chain: when stratigraphy is on, abraded rock is removed from the
layers it came from and re-deposited as a moraine split into the coarse/fine
fractions (carrying the abraded fine fraction), so the per-fraction solid budget
stays balanced. See `docs/ICE_SHEET_SUMMARY.md`.

### Flexure — existing machinery, no new code this cycle

No flexure source was modified (`addprocess.py` flexure routines,
`mesher/tectonics.py` are untouched). goSPL's **existing** flexural isostasy
(overview in the appendix) now receives a sediment load / surface elevation that
is lithology-aware and compaction-corrected (via §2), so accommodation responds
to a more realistic deposit — but there is no new flexure mechanism.

**Flexure–lithology coupling (scope choice):** the flexural load uses a single
fill density `flex_rhos` (`rhoc`) for all sediment. Making it
composition-dependent (coarse vs fine bulk density, which differ through
compaction) would let a mud-rich depocenter load the plate differently from a
sand-rich one — a second-order refinement on the load budget that is not
required at the long-wavelength scale flexure operates on here, where the
deposit's total mass dominates the response.

---

## Scope note — composition-resolved, not geometry-resolved deposition

This is **not a missing feature**; it is a modelling choice matched to goSPL's
spatial and temporal scope. Sediment partitioning fully happens (see §3): coarse
and fine are split, transported, deposited, compacted, recorded and conserved.
What goSPL deliberately does **not** do is treat coarse and fine as two
independently transported fields that build *geometrically distinct deposit
bodies*.

**How the current model treats it.** Erodibility and diffusivity are
grain-size-aware (they affect geometry). The depositional *geometry* — where
sediment piles, the elevation change at each node — is computed once for the
**total** sediment; the coarse/fine composition is then distributed over that
geometry (fine to the depocenter/distal water, coarse proximal, plus
fine-enriched overspill). So the deposit *shape* equals a single-fraction run,
while the *composition* over it is grain-size resolved and conserved.

**What a two-field differential-transport solve would add** (each fraction
transported independently, the deposit geometry emerging from their sum):

1. fine building its own distal bodies (mud lobes / basin-floor drapes) coarse
   never reaches — topography the total-sediment solve cannot produce;
2. feedback on flow and accommodation, because the *elevation* would differ
   (the composition-only re-partition leaves elevation identical to the total
   run, so these feedbacks are absent);
3. emergent sorting (continuous downstream fining, shoreline segregation,
   hyper-/hypopycnal partitioning) rather than the prescribed depth-bias;
4. resolved bypass/competition for accommodation between the two fractions;
5. interfingering sand/mud bodies and clinoform foreset/bottomset architecture.

**Why this is not required for goSPL's case of interest.** goSPL targets
landscape-to-global domains over 10⁴–10⁸ yr, at mesh resolutions where a single
cell aggregates many channel and lobe elements. Resolving *where each grain size
piles within a delta or lobe* — and the autogenic, sub-grid feedbacks that drive
it — is a basin-/reservoir-scale, higher-resolution, shorter-timescale problem
(the province of dedicated stratigraphic forward models). At goSPL's
source-to-sink, long-term mass-and-composition budget scale, tracking the
composition of what is eroded, transported and deposited — which the model does,
conservatively — is the signal of interest; a two-field solve would resolve
detail below the model's effective resolution at a substantial cost (roughly
doubling the transport/deposition machinery, harder MPI conservation, nonlinear
coupling of the two fields). It would become relevant only *outside* this scope —
e.g. strongly grain-size-partitioned, higher-resolution basin studies
(mud-dominated deltas/shelves, turbidite lobes) where the deposit geometry per
grain size and its feedback are the object of study.

---

## Appendix — flexure approaches already in goSPL

Flexural isostasy lives in `GridProcess` (`gospl/tools/addprocess.py`), driven
by `Model.applyFlexure` each flexure step. It corrects topography for the
flexural rebound/subsidence of erosional unloading / depositional loading by
solving the thin-elastic-plate biharmonic equation

> `D ∇⁴w + Δρ g w = q`,  with `D = E·Te³ / (12(1−ν²))`, load `q = ρ_s g Δz`

(`D` flexural rigidity, `w` deflection, `Δρ = ρ_m − ρ_fill`). The YAML
`flexure.method` selects one of two families:

### Planar / Cartesian — `method: fem` (flat models)

- Parallel mixed finite-volume biharmonic solve **directly on the DMPlex**
  (`_cmptFlexFEM`): the FV negative-Laplacian `Lm` applied twice gives the
  single-field system `[Lm·diag(D)·Lm + Δρg·I] w = q`. No gather-to-root, no
  regular grid, no external dependency. (This replaced the former gFlex `FD` and
  FFT solvers.)
- Spatially-variable elastic thickness `Te` via a `temap` time series
  (`_updateTe`, per-node) or a uniform `thick` — a single linear solve either
  way (the rigidity goes into `diag(D)`; no iteration over the contrast).
- The operator + factorisation are cached and reused each step (serial PETSc LU
  / parallel MUMPS); only the RHS changes.
- Per-edge boundary conditions `bcN/S/E/W`: `0Slope0Shear` and `Mirror` (the
  natural zero-flux FV boundary) and `0Displacement0Slope` (clamped, `w=0`).

### Global / spherical — `method: global`

- Uses **pyshtools** spherical-harmonic transforms (`_cmptFlexGlobal`).
- Solves thin-elastic-**shell** flexure on a Driscoll–Healy (DH2) grid; the mesh
  ↔ DH interpolation weights are precomputed once (`_buildDHGrid`,
  inverse-distance forward + bilinear backward).
- Per-degree spectral filter `w_lm = q_lm / (Δρ g + D·P_l)`; degrees 0 and 1
  (mean & centre-of-mass drift) are dropped.
- Constant `Te` → a single spectral solve; varying `Te` → an iterative
  Picard/Anderson scheme (reference rigidity from max-`Te` for contraction).

### Common parameters (`flexure` YAML block)

`method`, `rhoa` (mantle density, default 3300), `rhoc` (load/fill density →
`flex_rhos`, default 2300), `young` (Young's modulus, 65 GPa), `nu` (Poisson,
0.25), `thick` (uniform `Te`, m), `regdx` (Cartesian grid spacing),
`res_deg` (DH grid resolution), `bcN/S/E/W`. Ice load is converted to an
equivalent sediment thickness before loading. The deflection is accumulated into
`localFlex` (an output field).
