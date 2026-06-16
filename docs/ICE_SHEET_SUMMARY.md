# goSPL ice-sheet model — change summary

goSPL's glacial capability has been replaced by a true **Shallow-Ice-
Approximation (SIA) ice-sheet model**: an explicit, mass-conserving non-linear
diffusion of the ice thickness that drives glacial abrasion, till transport and
flexural loading. It supersedes the previous MFD flow-routing proxy, which is
removed — when an `ice` section is present, goSPL now runs the SIA model.

> **Note (2026-06-16).** The thickness solve was changed from implicit to
> **explicit with a mass-conserving flux limiter** (`ice_flux_limiter` +
> `ice_flux_rscaled`): the
> implicit `SNES` diverged on the `H≥0` free boundary (an obstacle problem),
> and the explicit CFL turned out benign at goSPL scales. See
> `docs/DESIGN_ICE_SHEET.md` §9. The descriptions of the *solver* below are
> historical; the dynamics, abrasion, till and loading are unchanged.

The work was delivered across the design doc and PRs #423–#425 (dynamics +
basal velocity, abrasion + till, flexural loading + output), then refined:

- **#426** — the MFD proxy and the explicit SIA reference scheme removed, leaving
  the implicit SIA as the single model; glacial till coupled to the
  stratigraphic / dual-lithology fractions; `iceMelt` / `iceAbr` outputs added;
  per-vertex / time-series ELA maps (the `glaciers` series) for global runs;
  optional `hinit` seed-and-evolve initial ice.
- **#427** — terminus floor = `max(hterm, sea level)`, defaulting to the
  sea-level position.
- **#428** — opt-in catchment-aware till routing (`till.route`) for
  high-resolution regional runs.
- **#429** — `scripts/ela_from_temperature.py` helper deriving ELA maps from
  paleo-climate temperature.

Mapped to the roadmap item *"glacial and ice-sheet model and associated
erosion / deposition and loading"*:

| Roadmap sub-item | Status |
|---|---|
| Ice dynamics (real flow, not a routing proxy) | ✅ implicit SIA |
| Glacial erosion | ✅ velocity-based abrasion |
| Glacial deposition | ✅ till → moraine (bulk + stratigraphic / dual lithology) |
| Loading | ⚠️ existing flexure leveraged, **no new flexure code** (see §4) |

Full design and rationale: `docs/DESIGN_ICE_SHEET.md`. User/technical/API docs:
`docs/user_guide/surfproc.rst`, `docs/tech_guide/ice.rst`,
`docs/api_ref/ice_ref.rst`.

---

## 1. Ice dynamics — implicit SIA

The ice thickness `H` evolves as a non-linear diffusion of the ice surface
`s = z_bed + H`:

> `∂H/∂t = ṁ − ∇·q`, with `q = −[2A(ρ_i g)ⁿ/(n+2)·H^(n+2) + a_s·H^n]·|∇s|^(n−1)·∇s`

combining Glen's-law internal **deformation** (`Aglen`, exponent `glen`) and
basal **sliding** (`slide`). The flux divergence is evaluated on the
unstructured finite-volume mesh by a dormant-but-now-wired `ice_flux` Fortran
kernel, using the same mass-conservative edge fluxes as the hillslope diffusion
operator (so it runs unchanged on the spherical/planar Voronoi mesh and fills /
overflows closed basins because it acts on the total bed + ice surface).

Because the SIA diffusivity scales as `H^(n+2)|∇s|^(n−1)`, an explicit scheme
would need a CFL step orders of magnitude below goSPL's km / 10²–10⁴ yr
resolution. The thickness is therefore advanced **implicitly**, with a cached
Jacobian-free PETSc `SNES` (`ngmres` + CG / HYPRE BoomerAMG) on the
backward-Euler residual `F(H) = H − H_old − Δt(ṁ − ∇·q(H))` — the same
non-linear-diffusion machinery as the non-linear hillslope solver. It is
unconditionally stable and takes the **full goSPL time step in one solve**. The
free boundary `H ≥ 0` is a post-clamp; ice is removed below the terminus
`hterm`.

The **surface mass balance** keeps the ELA proxy: `ṁ = P·min(1,(η−hela)/(hice−hela))`
(accumulation above the ELA, ablation below). The **basal sliding speed** `u_b`
is derived from the converged thickness (`ice_velocity` kernel) and exported.

## 2. Glacial erosion — velocity-based abrasion

Sliding ice abrades the bed at `E_g = K_g·|u_b|^l` (`abrasion.Kg`, `abrasion.l`;
off by default). It is masked to subaerial, ice-covered cells. With till
handling off, `E_g` is an incision added to the erosion–deposition rate and
flows into the fluvial sediment system through the standard erosion path in all
three SPL flavours (`SPL`, `nlSPL`, `soilSPL`).

This **replaces** the previous stream-power glacial term `Ki·Fⁱᶜᵉ^m`, which is
removed along with the MFD ice-flow field it consumed.

## 3. Glacial deposition — till and moraine

With `till.on`, abraded rock is carried as **till** by the ice and deposited
where the ice melts out — the ablation zone / terminal moraine — weighted by the
local meltwater rate. Two conservation-guarded modes:

- **Bulk bed** (no stratigraphy): the abraded volume is redistributed as a
  bed-to-bed transport, so the net bed-volume change is exactly zero (rock
  moved, not created).
- **Stratigraphic / dual lithology**: the till is removed from the stratigraphic
  layers it was abraded from (`erodeStrat`) and re-deposited as a fresh moraine
  layer in the ablation zone (`deposeStrat`), split into coarse/fine carrying the
  abraded (ice-mixed) fine fraction. Conservation is on the *solid* phase per
  fraction — the fine deposited equals the fine eroded, keeping the
  dual-lithology budget balanced — and the bed bulks up by the
  uncompacted-till vs compacted-rock porosity contrast.

**Meltwater** generated by ablation is captured (`m = max(−ṁ,0)·area·[H>0]`) and
re-injected into the river flow accumulation, so cells downstream of glacier
termini see the corresponding discharge instead of losing it.

## 4. Loading — existing flexure leveraged

No flexure source was modified. goSPL's **existing** flexural isostasy (FD/FFT
gFlex for planar models, spherical-harmonic shell for global) now receives the
SIA ice thickness as the ice load: the change in thickness between steps is
converted to an equivalent load (scaled by `ρ_i/ρ_c`) and applied through the
existing path, so a growing ice sheet drives isostatic subsidence and
deglaciation drives rebound. There is no new flexure mechanism — see the
appendix of `DUAL_LITHOLOGY_SUMMARY.md` for the flexure-method overview.

---

## Outputs

| Field | File / attribute | Meaning |
|---|---|---|
| `iceH` | mesh HDF5 / XDMF | ice thickness (m); restored on restart |
| `iceUb` | mesh HDF5 / XDMF | basal sliding speed (m/yr); abrasion driver and dynamics diagnostic |
| `iceMelt` | mesh HDF5 / XDMF | ablation meltwater (m³/yr) re-injected into the rivers — glacial discharge contribution |
| `iceAbr` | mesh HDF5 / XDMF | glacial abrasion rate `E_g = Kg·|u_b|^l` (m/yr); zero when abrasion is off |

The old `iceFA` (MFD ice flow accumulation) output is removed.

## Configuration

The `ice` section is the opt-in. SIA is the only model, so the previous
`flow_model` / `solver` selectors and the MFD-only keys (`Ki`, `icedir`, `melt`,
`diff`, `fwidth`, `eheight`) are gone:

```yaml
ice:
    hela: 1850.0          # equilibrium-line altitude (m)
    hice: 2100.0          # ice-cap altitude (m)
    hterm: 1700.0         # glacier terminus (m); floor = max(hterm, sea level),
                          #   defaults to sea level if omitted  [or evol: <csv>]
    sia:
        Aglen: 1.0e-16    # Glen rate factor
        slide: 1.0e-3     # basal sliding coefficient
        glen:  3.0        # Glen exponent n
    abrasion:
        Kg: 1.0e-4        # abrasion coefficient (0 = off)
        l:  1.0           # sliding-velocity exponent
    till:
        on: True          # carry abraded rock as till / moraine
        route: False      # True = route to termini (high-res); False = melt-spread
```

`hela`/`hice`/`hterm` accept a uniform scalar (as above), a per-vertex map
`[file, key]`, or a `glaciers` time series of either — see the user guide for
the spatial/time-varying ELA syntax used by global models. An optional `hinit`
(scalar or map) seeds a pre-existing ice thickness on a fresh start, which the
SIA solve then evolves (the ELA geometry still drives the evolution); a restart
restores the evolved thickness instead.

## Validation

Seven regression tests guard the model (`tests/test_regression.py`,
`tests/fixtures/minimal_ice_*.yml`):

| Test | Protects |
|---|---|
| `test_ice_opt_in` | parsing: ice off without the section, SIA params with it |
| `test_ice_sia_runs_and_invariants` | end-to-end run; H ≥ 0, finite, no ice below terminus/ELA; valid basal velocity |
| `test_ice_sia_volume_conservation` | SIA flux conserves ice volume under zero mass balance (Halfar-style) |
| `test_ice_glacial_abrasion` | `E_g = −K_g·u_b` exactly where ice slides; no-op when `Kg=0` |
| `test_ice_glacial_till_conserves` | bulk till: eroded == deposited, net bed-volume change ≈ 0 |
| `test_ice_glacial_till_dual_lithology` | stratigraphic till: per-fraction solid conserved; moraine carries fine |
| `test_ice_sia_flexure_loading` | SIA ice load drives isostatic subsidence; `iceUb` written |

---

## Scope notes

- **ELA mass balance retained.** Accumulation/ablation is the elevation-ramp ELA
  proxy, not an energy-balance or PDD melt model. It is the established goSPL
  forcing and is adequate for the long-term, landscape-to-global mass budget
  the model targets; a climate-driven SMB would be the upgrade for studies where
  the melt physics itself is the object of interest.
- **Spatially / temporally varying ELA.** For global runs, `hela`, `hice` and
  `hterm` can each be a per-vertex map and/or a `glaciers` time series (like the
  precipitation `climate` block), so the ELA varies with latitude (tropical vs
  polar) and through time — a single global scalar cannot represent both. The
  mass-balance ramp is evaluated per node with its local ELA. Uniform scalars
  and the `evol` CSV remain unchanged. The `scripts/ela_from_temperature.py`
  helper derives `hela`/`hice` maps from a paleo-climate temperature history by
  lapse-rate inversion (it sets the ELA *position*; ablation stays
  precipitation-scaled).
- **Till routing.** By default till is deposited across the ablation zone
  weighted by the meltwater rate — matched to coarse resolution where a cell
  aggregates many glacial elements. For high-resolution (sub-km) regional runs,
  `till.route: True` instead routes the till down the ice-surface flow network
  and melts it out toward each catchment's terminus, building moraine at the
  actual ice margins. Both conserve mass.
- **Scalar SIA coefficients.** `Aglen` / `slide` are uniform (no spatial
  thermomechanical coupling). Spatially-variable softness/sliding would matter
  only for studies resolving thermal regime or till rheology, below the scope
  here.
