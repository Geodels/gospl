# goSPL ice-sheet model — change summary

goSPL's glacial capability is a cheap, robust **diagnostic glacial-erosion
model**. When an `ice` section is present, each goSPL step routes the ELA
accumulation downhill into an ice discharge, from which a Bahr ice thickness and
a bounded Glen-sliding velocity are derived (one linear solve, no ice-thickness
PDE time integration); these drive glacial abrasion, till/moraine deposition,
discharge-conserving meltwater and flexural ice loading.

> **Why a diagnostic, not full ice dynamics.** A true Shallow-Ice-Approximation
> (SIA) thickness solve was implemented and then removed. The ice margin `H ≥ 0`
> is a free-boundary (obstacle) problem on which the implicit `F(H)=0` thickness
> solve diverged on real continental runs, and an ice-dynamics solve is stiff and
> **over-thickens km-scale continental ice at coarse resolution** over goSPL's
> long (10²–10⁴ yr) timesteps. The diagnostic is physical and robust at any
> resolution and over those long steps, so it is the only model goSPL ships. See
> `docs/DESIGN_ICE_SHEET.md`.

Status: **implemented / shipped**. Full design and rationale:
`docs/DESIGN_ICE_SHEET.md`. User/technical/API docs:
`docs/user_guide/surfproc.rst`, `docs/tech_guide/ice.rst`,
`docs/api_ref/ice_ref.rst`.

Mapped to the roadmap item *"glacial and ice-sheet model and associated
erosion / deposition and loading"*:

| Roadmap sub-item | Status |
|---|---|
| Ice flow (routed discharge → Bahr thickness → Glen sliding velocity) | ✅ diagnostic |
| Glacial erosion | ✅ velocity-based abrasion (vertical + optional lateral / U-shaping) |
| Glacial deposition | ✅ till → moraine (bulk + stratigraphic / dual lithology), volume-conserving |
| Loading | ⚠️ existing flexure leveraged, **no new flexure code** (see §4) |

---

## 1. Ice flow — routed discharge, Bahr thickness, Glen sliding velocity

Each step, in `gospl/flow/iceplex.py`:

1. The ELA surface mass balance `mdot` is computed (accumulation above the ELA
   `hela`, full accumulation at the ice-cap altitude `hice`, ablation below),
   with the accumulation scaled by `accum_factor` / capped by `accum_max`.
2. The **net mass balance** — accumulation minus the `melt`-scaled ablation
   `(mdot⁺ − melt·mdot⁻)` — is **routed downhill** on a terminus-anchored,
   drainage-conditioned bed (filled in parallel, no serial gather/epsfill) by a
   multiple-flow-direction (MFD) algorithm (`icedir` directions) — the same
   flow-matrix / KSP machinery as the river flow accumulation — into an **ice
   discharge** `Q` (m³/yr). The tongue ends where downstream ablation eats the
   upstream accumulation, so `melt` (default 0 = accumulation-only, 1 = true net,
   >1 = shorter tongues) controls glacier extent; the raw ablation still drives
   the till melt-out and meltwater. One linear solve, no time integration.
3. **Ice thickness** from a Bahr discharge scaling: `H = eheight·fwidth·Q^0.3`.
4. **Basal sliding velocity** `u_b ∝ H^(n-1)|∇s|^(n-1)∇s` from Glen's sliding law
   on that thickness and the bed-surface slope (`slide`, `glen`; the
   `ice_velocity` Fortran kernel) — physically bounded, exported as `iceUb`.

## 2. Glacial erosion — velocity-based abrasion

Sliding ice abrades the bed (deepening valleys) at `E_g = Kg·|u_b|^l`
(`abrasion.Kg`, `abrasion.l`; off by default). It is masked to subaerial,
ice-covered cells. An optional **lateral (valley-wall) erosion** term
`Kl·u_b_neighbour^lat_l` (`abrasion.Kl`, `abrasion.lat_l`; off by default) erodes
the walls flanking fast ice, **widening valleys toward a U-profile**; the eroded
wall rock joins the same conserved till budget. With till handling off, the
abraded material is an incision added to the erosion–deposition rate and flows
into the fluvial sediment system through the standard erosion path in all three
SPL flavours (`SPL`, `nlSPL`, `soilSPL`).

## 3. Glacial deposition — till and moraine

With `till.on` (default `True`), abraded rock is carried as **till** by the ice
and deposited as **moraine** where the ice melts out. Two conservation-guarded
modes:

- **Bulk bed** (no stratigraphy): the abraded volume is redistributed as a
  bed-to-bed transport, so the net bed-volume change is exactly zero (rock moved,
  not created).
- **Stratigraphic / dual lithology**: the till is removed from the stratigraphic
  layers it was abraded from (`erodeStrat`) and re-deposited as a fresh moraine
  layer in the ablation zone (`deposeStrat`), split into coarse/fine carrying the
  abraded (ice-mixed) fine fraction. Conservation is on the *solid* phase per
  fraction, and the bed bulks up by the uncompacted-till vs compacted-rock
  porosity contrast.

Independently, `till.route` (default `True`) routes the till down the ice-surface
flow network and melts it out toward each catchment's terminus (correct on
multi-glacier / global domains); `False` spreads the global abraded volume across
the ablation zone weighted by the meltwater rate. Both conserve mass.

**Meltwater** is **discharge-conserving** by default (`melt_conserve: True`): the
precipitation that fell as ice above the ELA is released as liquid meltwater
where the ice melts out, so total meltwater == total accumulation, and it is
re-injected into the river flow accumulation. With `melt_conserve: False` the
local precipitation-scaled ablation rate `m = max(−mdot,0)·area·[H>0]` is used
instead.

## 4. Loading — existing flexure leveraged

No flexure source was modified. goSPL's flexural isostasy (the parallel FV
biharmonic for planar models, spherical-harmonic shell for global) receives the
diagnostic ice thickness as the ice load: the change in thickness between steps
is converted to an equivalent load (scaled by `ρ_i/ρ_c`) and applied through the
existing path, so a growing ice sheet drives isostatic subsidence and
deglaciation drives rebound.

---

## Outputs

| Field | File / attribute | Meaning |
|---|---|---|
| `iceH` | mesh HDF5 / XDMF | ice thickness (m); restored on restart |
| `iceUb` | mesh HDF5 / XDMF | basal sliding speed (m/yr); abrasion driver |
| `iceMelt` | mesh HDF5 / XDMF | ablation meltwater (m³/yr) re-injected into the rivers |
| `iceAbr` | mesh HDF5 / XDMF | glacial abrasion rate `E_g = Kg·|u_b|^l` (m/yr); zero when abrasion is off |

## Configuration

The `ice` section is the opt-in. The valid keys are:

```yaml
ice:
    hela: 1850.0          # equilibrium-line altitude (m)
    hice: 2100.0          # ice-cap altitude (m)
    hterm: 1700.0         # glacier terminus (m); floor = max(hterm, sea level),
                          #   defaults to sea level if omitted  [or evol: <csv>]
    icedir: 1             # MFD flow directions for the ice routing
    eheight: 0.25         # Bahr thickness factor
    fwidth: 1.5           # Bahr width factor
    melt: 1.0             # ablation in net balance (0 = accumulation-only; 1 = true net; >1 shorter tongues)
    slide: 1.0e-3         # basal sliding coefficient (Glen sliding law)
    glen:  3.0            # Glen sliding exponent n
    accum_factor: 1.0     # precipitation -> ice accumulation fraction
    accum_max: 2.0        # (optional) cap on accumulation rate (m ice/yr)
    melt_conserve: True   # discharge-conserving river meltwater (default)
    abrasion:
        Kg: 1.0e-4        # vertical abrasion coefficient (0 = off)
        l:  1.0           # sliding-velocity exponent
        Kl: 0.0           # lateral (valley-wall) erosion coefficient (0 = off)
        lat_l: 1.0        # lateral velocity exponent (defaults to l)
    till:
        on: True          # carry abraded rock as till / moraine (default True)
        route: True       # catchment-routed (default True); False = melt-spread
```

`hela`/`hice`/`hterm` accept a uniform scalar (as above), a per-vertex map
`[file, key]`, or a `glaciers` time series of either — see the user guide for
the spatial/time-varying ELA syntax used by global models.

## Validation

Regression tests guard the model (`tests/test_regression.py`,
`tests/fixtures/minimal_ice_*.yml`): ice off/on parsing; end-to-end run
invariants (H ≥ 0, finite, no ice below terminus/ELA, valid basal velocity);
glacial-abrasion sign (`E_g = Kg·u_b`) and no-op when `Kg=0`; bulk till
conservation (eroded == deposited, net bed-volume change ≈ 0); stratigraphic /
dual-lithology till (per-fraction solid conserved, moraine carries fine); and
ice-load flexural subsidence with `iceUb` written.

---

## Scope notes

- **ELA mass balance.** Accumulation/ablation is the elevation-ramp ELA proxy,
  not an energy-balance or PDD melt model — the established goSPL forcing,
  adequate for the long-term landscape-to-global mass budget the model targets.
- **Spatially / temporally varying ELA.** For global runs, `hela`, `hice` and
  `hterm` can each be a per-vertex map and/or a `glaciers` time series (like the
  precipitation `climate` block), so the ELA varies with latitude and through
  time. Uniform scalars and the `evol` CSV remain unchanged. The
  `scripts/ela_from_temperature.py` helper derives `hela`/`hice` maps from a
  paleo-climate temperature history by lapse-rate inversion (it sets the ELA
  *position*; ablation stays precipitation-scaled).
- **Till routing.** `till.route: True` (default) routes till down the ice-surface
  flow network to each catchment's terminus; `False` spreads it across the
  ablation zone weighted by the meltwater rate. Both conserve mass.
