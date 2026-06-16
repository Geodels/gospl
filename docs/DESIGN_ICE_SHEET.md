# DESIGN: diagnostic glacial model — flow proxy, glacial erosion, till deposition, loading

Status: **implemented / shipped** on `dev`.
Last updated: 2026-06-16.

Design of goSPL's glacial capability: a cheap, robust **diagnostic
glacial-erosion model**. Each goSPL step routes the ELA accumulation downhill
into an ice discharge, derives an ice thickness and a basal sliding velocity from
that discharge (one linear solve, no ice-thickness PDE time integration), and
from those drives velocity-based glacial erosion (vertical abrasion + optional
lateral valley-wall erosion / U-shaping), ice-transported till / moraine
deposition, discharge-conserving meltwater, and flexural ice loading. Read
alongside `AGENTS.md` and `docs/ICE_SHEET_SUMMARY.md`.

---

## 1. Background — why not full ice dynamics

A true **Shallow-Ice-Approximation (SIA)** ice-thickness solve
(`∂H/∂t = ṁ − ∇·q`, Glen's-law deformation + sliding) was implemented and then
**removed**. Two facts drove the decision:

1. **The `H≥0` free boundary is an obstacle problem, not a root-find.** At the
   ice margin `H=0` but the residual is `≠0`, so the implicit `F(H)=0` thickness
   solve has no root there; it **diverged on real continental runs** (and clamping
   `max(H,0)` injects mass). An explicit flux-limited scheme handles the boundary
   but adds substep machinery and CFL bookkeeping.
2. **An ice-dynamics solve over-thickens km-scale continental ice at coarse
   resolution.** The SIA diffusivity `D ∝ H^(n+2)|∇s|^(n-1)` is stiff; over
   goSPL's long (10²–10⁴ yr) landscape timesteps and km-scale Voronoi cells it
   produces unphysically thick continental ice — wrong for the *morphology of
   glacial erosion* goSPL targets, and expensive.

The diagnostic model below is **physical and robust at any resolution and over
goSPL's long timesteps**, captures the erosional morphology (deepened and
U-shaped valleys, moraines), and conserves rock and water budgets — which is why
it is the only glacial model goSPL ships.

## 2. The algorithm (per goSPL step, `gospl/flow/iceplex.py`)

1. **ELA surface mass balance** `ṁ` (m ice/yr): accumulation above the ELA
   `hela`, full accumulation reached at the ice-cap altitude `hice`, ablation
   below — `ṁ = P·min(1, (η − hela)/(hice − hela))`. Accumulation is scaled by
   `accum_factor` and optionally capped at `accum_max`; ablation amplified by
   `melt`. Degenerate-config guard zeroes ice when `hice ≤ hela` or the max
   surface lies below the ELA.
2. **Route accumulation downhill** on the epsilon-filled bed via a
   multiple-flow-direction (MFD) algorithm (`icedir` directions) — the same
   flow-matrix / KSP machinery as the river flow accumulation — into an **ice
   discharge** `Q` (m³/yr). One linear solve; no time integration.
3. **Ice thickness** from a Bahr discharge scaling:
   `H = eheight · fwidth · Q^0.3`.
4. **Basal sliding velocity** from Glen's sliding law on that thickness and the
   bed-surface slope: `u_b ∝ H^(n-1)·|∇s|^(n-1)·∇s` (`s = z_bed + H`; `slide`,
   `glen`; the `ice_velocity` Fortran kernel) — physically bounded.
5. **Velocity-based glacial abrasion** `E_g = Kg·|u_b|^l` (vertical, deepens
   valleys) plus an optional **lateral valley-wall erosion**
   `Kl·u_b_neighbour^lat_l` (widens valleys toward a U-profile; the wall cell is
   tapered by its contact with the neighbouring ice column). Off by default
   (`Kg=0`, `Kl=0`).
6. **Glacial till**: abraded rock carried by the ice and deposited as **moraine**
   where the ice melts out; conserves rock volume; couples to stratigraphy and
   dual lithology (coarse/fine). `till.on` and `till.route` both default `True`.
7. **Glacial meltwater** to rivers: discharge-conserving by default
   (`melt_conserve: True`) — accumulation is released as meltwater where the ice
   melts out so total meltwater == total accumulation; re-injected into the river
   flow accumulation.
8. **Ice loading**: the ice thickness feeds the existing flexural isostasy.

## 3. State

| Field | Type | Holds |
|---|---|---|
| `self.iceH` / `iceHLocal` | global/local Vec | Bahr ice thickness `H` (restored on restart) |
| `self.iceUbL` | local array | basal sliding speed `\|u_b\|` for abrasion / output |
| `iceMeltL` | local | meltwater re-injected into `flowAccumulation` |
| `iceFlex` / `iceFAL` | — | flexure snapshot / ice flow-accumulation field |

Persistent Vecs are registered in `destroy_DMPlex` (`unstructuredmesh.py`).

## 4. Invariants / contracts (AGENTS.md)

- **Rock-volume conservation (till).** Abraded rock is moved, not created: bulk
  mode nets to zero bed-volume change; stratigraphic / dual-lithology mode
  conserves the *solid* phase per fraction (fine deposited == fine eroded), with
  the bed bulking up only by the uncompacted-till vs compacted-rock porosity
  contrast. Guarded by dedicated tests (the dual-lithology lesson: a per-fraction
  leak the total budget can't see must have its own guard).
- **Water-budget closure (meltwater).** With `melt_conserve: True`, total
  re-injected meltwater equals total ELA accumulation — the precipitation that
  fell as ice is released as liquid where the ice melts out, so downstream basins
  don't under-predict discharge.
- **MPI.** Routing and till-routing use the MPI-correct flow-matrix / KSP
  machinery; abrasion/velocity use rank-local `FVmesh` neighbour arrays with
  halo-synced fields. All reductions collective.
- **Soil-SPL compatibility.** When till handling is off, `E_g` enters the
  erosion–deposition rate as an incision through the standard erosion path in all
  three SPL flavours (`SPL`, `nlSPL`, `soilSPL`) with the existing `Eb`
  thickness-rate sign convention.
- **Ice-off byte-identical.** Without the `ice` section goSPL is unchanged;
  guarded by parity tests.

## 5. Config (`ice` YAML block — the only valid keys)

```yaml
ice:
    hela: 2000.          # equilibrium-line altitude (m); scalar, [file,key] map,
                         #   or via `glaciers` time series / `evol` CSV
    hice: 3000.          # ice-cap altitude (m): full accumulation above it
    hterm: 1500.         # terminus floor (m); effective floor = max(hterm, sea level);
                         #   omitted -> sea level
    icedir: 1            # MFD flow directions for the ice routing
    eheight: 0.25        # Bahr thickness factor
    fwidth: 1.5          # Bahr width factor
    melt: 10.            # ablation amplifier
    slide: 1.0e-3        # basal sliding coefficient (Glen sliding law)
    glen: 3.0            # Glen sliding exponent n
    accum_factor: 1.0    # precipitation -> ice accumulation fraction
    accum_max: 2.0       # (optional) cap on accumulation rate (m ice/yr)
    melt_conserve: True  # discharge-conserving river meltwater (default)
    abrasion:
        Kg: 1.0e-4       # vertical abrasion coefficient (0 = off)
        l: 1.0           # sliding-velocity exponent
        Kl: 0.0          # lateral (valley-wall) erosion coefficient (0 = off) -> U-shaping
        lat_l: 1.0       # lateral velocity exponent (defaults to l)
    till:
        on: True         # carry abraded rock as till -> moraine (default True)
        route: True      # catchment-routed deposition (default True);
                         #   False = global melt-weighted spreading
    # `evol` / `glaciers` time series and per-vertex ELA maps are also supported
```

Parsed in `inputparser._extraIce` (mandatory-continuation contract).

## 6. Spatially / temporally varying ELA

For global runs a single ELA is unphysical (~5000–6000 m in the tropics, near
sea level at the poles). `hela`, `hice` and `hterm` can each be a per-vertex map
`[file, key]` and/or a `glaciers` time series (mirroring the precipitation
`climate` block); the mass-balance ramp is evaluated per node with its local ELA.
Uniform scalars and the `evol` CSV remain available. The
`scripts/ela_from_temperature.py` helper derives `hela`/`hice` maps from a
paleo-climate temperature history by lapse-rate inversion (it sets the ELA
*position*; ablation stays precipitation-scaled). See `docs/user_guide/surfproc.rst`.

## 7. Outputs

| Field | Meaning |
|---|---|
| `iceH` | ice thickness (m); restored on restart |
| `iceUb` | basal sliding speed (m/yr); abrasion driver and dynamics diagnostic |
| `iceMelt` | ablation meltwater (m³/yr) re-injected into the rivers |
| `iceAbr` | glacial abrasion rate `E_g = Kg·\|u_b\|^l` (m/yr); zero when abrasion is off |

## 8. Validation

Regression tests (`tests/test_regression.py`, `tests/fixtures/minimal_ice_*.yml`):
ice off/on parsing; end-to-end run invariants (H ≥ 0, finite, no ice below
terminus/ELA, valid basal velocity); glacial-abrasion sign (`E_g = Kg·u_b`) and
no-op when `Kg=0`; bulk till conservation (eroded == deposited, net bed-volume
change ≈ 0); stratigraphic / dual-lithology till (per-fraction solid conserved,
moraine carries fine); ice-load flexural subsidence with `iceUb` written.
