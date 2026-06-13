# DESIGN: SIA ice-sheet model — dynamics, glacial erosion, till deposition, loading

Status: **design / not yet implemented**. Target branch: `feature/ice-sheet-sia` → `dev`.
Last updated: 2026-06-13.

Implementation plan for upgrading goSPL's glacial capability from the current
**MFD flow-routing proxy** to a physically-based **Shallow-Ice-Approximation
(SIA) ice-sheet model**, with velocity-based glacial erosion, ice-transported
till / moraine deposition, and ice loading. Read alongside `AGENTS.md` and the
existing ice audit. Mirrors the incremental, opt-in, conservation-guarded
delivery used for dual lithology (`docs/DESIGN_DUAL_LITHOLOGY.md`).

---

## 1. Locked decisions

| Decision | Choice |
|---|---|
| Ice dynamics | **True SIA**, solved **implicitly** (semi-implicit SNES) — ice-thickness evolution `∂H/∂t = ṁ − ∇·q`, wiring the dormant Glen's-law `ice_flux` Fortran kernel (deformation + sliding). Implicit is required at goSPL's km / 10²–10⁴ yr scale; explicit is reference-only (see §3). |
| Glacial erosion | **Velocity-based abrasion** `E = Kg·|u_b|^l` using the SIA basal/sliding velocity. |
| Glacial deposition | **Ice-transported till** — abrasion produces till carried with the ice flux, deposited where ice thins/melts (terminus moraine, basal till), then reworked by meltwater. |
| Surface mass balance | **Keep the ELA model** (`hela`/`hice` elevation bands; existing `iceAccumulation` source). |
| Delivery | **Opt-in** `ice: flow_model: sia` (default `mfd` = current proxy). Existing ice runs byte-identical; SIA selectable, validated incrementally, can become default later. |

## 2. What already exists (reuse, don't rebuild)

- **`ice_flux(nb, hice, zbed, ad, as, glen, val)`** (`fortran/functions.F90`, bound in `functions.pyf:16` but **never imported**): per-cell SIA flux divergence from deformation velocity `u_d = ad·H^(n+1)·|∇s|^(n-1)∇s` and sliding `u_s = as·H^(n-1)·|∇s|^(n-1)∇s`, `s = zbed + H`, Glen exponent `glen`. This is the SIA core — wire it in.
- **ELA mass balance**: `iceplex.iceAccumulation` already builds the elevation-dependent accumulation/ablation source (`hela`/`hice`); reuse it as the SMB `ṁ`.
- **Nonlinear implicit-diffusion solver pattern**: SIA is nonlinear diffusion of `H` (`∂H/∂t = ∇·(D∇s)+ṁ`, `D ∝ H^(n+2)|∇s|^(n-1)`), structurally identical to the non-linear hillslope (`hillslope._hillSlopeNL` SNES `ngmres`, `_diffuseImplicit` TS `rosw`, `_form_residual_nl_hillslope`/`fctcoeff`). Reuse this CACHED-solver pattern (AGENTS.md > KSP/SNES/TS lifecycle), registering any new cached solver in `destroy_DMPlex`.
- **Ice loading**: `addprocess.applyFlexure` already converts ice-thickness change to an equivalent sediment load (×910/`rhoc`). SIA just supplies a better `iceHL`.
- **Meltwater coupling**: `iceMeltL` re-injection in `flowplex.flowAccumulation` — reuse for till reworking.

## 3. Physics

**Ice-thickness evolution** (per cell, finite-volume on the Voronoi mesh):

> `∂H/∂t = ṁ − ∇·q`,  `q = (u_d + u_s)·H`,  `s = z_bed + H`

`ice_flux` returns `∇·q` per cell. The diffusivity `D ∝ H^(n+2)|∇s|^(n-1)` blows up where ice is thick/steep, so an **explicit** scheme has a severe CFL limit `Δt ≲ Δx²/(2 D_max)`.

**The solve MUST be implicit — this is decided, not optional.** At goSPL's
scales (km-scale Voronoi cells, landscape timesteps of 10²–10⁴ yr) the stable
*explicit* ice timestep for thick/steep continental ice is years–decades, so an
explicit subcycled SIA would need **10²–10³ substeps per goSPL step** — correct
but prohibitive. An **implicit** nonlinear-diffusion solve is unconditionally
stable in `dt` and takes the landscape timestep in one (iterated) solve.

- **Primary: semi-implicit / SNES.** Solve `H^{t+dt}` implicitly with the
  residual `F(H) = (H − Hₒₗₐ)/dt − ṁ + ∇·q(H)`, where `∇·q(H)` is evaluated by
  `ice_flux`. This is the **same construction as the non-linear hillslope**
  (`_hillSlopeNL`: SNES `ngmres`, residual from a Fortran kernel `hillslp_nl`,
  Jacobian-free — no analytic Jacobian needed); SIA reuses that pattern with
  `ice_flux` in place of `hillslp_nl`. SIA is stiffer (the `H^{n+2}`
  dependence), so a **lagged-diffusivity Picard** iteration (or SNES line-search
  damping) is the robustness fallback if `ngmres` struggles — this is the one
  real implementation risk, retired by validating against the analytical dome.
- **Reference only: explicit subcycling.** Keep a CFL-subcycled explicit path
  *for validation on small/short cases* (it makes the `ice_flux` call trivial to
  check term-by-term), NOT for production.

**Unstructured-grid compatibility (confirmed).** `ice_flux` is written on the
`meshparams` finite-volume arrays (`FVarea`, `FVnNb`/`FVnID`, `FVeLgt`,
`FVvDist`) — the exact arrays the already-validated implicit diffusion kernels
(`hillslp_nl`, `fctcoeff`) use, populated once by `definetin`
(`unstructuredmesh.py:145`). It is a mass-conservative FV flux divergence
(`Σ_neighbours (u_d+u_s)·H̄·FVvDist / FVarea`, equal-and-opposite across shared
faces) on the **total surface** `∇(z_bed+H)`, so ice ponds/overflows closed
basins with no sink-filling. Two small gaps: (i) it returns only the divergence,
**not the basal velocity** — Phase 2 extends it (or adds a companion) to also
output per-cell `|u_s|`; (ii) `ad`/`as`/`glen` are **scalar** (uniform) — fine
for v1, a signature extension if spatially-variable sliding is wanted later.

**Velocity for erosion**: the sliding velocity `u_b = u_s` (basal) is the physically-relevant driver of abrasion. Expose a per-cell `|u_b|` from the SIA solve.

**Glacial abrasion**: `E_g = Kg·|u_b|^l` (default `l≈1`–2), applied where ice is present, added to the erosion step alongside fluvial/hillslope (sign/units per the `Eb` thickness-rate convention, AGENTS.md).

**Till**: abraded rock becomes a till volume that (i) is advected with the ice flux (same `ice_flux` direction field), (ii) deposits where ice thins/melts (∂H/∂t<0 or below terminus) as basal till / terminal moraine, (iii) the melt-zone till is then available to the meltwater/fluvial system. Conserve till mass per step (a glacial analogue of the dual-lithology fine-conservation guard).

## 4. State (allocated only when SIA on)

| Field | Type | Holds |
|---|---|---|
| `self.iceH` / `iceHLocal` | global/local Vec | SIA ice thickness `H` (replaces the Bahr-scaled proxy thickness when `flow_model: sia`) |
| `self.iceU` (or `iceUbL`) | local Vec/array | basal/sliding speed `|u_b|` for abrasion |
| `self.tillH` / `tillHLocal` | global/local Vec | mobile glacial till thickness (eroded-but-not-yet-deposited) |
| existing `iceMeltL`, `iceFlex`, `iceFAL` | — | reused (melt, flexure snapshot, optional proxy field) |

New persistent Vecs **must** be added to `destroy_DMPlex` (`unstructuredmesh.py`).

## 5. Config (`ice` YAML block — additive, defaults preserve current behaviour)

```yaml
ice:
    # existing: hela, hice, hterm, diff, Ki, melt, fwidth, eheight, icedir, evol
    flow_model: sia            # 'mfd' (default, current proxy) | 'sia'
    sia:
        Aglen: 1.0e-16         # Glen rate factor -> deformation coeff `ad`
        slide: 1.0e-3          # sliding coeff `as` (0 = no sliding -> no abrasion)
        glen:  3.0             # Glen exponent n
        cfl:   0.25            # ice-substep CFL number (scheme A)
        max_substeps: 500
    abrasion:
        Kg: 0.0                # glacial abrasion coeff (0 = off, like Ki today)
        l:  1.0                # sliding-velocity exponent
    till:
        on: True               # produce/transport/deposit till
```

Parsed in `inputparser._extraIce` (extend; mandatory-continuation contract).

## 6. Phased plan (each phase: opt-in-gated, tested, MFD-proxy & ice-off byte-identical)

| Phase | Scope | Files | Risk |
|---|---|---|---|
| **0** | Config: `flow_model`, `sia`/`abrasion`/`till` params; `iceSIA` flag. Default `mfd` → existing path bitwise-unchanged. | `inputparser.py`, AGENTS.md | low |
| **1** | ✅ DONE. SIA thickness evolution: `ice_flux` wired; **implicit `ngmres` SNES** (`_iceFlowSIAImplicit`, reusing the `_hillSlopeNL` pattern, `_snes_ice` in `destroy_DMPlex`) is the default/production solver; explicit CFL-subcycled reference (`_iceFlowSIAExplicit`) selectable via `sia.solver`. ELA SMB; free boundary `H≥0`; terminus clamp. Validated: implicit ≡ explicit on a flux-active 500 m dome (rel ~2e-5). Selected by `sia.solver` (default `implicit`). *Quantitative Halfar/Bueler dome remains a deeper validation.* | `iceplex.py`, Fortran import, `destroy_DMPlex` | **high** |
| **2** | Expose basal/sliding velocity `|u_b|` per cell from the SIA solve. | `iceplex.py`, maybe `ice_flux` companion | medium |
| **3** | Velocity-based abrasion `Kg·|u_b|^l` in the erosion step (gated; replaces/【supplements】the `Ki·F^m` proxy under `flow_model: sia`). | `eroder/SPL.py`/`nlSPL.py`/`soilSPL.py` | medium |
| **4** | Till: production from abrasion → advect with ice flux → deposit at melt/terminus (moraine) → hand melt-zone till to meltwater. Per-step till mass-conservation guard. | `iceplex.py`, `flowplex.py`, `sedplex` hand-off | **high** |
| **5** | ✅ DONE. Loading works through the **existing** `applyFlexure` unchanged: the `iceFlex` snapshot is taken before the SIA update, so `dIce = iceHL − iceFlex` (×910/`rhoc`) loads the plate with the SIA thickness. Verified: SIA+global flexure produces isostatic subsidence under ice. | `addprocess.py` (no change) | low |
| **6** | ✅ DONE. Output: `iceUb` (basal velocity) added (`iceH` already written). Tests: SIA invariants, implicit≡explicit, glacial-abrasion sign, till volume conservation, **ice-volume conservation under zero SMB (Halfar-style flux check)**, flexure loading, MFD/ice-off parity. *Quantitative transient Halfar dome on a flat-bed fixture remains a deeper future validation.* | `outmesh.py`, `tests/` | medium |

**Recommended PR grouping:** PR-A = Phases 0–2 (dynamics + velocity, the new ice engine, validated against the dome); PR-B = Phases 3–4 (erosion + till); PR-C = Phases 5–6 (loading + output + tests). Stack like the dual-lithology PRs.

## 7. Invariants / contracts to honour (AGENTS.md)

- **CACHED solver lifecycle**: any new SNES/TS (scheme B) → register in `destroy_DMPlex`. New persistent Vecs (`iceH`, `iceU`, `tillH`) → same.
- **MPI**: `ice_flux` uses `FVmesh` neighbour arrays (rank-local); halo-sync `H` (`localToGlobal`/`globalToLocal`) before each flux eval and across substeps. All reductions collective.
- **`rcvIDi`/`Eb` sign conventions** unchanged; glacial erosion enters with the thickness-rate sign.
- **`flow_model: mfd` (default) and ice-off MUST stay byte-identical** — parity tests.
- **Mass conservation**: ice volume (`∂H/∂t` vs SMB − flux on a closed domain) and till volume, each guarded by a dedicated test (the lesson from dual lithology: a fraction/volume leak the total budget can't see must have its own guard).

## 8. Open questions for implementation time
- Implicit solver robustness for stiff SIA: `ngmres` vs lagged-diffusivity Picard vs SNES line-search (decide empirically against the dome). Implicit is required; explicit is reference-only.
- Whether SIA-mode abrasion replaces or supplements the MFD `Ki·F^m` erosion.
- Bedrock vs ice-surface gradient handling at margins / floating ice (marine-terminating glaciers — likely out of scope v1, grounded ice only).
- Till: single lithology now; couple to dual-lithology fractions later (out of scope v1).
