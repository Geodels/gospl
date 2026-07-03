# DESIGN: groundwater solute transport (duricrust "Level B" geochemistry)

Status: **PROPOSED (design only, not implemented).** The conservative-geochemistry
extension of the water-table + duricrust module (`DESIGN_WATERTABLE_DURICRUST.md`
§3a "Level B", §15). It closes the duricrust's chemical budget: **dissolve →
transport the solute along the groundwater flux `q = −T∇h` → precipitate at the
capillary fringe → export the remainder via baseflow**, with a domain mass-balance
guard. Ships disabled by default; opt-in, byte-identical to the current code when
off. Honors the invariants in `AGENTS.md` (MPI contract, KSP lifecycle,
scratch-vector contract, `destroy_DMPlex` registration).

Companion to `DESIGN_WATERTABLE_DURICRUST.md` (the host module — the head solve,
recharge, fringe, baseflow, `stratDuri` all pre-exist here), `DESIGN_PROVENANCE.md`
and `DESIGN_DUAL_LITHOLOGY.md`.

> **Key framing decision.** The solute is an **`n_species` array from the start**
> (default `1`). **Single-tracer ships first**; **multi-tracer** (a few independent
> lumped solutes — e.g. carbonate/Ca, silica, Fe — for duricrust *typing* and
> per-species ocean fluxes) is then a **config change, not a rewrite**, because
> every tracer shares the same transport / precipitation / export template.
> **Coupled multi-species aqueous equilibrium chemistry** (speciation, activity
> models, pH, mineral saturation via a PHREEQC-class solver) is **explicitly OUT
> OF SCOPE** — see §9.

---

## 1. Scope & decisions (locked)

| # | Decision | Choice |
|---|---|---|
| 1 | What Level B adds over Level A | Level A (`_weatheringSupply`) is **supply-only, local, non-conservative** — the crust precipitates from a local rate with no accounting. Level B makes the crust a **closed, transported budget**: solute is dissolved, advected down the water-table gradient, precipitated where it saturates, and the remainder exported. |
| 2 | Tracer model | **`n_species` independent lumped conservative tracers** (default `1`). Each is a single scalar concentration with its own dissolution/saturation/precipitation parameters. **No coupled speciation** (§9). |
| 3 | Transport | **Quasi-steady advection-reaction** on the DMPlex per step, `∇·(q c) = D − P`, with `q = −T∇h` the Darcy flux from the existing head solve. Steady (not transient) — solute transit ≪ Δt (§2). Reuses the FV upwind-advection kernels (`getfacevelocity` / `advecupwind`), same as orography/marine. |
| 4 | Precipitation | At the capillary fringe (`Φ`, the existing favourability), when the transported concentration exceeds a per-tracer saturation `c_sat`. Feeds the **existing** `duriH` / `duriF` (Level B replaces Level A's `Ψ` as the crust source when on). |
| 5 | Dissolution | Solute source from chemical weathering (reuse the Level-A rate `W` / `prodSoil` / recharge drivers), **debiting a conserved source pool** so mass balances. |
| 6 | Export | Solute leaving at the seepage nodes rides the **baseflow** already routed to rivers (`conserve_baseflow`, `baseflowL`) → a **dissolved flux to the ocean** (headline output). |
| 7 | Conservation | Domain guard `Σ dissolved − Σ precipitated − Σ exported − ΔS_storage ≈ 0`, MPI-reduced (like the dual-lithology fine mass-balance / baseflow budget). |
| 8 | Opt-in | `groundwater: geochem:` block. Absent ⇒ inert ⇒ bit-identical. Requires the water table on (`gwOn`); the duricrust (`duriOn`) for the precipitation coupling. |
| 9 | Out of scope | Coupled multi-species equilibrium, pH/activity models, porosity→permeability feedback on the head solve, redox. §9. |

---

## 2. Why it is scale-appropriate (quasi-steady, not transient)

goSPL runs at Δt ≈ 10²–10³ yr, Δx ≈ 500 m–km. The groundwater Darcy velocity
`v = q/(S·b)` is ~10¹–10² m/yr, so over one step a solute parcel travels ~km —
**comparable to Δx**. The solute field therefore reaches a quasi-equilibrium
distribution *within* a step, exactly like the Dupuit head. So the transport is
solved as a **steady advection-reaction** each step (no CFL, one linear solve),
the same justification and the same solver class as the head solve. Lateral
source→sink transport over km **is** resolvable at this scale; sub-km chemistry is
not (§9), which is why the tracers are lumped and precipitation is parameterized.

---

## 3. Governing equation (per tracer)

For each lumped tracer with groundwater concentration `c` (mass solute per unit
water volume), the steady reactive-transport balance is

```
∇·(q c) = D(x) − P(x, c)                         q = −T ∇h   (Darcy flux)
```

- **Transport** `∇·(q c)` — FV upwind advection of `c` by the groundwater flux
  `q` (from the head `h`). Assembled with the existing `getfacevelocity` +
  `advecupwind` kernels (the operator orography/marine already use); solved with a
  cached `gw_solute_` KSP (fgmres + hypre, like the head).
- **Dissolution `D`** — the solute source (mass/vol/yr): the chemical-weathering
  rate driving material into solution, from the Level-A driver
  (`_weatheringSupply` / `prodSoil` / recharge `R`), scaled per tracer by a
  `weatherability`. Debits a conserved source pool (regolith/rock) so mass is
  tracked.
- **Precipitation `P`** — the sink at the fringe: `P = k_p · Φ · max(0, c − c_sat)`
  (precipitates the super-saturated fraction where the water table sits in the
  capillary fringe `Φ`). This is the crust source term — it feeds `duriH` (below).
- **Export** — at the seepage nodes the outflow carries `c` out of the aquifer;
  summed with the `baseflowL` discharge it is the **dissolved flux to the surface
  network / ocean**.

The precipitated mass converts to a crust-thickness increment (a per-tracer solid
molar volume), so Level B's `P` **replaces** Level A's supply `Ψ` as the
`_updateDuricrust` formation term when geochem is on:

```
dduriH/dt |Level B = Σ_species P_species · v_solid_species          (m/yr)
```

Everything downstream of `duriH` (induration `duriF`, K-armoring, the `stratDuri`
record, exhumation) is **unchanged** — Level B only changes *where the crust
material comes from* (transported + conserved) not *what the crust does*.

---

## 4. Multi-tracer (the `n_species` array)

`c` is stored `(lpoints, n_species)` (default `n_species = 1`). Each species has
independent `weatherability_k`, `c_sat_k`, `k_p_k`, `v_solid_k`, so the same
transport/precip/export code runs per column of the array (a Python loop over
species, or a vectorised solve). Multi-tracer then gives, at ~`n_species`× the
single-tracer cost (a few extra advection solves, **no** coupled equilibrium):

- **Duricrust typing** — calcrete (carbonate/Ca tracer, arid discharge zones),
  silcrete (silica tracer), ferricrete (Fe tracer, humid Fe-mobilising) emerging
  from climate + lithology, recorded as a dominant-type field (like provenance's
  `dominant`).
- **Per-species dissolved ocean flux** — a weathering-derived alkalinity/Ca flux
  to the sea (weathering–climate–carbonate coupling; cf. the group's PNAS-2025
  carbonate-burial work).

Single-tracer is the identical code path with `n_species = 1`; multi-tracer is a
YAML count + per-species parameter rows, **not** a code rewrite.

---

## 5. State & coupling points (all gated on the opt-in)

**New state** (allocated only when the geochem block is on; in `destroy_DMPlex`):

| Field | Meaning | Kind |
|---|---|---|
| `self.gwSoluteL[:, k]` / `G` | groundwater solute concentration per tracer | Vec(s) / `(lpoints, n_species)` |
| `self.gwSourcePool[:, k]` | remaining dissolvable source mass per tracer (conservation) | numpy |
| `self.gwOceanFlux[k]` | running dissolved export to the sea per tracer (diagnostic) | scalar(s) |
| cached `self._ksp_solute`, `self._soluteMat` | the advection-transport solver | — |

**Reuses (no new machinery):** `q` from the head `h` (the `−(L·h)` flux already
computed for the lake coupling generalises to face fluxes via `getfacevelocity`);
the fringe `Φ` and `duriH`/`duriF` from `_updateDuricrust`; the seepage set and
`baseflowL` from `_baseflowClosure`; the FV advection kernels; the per-layer
record slot alongside `stratDuri` for **solute-source provenance** (§ unblocks
`DESIGN_WATERTABLE_DURICRUST.md` §11).

---

## 6. YAML opt-in

```yaml
groundwater:
    # ... existing hydrology + duricrust keys ...
    geochem:                 # absent ⇒ off ⇒ byte-identical
        species:             # 1 entry = single tracer (default); N = multi-tracer
          - name: carbonate
            weatherability: 1.0     # dissolution scaling (per tracer)
            c_sat: 1.0              # saturation concentration (precip threshold)
            precip_rate: 1.0        # k_p
            solid_volume: 1.0       # m crust per unit precipitated mass
        # - name: silica   { ... }   # add rows for multi-tracer typing
        # - name: iron     { ... }
        conserve: True       # domain mass-balance guard on (dissolve=precip+export+ΔS)
```

---

## 7. Outputs
`solute` (per-tracer groundwater concentration), `crust_type` (dominant tracer
where indurated, multi-tracer), and the per-tracer **`ocean_solute_flux`**
(dissolved export, m³·conc/yr). The precipitated mass already shows through the
existing `duricrust`/`induration` fields; solute provenance rides the `stratDuri`
record.

---

## 8. Conservation & tests
Domain guard (owned nodes, `Allreduce`): `Σ D·A − Σ P·A − Σ export − ΔS_pool ≈ 0`
per tracer, per step. Guard tests: `test_geochem_opt_in` (off ⇒ inert /
byte-identical); `test_geochem_conserves` (closed budget); `test_geochem_transport`
(solute moves downgradient — a source upstream precipitates a crust downstream,
which the local Level-A supply cannot do); `test_geochem_multitracer` (two tracers
give distinct crust types); parallel np=1-vs-2 invariance of the ocean flux.

---

## 9. Explicitly OUT OF SCOPE (and why)
- **Coupled multi-species aqueous equilibrium** — speciation (carbonate system,
  complexation), activity coefficients, pH, mineral saturation via a
  PHREEQC-class per-node Newton solve. It is a *category jump* (embedding a
  reactive-transport geochemical engine, e.g. `phreeqcRM`), and it is
  **scale-inappropriate here**: at km / annual-mean / My resolution the controls
  on speciation (seasonal wet/dry, evaporative concentration, soil-CO₂,
  microscale reactive area) are all sub-grid, so the added chemical fidelity
  cannot be fed meaningful inputs — high-precision chemistry on order-of-magnitude
  data. It also costs a per-node Newton solve × millions of nodes × 10⁴–10⁵ steps.
- **Porosity → permeability feedback** on the head solve (precipitation clogging
  pores). Deferred — a stiff two-way coupling; Level B keeps transport one-way
  from the head.
- **Redox / temperature-dependent equilibria** beyond the optional Arrhenius
  already in `_weatheringSupply`.

The lumped multi-tracer model is chosen precisely because it captures the
*scale-appropriate, observable* signal (source→sink transport, crust typing,
ocean flux) without pretending to resolve chemistry the grid cannot support.

---

## 10. Compatibility
Independent of and composable with dual lithology (`stratHf`), provenance
(`stratP`) and the existing duricrust — the solute fields are new, transport is
composition-only on the head, and precipitation feeds the existing `duriH` hook.
No change to the conservation invariants the other modules are guarded by.

---

## 11. Phased plan

| Phase | Deliverable | Guard test |
|---|---|---|
| G0 | **DONE.** `geochem:` parser + `gwGeochemOn` flag + per-species param lists; `_GWMesh` state alloc (`gwSolute`/`gwSourcePool` `(lpoints, n_species)`, `gwOceanFlux`, scratch `soluteL`/`soluteG`, cached `_ksp_solute`/`_soluteMat`), registered in `destroy_DMPlex`. Inert — nothing solved. | `test_geochem_opt_in` (off ⇒ inert; on ⇒ n_species state; inert run byte-identical) |
| G1 | Steady solute transport `∇·(q c)=0` (no reactions yet): assemble the advection operator from `q=−T∇h`, cached `gw_solute_` KSP; seepage outflow BC. | bounded/finite; np-invariant |
| G2 | Dissolution source `D` (from the Level-A driver, debiting `gwSourcePool`) + precipitation sink `P` at the fringe; feed `duriH`; domain mass-balance guard. | `test_geochem_conserves`, `test_geochem_transport` |
| G3 | Baseflow export → per-tracer `ocean_solute_flux` output. | export ≈ dissolved − precipitated at steady state |
| G4 | **Multi-tracer** (`n_species>1`): per-species params, `crust_type` dominant field, per-species ocean flux. | `test_geochem_multitracer` |
| G5 | Solute-source provenance (attribute crust to source area, riding `stratDuri`); unblocks `DESIGN_WATERTABLE_DURICRUST.md` §11. | provenance sums close |
| G6 | Docs: `tech_guide/groundwater.rst` geochem section, `surfproc.rst` block, `api_ref` page; AGENTS milestone. | docs build green |

---

## 12. Estimate
Comparable in scope to dual-lithology or the duricrust itself: G0–G3 (single
tracer, the useful core incl. the ocean flux) is the bulk; G4 (multi-tracer) is
small given the array-from-the-start design; G5 (provenance) is optional. A
separate feature branch + PR, developed after the water-table/duricrust PR lands.
