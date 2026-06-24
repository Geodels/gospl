# goSPL refactor audit (read-only)

Scope: `gospl/` (~10.3k LOC across 22 Python files). I read every Python file end-to-end. The report cites files relative to the repo root.

## 1. INCONSISTENT PATTERNS

### 1.1 `MPIcomm` source: `MPI.COMM_WORLD` vs `petsc4py.PETSc.COMM_WORLD`

The two are *not* literally the same object even though they wrap the same underlying communicator — collectives invoked on each go through different code paths, and mixing them in the same module masks where MPI is actually being driven.

- **Pattern A — `MPIcomm = petsc4py.PETSc.COMM_WORLD`** (10 files):
  - `gospl/flow/flowplex.py:16`
  - `gospl/flow/iceplex.py:18`
  - `gospl/flow/pitfilling.py:30`
  - `gospl/sed/sedplex.py:16`
  - `gospl/sed/seaplex.py:23`
  - `gospl/sed/hillslope.py:18`
  - `gospl/sed/stratplex.py:16`
  - `gospl/eroder/SPL.py:16`
  - `gospl/eroder/nlSPL.py:17`
  - `gospl/eroder/soilSPL.py:18`
- **Pattern B — `MPIcomm = MPI.COMM_WORLD`** (4 files):
  - `gospl/tools/addprocess.py:20`
  - `gospl/tools/outmesh.py:16`
  - `gospl/mesher/unstructuredmesh.py:28`
  - `gospl/mesher/tectonics.py:22`
- `gospl/tools/inputparser.py:14` defines `MPIrank` only — no `MPIcomm` at all (yet collectives are called elsewhere).
- `gospl/model.py:90` does `MPIrank = MPI.COMM_WORLD.Get_rank()` (no PETSc).

Three of the files (`pitfilling.py`, `seaplex.py`, `stratplex.py`) **also** define `MPIsize` while most others don't, even though they use the same MPI patterns.

### 1.2 `petsc4py.init(sys.argv)` is called once per module, not once per program

Repeated calls happen at import time in **13** files (`gospl/tools/inputparser.py:13`, `gospl/tools/outmesh.py:13`, `gospl/tools/addprocess.py:17`, `gospl/mesher/unstructuredmesh.py:25`, `gospl/mesher/tectonics.py:19`, `gospl/flow/flowplex.py:14`, `gospl/flow/iceplex.py:16`, `gospl/flow/pitfilling.py:28`, `gospl/sed/*.py` (all four), `gospl/eroder/*.py` (all three)). PETSc tolerates repeats but the boilerplate hides where state actually gets created, and `sys.argv` is inspected by every module — switching to a single init in `gospl/__init__.py` would be safer.

### 1.3 YAML key reading: nested `try/except KeyError` blocks

There are **130** `except KeyError` blocks in `gospl/tools/inputparser.py` alone. Three different styles:

- **Pattern A — `try/except` per key with default** (most common, e.g. `_readSPL` 383-414):
```python
try:
    self.K = splDict["K"]
except KeyError:
    print("...", flush=True)
    raise ValueError("Surface Process Model: K coefficient not found.")
```
- **Pattern B — `try/except` returning a silent default** (e.g. `_extraDomain` 196-216):
```python
try:
    self.seaDepo = domainDict["seadepo"]
except KeyError:
    self.seaDepo = True
```
- **Pattern C — outer `try/except` wrapping the whole section with parallel defaults inside the except** (e.g. `_readHillslope` 423-470, `_readSoilInfo` 546-644). The defaults are written twice (once inside the inner `except KeyError`, once inside the outer one), so changing a default means editing two places. This is where `_readSoilInfo` currently has a real divergence: line 564 defaults `self.Hs` to `0.0` in the populated branch but never sets it in the outer `except` (line 634 sets `self.Hs = 0.0` too — same here, but `H0` is `0.7` in one branch and `1.0` in the other: lines 578 vs 636).

### 1.4 KSP solver creation: cached vs ad-hoc

- **Pattern A — module-cached and reused** (`gospl/flow/flowplex.py:169-186` `_ksp_main`, `:123-149` `_ksp_fallback`; `gospl/eroder/nlSPL.py:214-232`; `gospl/sed/hillslope.py:215-233`):
```python
if self._ksp_main is None:
    ksp = petsc4py.PETSc.KSP().create(...)
    ...
    self._ksp_main = ksp
```
- **Pattern B — created+destroyed per call** (`gospl/eroder/SPL.py:192-243` `_coupledEDSystem`; `gospl/sed/seaplex.py:278-329` `_depMarineSystem`):
```python
ksp = petsc4py.PETSc.KSP().create(petsc4py.PETSc.COMM_WORLD)
...
ksp.destroy()
```

Both patterns assemble a fresh `Mat().createNest(...)` per call. The cached pattern is newer; the per-call pattern in the nested-matrix solves is a real performance asymmetry, not just a style difference.

### 1.5 `forcing-block` time-window lookup: three near-identical bodies

The `tecdata` / `raindata` / `sedfacdata` / `tedata` time-window advance logic is duplicated. Examples that should be a single helper:

- `gospl/mesher/tectonics.py:62-76` (`tecNb` advance + load `iloc` columns):
```python
nb = self.tecNb
if nb < len(self.tecdata) - 1:
    if self.tecdata.iloc[nb + 1, 0] < self.tNow + self.dt:
        nb += 1
if nb > self.tecNb or nb == -1:
    if nb == -1: nb = 0
    self.tecNb = nb
    ...
```
- `gospl/mesher/unstructuredmesh.py:665-693` (`rainNb` advance, same shape):
```python
nb = self.rainNb
if nb < len(self.raindata) - 1:
    if self.raindata.iloc[nb + 1, 0] <= self.tNow:
        nb += 1
if nb > self.rainNb or nb == -1: ...
```
- `gospl/mesher/unstructuredmesh.py:715-732` (`sedfactNb`, same shape).
- `gospl/tools/addprocess.py:150-176` (`teNb`, same shape).

Four copies × ~25 lines each. They differ subtly: `tectonics.py:64` uses `< self.tNow + self.dt`, the rain version (`unstructuredmesh.py:667`) uses `<= self.tNow`, and the Te version (`addprocess.py:152`) uses `<= self.tNow`. The drift in comparison operator could shift the timing of forcing events by one `dt` between subsystems.

### 1.6 Receiver-array snapshotting: `*i` (initial) duplicates

In `gospl/flow/flowplex.py:421-426`, `flowAccumulation` snapshots six arrays after the first `_buildFlowDirection` call so they remain available even after downstream routing rebuilds the matrix:
```python
self._buildFlowDirection(hl, False)
self.wghtVali = self.wghtVal.copy()
self.rcvIDi   = self.rcvID.copy()
self.distRcvi = self.distRcv.copy()
self.fMati    = self.fMat.copy()
self.lsinki   = self.lsink.copy()
```

The `i` suffix is the only thing telling consumers (`SPL.py:79`, `SPL.py:88`, `nlSPL.py:64`, `soilSPL.py:101`, `sedplex.py:67`) that this is "initial / unfilled-elevation" and not the current matrix. There's no docstring or constant for it.

### 1.7 Pit/sink ID transfer pattern is hand-rolled, slightly different each time

Several functions construct a per-pit float array, gather it via `Allreduce`, and inverse-index it back:

- `gospl/sed/sedplex.py:264-283` (`_spillCoords`, SUM reduce).
- `gospl/sed/sedplex.py:295-347` (`_addPitMicroTilt`, two SUM reduces for sum and count).
- `gospl/sed/sedplex.py:382-561` (`_diffuseLargePit`, multiple SUM/MAX reduces, including the two-step trick at `:442-471`).
- `gospl/flow/pitfilling.py:299-333` (`_getPitParams`, three Allreduce calls SUM/MAX/MAX).
- `gospl/flow/flowplex.py:330` and `gospl/sed/sedplex.py:112` (both Allreduce per-pit volume).

Same shape, hand-coded variants. Worth a `_perPitReduce(values, op)` helper.

### 1.8 Sentinel value for "no data" before `Allreduce(MAX)`

Three different magic numbers used as the "no data" floor for a MAX reduction:

- `-1.0e8` (most common): `gospl/tools/addprocess.py:509`, `:588`, `gospl/mesher/tectonics.py:430`, `:436` (mixed with `-1.0e10`).
- `-1.0e6` (boundary marker): `gospl/sed/seaplex.py:153`, `gospl/flow/iceplex.py:76`.
- `-1.0e10`: `gospl/mesher/tectonics.py:436` (gED), `:442` (gFI), `:449` (gSL).

The three pre-fill values are NOT interchangeable — `_distribute*` paths sometimes treat `-1e8` as "data is missing" — but the script-by-script choice is opaque.

### 1.9 Convergence threshold: hand-tuned per routine

- `gospl/sed/sedplex.py:184` `max_iters = 5000`, threshold = `1.0e-3` m³.
- `gospl/sed/seaplex.py:360-363` `max_iters = 5000`, threshold = `max(1.0, 1e-6 * initial_total)`.
- `gospl/flow/flowplex.py:363` threshold = `> 1.0e-3` and `gospl/flow/flowplex.py:375` skip-condition is `self.maxarea[0]`.
- `gospl/flow/pitfilling.py:217-223` `stp > 1000` fall-through.

Each value is a comment-explained heuristic; none are exposed in the YAML.

## 2. IMPLICIT CONVENTIONS

### 2.1 The `Model` god-class assembles 16 mixins; init order is load-bearing

`gospl/model.py:93-110` declares `Model` as a multi-inheritance of 16 mixin classes. `gospl/model.py:126-198` `__init__` calls each parent's `__init__` *by name*, in a specific order:

```
_ReadYaml → _STRAMesh → _VoroBuild → _UnstMesh → _WriteMesh
→ _FAMesh → _IceMesh → _SPL → _nlSPL → _soilSPL → _PITFill
→ _SEDMesh → _hillSLP → _SEAMesh → _GridProcess → _Tectonics
```

Convention: each class's `__init__` only allocates Vec/Mat from `self.hGlobal/hLocal` and assumes anything earlier in the list has populated `self.*`. **What goes wrong if you don't know**: re-ordering mixins in the class declaration looks safe (Python's MRO is consistent), but the `__init__` calls are explicit-by-name; permuting one line silently uses an unallocated attribute (e.g. `_FAMesh.__init__` reads `self.hLocal` set by `_UnstMesh._buildMesh`; if `_UnstMesh` ran after, you'd `AttributeError`).

### 2.2 PETSc Vec ownership is shared, not encapsulated

Every mixin's `__init__` does `self.foo = self.hGlobal.duplicate()`. Names like `tmp`, `tmpL`, `tmp1`, `dh`, `h`, `hl`, `newH` are **scratch vectors reused by every module**. `gospl/sed/sedplex.py:30-35` allocates `tmp`/`tmpL`/`tmp1`/`Qs`/`QsL`/`nQs`. Once that runs, *everyone* (eroder, sed, flow) reuses them. **What goes wrong if you don't know**: if you add a new method that uses `self.tmpL` and call it between two consumers that also use `self.tmpL`, you silently overwrite their state. There is no convention for "your code may zero this" vs "this is sacred".

### 2.3 `seaID` is recomputed at multiple entry points

- `gospl/flow/flowplex.py:252` (every flow accumulation).
- `gospl/sed/hillslope.py:274` (every hillslope step).
- `gospl/sed/seaplex.py:431` (`self.sinkIDs = self.lFill <= self.sealevel` — close kin but different array).

**What goes wrong**: a function that runs between an erosion step and a deposition step and looks at `self.seaID` will see *stale* values relative to the current `self.hLocal`. The class fields read like global state but are actually written cooperatively.

### 2.4 `idBorders` zeros out border arrays everywhere; `flatModel` guards it

Every numerical kernel checks `if self.flatModel:` before forcing border indices to zero or `-1e8`. `unstructuredmesh.py:423` sets `flatModel = True` only when there are boundary points anywhere across all ranks. **What goes wrong if you don't know**: a kernel written for the spherical case will quietly produce wrong values at the domain edges on 2D grids because the flatness clamp wasn't applied. There's no decorator or hook — just convention.

### 2.5 The forcing `DataFrame` row layout is positional and undocumented

`gospl/mesher/tectonics.py:73-117` does `self.tecdata.iloc[nb, 1]`, `.iloc[nb, 2]`, `.iloc[nb, 3]`, `.iloc[nb, -1]`. Column names exist (`"start", "end", "tMap", "zMap", "hMap"`) and are set in `inputparser._storeTectonics` (`gospl/tools/inputparser.py:788`), but the *consumer* uses integer positions. **What goes wrong**: add a column in `_storeTectonics` and every consumer silently picks up the wrong column. Same issue in `unstructuredmesh._updateRain:677-688` (`.iloc[nb, 4]` is `rMap`, `.iloc[nb, 5]` is `rKey` — but only if `_defineRain` appends columns in that order, which it does because line 1239 hard-codes it).

### 2.6 `-1e8` is "uninitialised"; `1.0e-8` is the discharge floor; `1.0e-3` is the deposit floor

Recurring magic numbers with no shared constant:
- `-1e8` — missing-data sentinel for `Allreduce(MAX)`.
- `1.0e-8` — discharge floor for output (`outmesh.py:288, :299, :328, :357`).
- `1.0e-3` (1 mm) — "drop sub-mm deposit" rule (`seaplex.py:465`, `sedplex.py:134`).
- `1.0e-1` (10 cm) — "bedrock exposed" threshold (`soilSPL.py:166`).
- `1.0e6` — bedrock-thickness sentinel (`stratplex.py:99`, `stratplex.py:194-225`).

The bedrock sentinel `1e6` in particular is documented in `stratplex.py:96-101` but used as a literal `1.0e6` four times in `erodeStrat` (`stratplex.py:194, 225`) with no constant. Renaming it requires grepping for `1.0e6` everywhere.

### 2.7 N/S boundary swap in gFlex setup — RESOLVED (gFlex removed)

This item is obsolete. gFlex and its regular-grid flat-model flexure were removed; the flat model now uses the parallel FV biharmonic solver (`method='fem'`), which maps the boundary sides geographically (N=ymax, E=xmax, S=ymin, W=xmin) with no swap. There is no longer a `simflex.BC_S/BC_N` assignment in `addprocess.py`.

### 2.8 KSP/SNES failure is logged but execution continues with a zero or stale solution

- `gospl/flow/flowplex.py:144-147` returns `vector2.set(0.0)` after divergence.
- `gospl/sed/seaplex.py:296-313` "Marine deposition may be inaccurate; operators should monitor this warning".
- `gospl/eroder/SPL.py:210-228` same: "downstream newH may be inaccurate".

**What goes wrong if you don't know**: a long run silently produces zero-discharge or stale-deposit fields after a divergence; downstream stratigraphy/flexure compounds the error.

### 2.9 Lazy import via `if "READTHEDOCS" not in os.environ`

Every module wraps its Fortran/heavy-deps imports in this guard, and `gospl/__init__.py` does it for the entire model class. **What goes wrong**: code that imports `gospl` in a docs CI job gets a no-op `Model` class. A test runner that sets that env var by mistake will silently fail to test anything. There is no warning.

### 2.10 Matrix assembly: which rows a rank may set (ghost-row `INSERT` hazard)

`Mat.setValuesLocal(row, …)` / `setValuesLocalCSR` map the *local* row index through the row local-to-global map. A **ghost** local row (a halo node owned by another rank) maps to that owner's global row, so setting it ships the value off-process. Whether that is correct depends on the stencil and the insert mode — and getting it wrong does **not** hang: it silently produces a partition-dependent operator (serial always passes; np>1 gives a wrong-but-finite result). There are three safe patterns in the codebase and one trap:

- **Diagonal-only** (`flowplex._matrix_build_diag`, the `seaplex` `createNest` sub-blocks): row `i` has only entry `(i,i)` with a per-node value; ghost rows carry the same synced value, so it is consistent under any mode.
- **Receiver-based** (`flowplex.matrixFlow`, `eroder/SPL.py`, `eroder/nlSPL._form_jacobian`, `flow/iceplex`, `seaplex` off-matrix): row `i` has entries at `i` and its flow **receivers**, with values built from per-node-`i` quantities (`rcvIDi`/`wghtVali`/`distRcvi`/slope) that are **partition-invariant** (the validated drainage arrays — see AGENTS.md "Mechanism #2"). The diagonal is *not* a sum over the full neighbour ring, so a ghost row computes the **same** values as the owner; assembling over all `lpoints` is therefore safe (and `matrixFlow` is the most-validated example).
- **Additive FV-Laplacian** (`hillslope._assembleDiffMat` via per-direction `axpy`, `_assembleDiffMatCSR` with `addv=ADD_VALUES`, `seaplex` off-matrix `ADD_VALUES`): the diffusion stencil's diagonal *is* `−Σ` neighbour terms, but `ADD_VALUES` makes assembly **additive** — each rank contributes the edges incident on the nodes it can see and PETSc sums the off-process partial rows into the complete row at `MatAssemblyEnd`. Safe over all `lpoints`.
- **The trap — FV-Laplacian over all `lpoints` with `INSERT_VALUES`.** A neighbour-summing row needs the node's *complete* neighbour ring; a ghost node's ring is **incomplete** on the ghosting rank (some neighbours are 2-ring, absent), so its locally-computed row is wrong, and under `INSERT` it **races** with the owner's correct row → partition-dependent boundary rows. This was the **2026-06-24 soil-diffusion spike**: `_evalJacobianSoil` looped `range(self.lpoints)` and `INSERT`ed every row including ghosts, producing an isolated elevation spike at sub-domain boundaries (np=1 clean — a single partition has no ghost rows — np>1 spiked by the 2nd output). The diagnostic tell was that an **exact** MUMPS solve still spiked, proving the *matrix* not the solver was partition-dependent. Fix: restrict to owned rows (`for i, row in enumerate(self.glIDs)`), which the marine counterpart `_evalJacobianMardDiff` already did. The owner always has the complete ring, so owned-row-only `INSERT` is correct.

**Rule**: a new matrix assembly must use one of the three safe patterns. If the row is a neighbour-summing (Laplacian) stencil, either assemble it **owned-rows-only** (`self.glIDs`) with `INSERT`, or **additively** (`ADD_VALUES`) over all `lpoints` — never `INSERT` a Laplacian row over all `lpoints`. A full `setValues*` sweep on 2026-06-24 confirmed the soil Jacobian was the only violation.

## 3. NAMING INCONSISTENCIES

### 3.1 Class naming convention is incoherent

Six different conventions in one package:

| Class | File | Convention |
|-------|------|------------|
| `FAMesh` | flow/flowplex.py | PascalCase + Mesh suffix |
| `IceMesh` | flow/iceplex.py | PascalCase + Mesh suffix |
| `SEDMesh` | sed/sedplex.py | ALL-CAPS + Mesh suffix |
| `SEAMesh` | sed/seaplex.py | ALL-CAPS + Mesh suffix |
| `STRAMesh` | sed/stratplex.py | ALL-CAPS + Mesh suffix |
| `UnstMesh` | mesher/unstructuredmesh.py | Camel + Mesh |
| `SPL` | eroder/SPL.py | ALL-CAPS |
| `nlSPL` | eroder/nlSPL.py | mixedCase |
| `soilSPL` | eroder/soilSPL.py | mixedCase |
| `hillSLP` | sed/hillslope.py | mixedCase (note: `SLP` not `SLOPE`) |
| `PITFill` | flow/pitfilling.py | ALL-CAPS + Pascal |
| `VoroBuild` | mesher/meshfunc.py | PascalCase verb |
| `Tectonics` | mesher/tectonics.py | PascalCase |
| `GridProcess` | tools/addprocess.py | PascalCase |
| `WriteMesh` | tools/outmesh.py | Pascal verb-noun |
| `ReadYaml` | tools/inputparser.py | Pascal verb-noun |

### 3.2 Entry-point method naming is incoherent

The "do one step of X for the runProcesses loop" verb varies:

| Method | File | Verb form |
|--------|------|-----------|
| `flowAccumulation` | flow/flowplex.py:394 | noun |
| `iceAccumulation` | flow/iceplex.py:129 | noun |
| `fillElevation` | flow/pitfilling.py:549 | verb-noun |
| `erodepSPL` | eroder/SPL.py:334 | verb (no get) |
| `erodepSPLnl` | eroder/nlSPL.py:387 | verb + suffix |
| `erodepSPLsoil` | eroder/soilSPL.py:309 | verb + suffix |
| `sedChange` | sed/sedplex.py:628 | noun |
| `seaChange` | sed/seaplex.py:419 | noun |
| `getHillslope` | sed/hillslope.py:268 | get-noun |
| `getTectonics` | mesher/tectonics.py:49 | get-noun |
| `applyTectonics` | mesher/unstructuredmesh.py:634 | apply-noun |
| `applyForces` | mesher/unstructuredmesh.py:577 | apply-noun |
| `applyFlexure` | tools/addprocess.py:492 | apply-noun |
| `cptOrography` | tools/addprocess.py:561 | cpt-noun (abbrev) |
| `updatePaleoZ` | mesher/tectonics.py:607 | update-noun |
| `updateSoilThickness` | eroder/soilSPL.py:292 | update-noun |
| `getCompaction` | sed/stratplex.py:299 | get-noun |
| `diffuseSoil` | eroder/soilSPL.py:441 | verb-noun |

`getTectonics` (computes the current advection/uplift) and `applyTectonics` (adds `upsub * dt`) sound like aliases but do **different** things; `getHillslope` is closer in spirit to `applyFlexure` than to `getCompaction`.

### 3.3 Hillslope/erodibility coefficient prefixes

The same kind of physical quantity gets different prefixes by module:

| Coefficient | YAML key | Python attribute | Used in |
|-------------|----------|------------------|---------|
| Hillslope diffusion (aerial) | `hillslopeKa` | `self.Cda` | hillslope.py |
| Hillslope diffusion (marine) | `hillslopeKm` | `self.Cdm` | hillslope.py |
| Non-linear hillslope strength | `hillslopenl` | `self.K_nl` | hillslope.py |
| Critical slope | `hillslopeSc` | `self.K_sc` | hillslope.py |
| Iteration count | `hillslopeNb` | `self.K_nb` | hillslope.py |
| Non-linear coef (lake/marine) | `nonlinKm` | `self.nlK` | hillslope.py |
| SPL erodibility | `K` | `self.K` | SPL/nlSPL/soilSPL |
| Soil SPL | `soilK` | `self.Ksoil` | soilSPL.py |
| Ice erodibility | `Ki` | `self.Kice` | SPL/nlSPL |

`Cd*` vs `K_*` vs `nlK` vs `Ksoil` vs `Kice` — five conventions for "an erodibility constant".

### 3.4 `mfdreceivers` vs `mfdrcvrs` — two Fortran functions, indistinguishable by name

- `from gospl._fortran import mfdreceivers` — `flow/flowplex.py:12`, `flow/iceplex.py:13`, `eroder/SPL.py:12`.
- `from gospl._fortran import mfdrcvrs` — `sed/seaplex.py:18` only.

The latter is the **12-direction marine variant**. Calling the wrong one silently changes the flow-direction count from `flowDir` (default 8) to a hard-coded 12. No comment in either module explains the split.

### 3.5 Boolean feature flags

| Flag | Type | Meaning |
|------|------|---------|
| `self.flexOn` | bool | flexural isostasy on |
| `self.iceOn` | bool | glacial dynamics on |
| `self.oroOn` | bool | orographic rain on |
| `self.cptSoil` | bool | soil production on |
| `self.fast` | bool | skip flow accumulation |
| `self.nodep` | bool | skip deposition |
| `self.seaDepo` | bool | enable marine deposition |
| `self.flatModel` | bool | 2D flat (not spherical) |
| `self.memclear` | bool | aggressive GC |

`*On` (3), `cpt*` (1), `*Depo` (1), plus four naked flags. Trivial to forget which form to test for.

### 3.6 Decay-depth parameters in soil

`gospl/tools/inputparser.py:565-578` defines, in order: `Hs` (production decay), `h_star` (roughness length), `H0` (transport decay), `Sperc` (bedrock conversion). Three different decay-depth-like quantities (`Hs`, `h_star`, `H0`) with three different naming conventions in 13 lines.

### 3.7 `iloc[nb, -1]` — magic-index dataframe access

`gospl/mesher/tectonics.py:78`:
```python
if self.tecdata.iloc[nb, -1] != "empty":
```
The `-1` index relies on `hMap` being the last column written by `_storeTectonics` in `inputparser.py:788`. The same dataframe is also accessed by name in `inputparser.py:869` (`tecdata["start"]`). Two access patterns coexist in the same dataframe — risky.

## 4. HIGH-RISK MODULES

Ranked by likelihood that an edit there breaks something elsewhere.

### 4.1 `gospl/mesher/unstructuredmesh.py` (828 lines) — risk: very high

Owns the PETSc DMPlex (`self.dm`), the local/global mapping (`self.locIDs`, `self.glbIDs`, `self.lpoints`, `self.mpoints`, `self.lcoords`, `self.mCoords`), boundary IDs (`self.idBorders`, `self.northPts`, etc.), the FV connectivity (`self.FVmesh_ngbID`, `self.larea`, `self.maxnb`), and the forcing dispatch (`applyForces` calls `cptOrography` or `_updateRain`, plus `_updateEroFactor` and `applyTectonics`). Every other module reads `self.dm`, `self.locIDs`, `self.lpoints`, `self.idBorders`, `self.sealevel`, `self.bL`/`self.bG`. Plus `destroy_DMPlex` lines 738-799 names every Vec/Mat in the project by hand — adding a new Vec elsewhere requires editing this list or the cleanup is incomplete.

### 4.2 `gospl/flow/flowplex.py` (495 lines) — risk: very high

Owns `_solve_KSP` and `_solve_KSP2`, used by `nlSPL`, `soilSPL`, `SPL`, `seaplex`, `tectonics`, `sedplex`, `hillslope`. Also owns `_matrix_build` and `_matrix_build_diag` — every PETSc-Mat construction in the project. Also runs `_buildFlowDirection`, which writes to `self.rcvID`, `self.distRcv`, `self.wghtVal`, `self.lsink`, `self.fMat`, `self.seaID` — consumed by everyone. Changing the KSP type or the matrix-building convention silently changes solver behaviour for ten downstream modules.

### 4.3 `gospl/tools/inputparser.py` (1704 lines) — risk: high

Single source of every model parameter. The `_readErofactor`/`_readTeMap`/`_readRain`/`_readTectonics` set of methods build the forcing dataframes whose **column order** is consumed by `iloc[nb, k]` calls in `unstructuredmesh.py:677-688`, `tectonics.py:78-117`, `addprocess.py:160-174`. Adding a column anywhere breaks consumers silently.

Also contains the **`rUni`/`sUni` mismatch bug** at line 915: in `_defineErofactor` the dict is written with key `"rUni"` (line 915) but the DataFrame is constructed with column `"sUni"` (line 929), so the `sUni` column ends up NaN when `sMap is None`. This propagates to `_updateEroFactor` at `unstructuredmesh.py:725`:
```python
if pd.isnull(self.sedfacdata["sUni"][nb]):  # always True when sMap is None
    loadData = np.load(self.sedfacdata.iloc[nb, 2])  # tries to load None
```
i.e. uniform-erodibility events crash. Worth confirming with a regression test before any refactor.

### 4.4 `gospl/sed/sedplex.py` (653 lines) — risk: high

Three deposition paths inside `_updateSinks` (full-fill vs bottom-up vs diffuse-large-pit) with shared per-pit reductions. `_diffuseLargePit` calls `_diffuseImplicit` (from `hillslope.py`) and mutates `self.Cd`/`self.minDiff_vec` in-place. The bottom-up vs diffuse-large-pit path classifier (`gospl/sed/sedplex.py:585-593`) controls which method runs per-pit; subtle errors here are basin-specific and not caught by simple smoke tests.

### 4.5 `gospl/flow/pitfilling.py` (652 lines) — risk: high

Twelve Fortran kernels, an MPI graph algorithm, and a per-iteration pit-ID merge loop (`_transferIDs` 1000-iteration fall-through at line 217-223). All other modules consume `self.pitIDs`, `self.pitInfo`, `self.pitParams`, `self.lFill`, `self.flatDirs`. The `_performFilling` algorithm is hard to test in isolation because it relies on inter-rank halo exchange.

### 4.6 `gospl/mesher/tectonics.py` (640 lines) — risk: medium-high

`getTectonics` mutates `self.hdisp`, `self.upsub`, `self.paleoZ`, `self.plateStep`, `self.tec_IDs`, etc. The advection chooses among IIOE1/IIOE2/upwind/plate-interp using `self.advscheme` integer codes — a quiet way to break the integration is to advect elevation but forget to advect `cumED` (lines 274-300) or `Lsoil` (329-358) consistently.

### 4.7 `gospl/sed/stratplex.py` (393 lines) — risk: medium

`stratK`/`stratH`/`stratZ`/`phiS` arrays use a per-node "column" scheme with a `1.0e6` sentinel for layer 0 representing infinite bedrock. The arithmetic in `erodeStrat:194-225` adds and subtracts `1.0e6` to make the cumsum work. Refactoring this requires understanding why the sentinel cancels out (because cumThick and eroVal share the offset). Stratigraphy is also touched from advection (`tectonics._advectStrati`), compaction (`getCompaction`), and three SPL flavours (`erodeStrat`/`deposeStrat`).

### 4.8 `gospl/tools/outmesh.py` (713 lines) — risk: medium

`_outputMesh` writes ~12 HDF5 datasets and the XMF schema mirroring them. `readData` reads them back. Adding a new field requires touching at least 3 places (write, read, XMF schema). Restart compatibility lives entirely here — silent breakage means existing checkpoint files become un-restartable.

### 4.9 `gospl/sed/hillslope.py` (539 lines) — risk: medium

Holds the cached PETSc TS (`_ts_marine`) shared by every "non-linear diffuse on a sub-domain" caller (marine, lake, soil). `_diffuseImplicit` mutates `self.Cd`/`self.minDiff_vec` in-place, called from `_diffuseOcean` (this file) and `_diffuseLargePit` (sed/sedplex.py:535). Concurrent calls would race — fine today because `runProcesses` is serial, but easy to forget if you parallelise.

## 5. YOUR CONFUSION LOG (interfaces I would likely misuse cold)

These are the things I'd genuinely get wrong if you handed me the codebase and a ticket.

### 5.1 `Model` is a 16-mixin god-class with no `super().__init__()`

`gospl/model.py:126-198` manually invokes `_ReadYaml.__init__(self, filename)`, `_STRAMesh.__init__(self)`, etc. by class name. Coming from typical Python I'd reach for `super().__init__()` and break the chain instantly. The clue that this is intentional (not just legacy) is that each parent's `__init__` takes `*args, **kwargs` but only some use them. I would also assume the inheritance order in the class declaration (line 93-110) matched the call order in `__init__` (line 136-184). **It does not** — `_VoroBuild` is third in inheritance but called 4th; `_Tectonics` is 6th in inheritance but called 16th. The MRO matters for *method resolution*, the init order matters for *attribute population*, and they are deliberately different.

### 5.2 `self.tmp`, `self.tmpL`, `self.tmp1`, `self.dh`, `self.h`, `self.hl`, `self.newH` are scratch

Every module borrows these. Reading any one of `seaplex.py:267-271`:
```python
self.tmpL.setArray(sedflux / self.dt)
self.dm.localToGlobal(self.tmpL, self.tmp1)
self.h.pointwiseMult(self.hGlobal, self.areaGlobal)
self.h.scale(1. / self.dt)
self.tmp1.axpy(1., self.h)
```
I'd assume `self.h` is the persistent state vector (it's named like `hLocal`/`hGlobal`). Actually `self.h` is allocated in `hillslope.py:39` as a scratch vector. Misuse risk: read `self.h` for the elevation, get whatever the last hillslope step wrote.

### 5.3 `rcvIDi` is not "the i-th receiver index" — it's "initial receiver array"

I would parse `rcvIDi[:, k]` as `rcvIDi` (indexed by something) at column `k`, with `i` being a placeholder. Actually `i` means "initial" — see `flowplex.py:421-426`. The `k` is the flow-direction index. So the array is **2-D**: `[node, direction]`, snapshotted before pit filling. Used 9 times in `SPL.py`/`nlSPL.py`/`soilSPL.py`/`sedplex.py`. The same `i`-suffix means "initial" everywhere (`wghtVali`, `distRcvi`, `lsinki`, `fMati`), which is a real convention — but only learnable by reading `flowAccumulation`.

### 5.4 `tecdata.iloc[nb, -1]` is the horizontal-displacement column

`tectonics.py:78`:
```python
if self.tecdata.iloc[nb, -1] != "empty":
```
Reading this cold I'd assume `iloc[nb, -1]` is "the last column whatever it is" and break the moment I add a new column anywhere. The convention is that **column order in `_storeTectonics` is part of the API contract** — adding columns must happen at the END, never in the middle. Not stated anywhere.

### 5.5 The `_extraX` methods are not optional — they are mandatory continuations

`gospl/tools/inputparser.py` has `_readDomain` that calls `_extraDomain` which calls `_extraDomain2`. `_readHillslope` calls `_extraHillslope`. `_readFlex` calls `_extraFlex`. `_readOrography` calls `_extraOrography(oroDict)`. `_readIce` calls `_extraIce(...)`.

Reading the name I'd assume `_extraDomain` parses some optional keys. It does — but it ALSO sets `self.advscheme` (line 252), which is required everywhere. Renaming or deleting `_extraDomain` would silently set `advscheme` to whatever was there before (or `AttributeError`). The "extra" prefix is misleading; they are continuations of `_readDomain`/`_readHillslope`/`_readFlex`/`_readIce`/`_readOrography`.

### 5.6 The `outmesh.py` XMF schema mirrors `_outputMesh` by hand

`_outputMesh` (`outmesh.py:202-398`) writes datasets with names like `elev`, `erodep`, `EDrate`, `waterFill`, `fillFA`, `FA`, `iceFA`, `iceH`, `flexIso`, `soilH`, `sedLoad`, `uplift`, `rain`. The XMF generator `_save_DMPlex_XMF` (line 494-683) hand-codes each `<Attribute>` block. The names there are aliased: HDF5 `elev` becomes XMF `Z`, HDF5 `erodep` becomes XMF `ED`, HDF5 `sedLoad` becomes XMF `SL`, HDF5 `uplift` becomes XMF `vTec`, HDF5 `rain` becomes XMF `Rain`. The renaming is silent and visible only in Paraview. Adding a field means: write HDF5, write XMF block, and choose whether to rename. I'd miss the rename half the time.

### 5.7 `_buildFlowDirection` is called from BOTH flow accumulation and downstream-sediment distribution and mutates `self.lsink`, `self.rcvID`, etc.

`flowplex.py:239-297` writes to six `self.*` attributes. `sedplex.py:147-149` calls `_buildFlowDirection(self.sedFilled)` (with the post-deposit topography) — overwriting `self.rcvID`/`self.lsink` mid-iteration. Then `_updateSinks` reads `self.lsink` and the SPL kernels read `self.rcvIDi` (the cached pre-fill version, separately preserved). So **the current state of `self.rcvID` is undefined depending on what just ran**. Reading any SPL code I'd assume `self.rcvID` was current; the convention is "always use the cached `rcvIDi` for erosion, the live `rcvID` only inside `_distributeSediment`/`_moveDownstream`".

### 5.8 `Eb` is "erosion rate", **negative** for deposition

`SPL.py:311-321`:
```python
E = -self.tmp.getArray().copy()
E = np.divide(E, self.dt)
self.Eb.setArray(E)
```
And later `tmp.setArray(-Eb * self.dt)`. So `self.Eb` is **positive** for incision, **negative** for deposition, and the conversion to thickness flips the sign again. I would track this wrong on first read; it took me 3 passes through `erodepSPL` (lines 334-383) to be sure. Compounding: `EbLocal.axpy(1.0, self.tmpL)` accumulates contributions from SEAMesh and hillSLP and soilSPL within a single timestep, so its sign is "net erosion rate this step including marine/hillslope deposits".

### 5.9 The deposition `fDep` coefficient is capped at 0.99 in three places

`SPL.py:298`, `nlSPL.py:210`, `soilSPL.py:202`: all cap `self.fDep[self.fDep > 0.99] = 0.99` because `fDep ≥ 1` makes the coupled `(I - W^T)Q + (1-fDep)h` block singular. Three places, three identical lines. Refactoring one and missing the others would produce three different convergence behaviours.

### 5.10 `meltfac` is "a numerical sink amplifier", not a melt physical rate

`flowplex.py:438-444` comment:
> `meltfac` is deliberately NOT applied here because it is a numerical sink-amplifier inside the implicit ice solver, not a physical melt multiplier.

Yet the YAML key is `melt` (`inputparser.py:1599`) and the docstring there says "Melting factor adjustment". Reading the YAML doc I would set `melt: 10` thinking it controls ablation; it actually controls the solver's sink-pull strength. The right place to extract melt is the new `iceMeltL` field (`iceplex.py:50-209`), wired into `rainA` separately.

### 5.11 `_diffuseImplicit` is one cached TS shared across "marine", "lake", and (independently) soil

`hillslope.py:471-487` caches `self._ts_marine`. It is reused for marine deposition AND lake (large-pit) deposition — the `label` argument is print-only. `soilSPL.diffuseSoil` uses a *separate* cached TS `self._ts_soil`. So three different non-linear diffusions, two caches. I would assume "one TS per use case" and either over-cache (breaking lake by sharing with marine) or under-cache (creating TSes per call).

### 5.12 The `tedata` index `nb` advance bug-shape

`addprocess.py:151`:
```python
nb = self.teNb
if nb < len(self.tedata) - 1:
    if self.tedata.iloc[nb + 1, 0] <= self.tNow:
        nb += 1
if nb > self.teNb or nb == -1:
    if nb == -1: nb = 0
    self.teNb = nb
```

Initialised to `-1` in `unstructuredmesh.py:483`, and the check `nb > self.teNb or nb == -1` is the trigger for the load. The first iteration: `self.teNb == -1`, so `nb < len(...) - 1` is true; if the start time is in the past, `nb += 1` → `nb = 0`, and the loader runs. If the start time is in the future, `nb` stays at `-1`, and the next branch `if nb > self.teNb or nb == -1` fires anyway with `nb == -1`, setting `nb = 0` and reading `self.tedata.iloc[0, ...]`. So **the first load always happens at index 0 regardless of timing** — which is fine, but the code reads like "load when the new event becomes current". Coupled with the per-block sentinel-prepending logic in `_readTeMap:1165-1178`, the actual time-zero behavior depends on whether the input data started after `tStart`. I'd misread the loop control and not realise the sentinel row exists.

---

**One concrete thing I'd fix first** if you wanted to derisk a refactor: extract the `tecdata`/`raindata`/`sedfacdata`/`tedata` "advance + load" pattern into a single helper that takes a name-keyed dataframe (instead of `iloc`-positional access). That collapses §1.5, §2.5, §3.7, §5.12 into one location; it surfaces the `rUni`/`sUni` bug in §4.3 as a typed key mismatch instead of a silent NaN; and it makes the existing intentional patterns (the `teNb`/`rainNb` first-load contract, the `i`-suffix snapshot convention) explicit rather than implicit.
