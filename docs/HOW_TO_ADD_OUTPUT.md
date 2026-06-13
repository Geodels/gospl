# How to add a new output field to goSPL

This runbook covers adding a new visualisable output field (a per-node scalar written every output timestep, viewable in Paraview). It assumes you have read [`AGENTS.md`](../AGENTS.md) — specifically the **Scratch vector contract** and **High-risk modules** sections.

Output is centralised in a single file: `gospl/tools/outmesh.py`. The hard part isn't writing the HDF5 dataset (3 lines); it's understanding which source-of-data pattern fits your field and whether you need a new Vec (which then has to be wired into `destroy_DMPlex` per AGENTS.md).

---

## Overview

goSPL writes outputs per timestep as **two coordinated artefacts**:

```
<outputDir>/
├── gospl.xdmf                # master time-series manifest (rank 0 only)
├── xmf/
│   ├── gospl0.xmf           # per-step XMF manifest (rank 0 only)
│   ├── gospl1.xmf
│   └── ...
└── h5/
    ├── topology.p0.h5       # mesh topology, written once at step 0 per rank
    ├── topology.p1.h5
    ├── gospl.0.p0.h5        # per-step, per-rank data file
    ├── gospl.0.p1.h5
    └── ...
```

Each rank writes its own partition of the data into `gospl.<step>.p<rank>.h5`. The XMF file (one per step, rank-0 only) tells Paraview which HDF5 dataset each visualisation attribute lives in. The XDMF file is a temporal collection of all the XMF files.

**Existing fields** (from `_outputMesh` at outmesh.py:203-399):

| HDF5 name | XMF name | Source | Conditional on |
|---|---|---|---|
| `elev` | `Z` | `self.hLocal.getArray()` (Vec) | always |
| `erodep` | `ED` | `self.cumEDLocal.getArray()` (Vec) | always |
| `EDrate` | `EDrate` | `self.EbLocal.getArray().copy()` (Vec, defensively copied) | always |
| `waterFill` | `waterFill` | `self.waterFilled` (numpy array) | `not self.fast` |
| `fillFA` | `fillFA` | `self.fillFAL.getArray().copy()` + floor + sea-cell clamp | always |
| `FA` | `FA` | `self.FAL.getArray().copy()` + floor + sea-cell clamp | always |
| `iceFA`, `iceH` | `iceFA`, `iceH` | `self.iceFAL/iceHL.getArray().copy()` + floor | `self.iceOn` |
| `flexIso` | `flexIso` | `self.localFlex` (numpy array) | `self.flexOn` |
| `soilH` | `soilH` | `self.Lsoil.getArray().copy()` (Vec) | `self.cptSoil` |
| `sedLoad` | `SL` | `self.vSedLocal.getArray().copy()` + floor | always |
| `sedLoadF` | `SLf` | `self.vSedFLocal.getArray().copy()` + floor (fine sub-flux; coarse = sedLoad − sedLoadF) | `self.stratLith` (dual lithology) |
| `uplift` | `vTec` | `self.upsub` (numpy array) | `self.upsub is not None` |
| `rain` | `Rain` | `self.rainVal` (numpy array) | `self.rainVal is not None` |
| `sea` | `sea` | constant function of `self.sealevel`, no HDF5 dataset | always (XMF only) |

Note the gap between **HDF5 dataset name** (lowercase verbose, e.g. `elev`, `erodep`, `EDrate`) and **XMF attribute name** (Paraview-facing, sometimes short or aliased: `Z`, `ED`, `SL`, `vTec`, `Rain`). When you add a field, pick both names and use each consistently.

---

## Files to touch (checklist)

```
gospl/tools/outmesh.py         — add field write in _outputMesh
                                  add XMF block in _save_DMPlex_XMF
                                  (optional) add restart-read block in readData
gospl/mesher/unstructuredmesh.py — destroy_DMPlex (ONLY if you allocate a new
                                  persistent Vec — see Step 2)
tests/test_regression.py       — (recommended) a slow-tier test asserting
                                  the field exists with correct shape
```

`unstructuredmesh.py` is HIGH-RISK per AGENTS.md. You only edit it if Step 2 forces a new persistent Vec.

---

## Step 1: Decide the source

Three cases. Pick the one that matches your field:

### Case A: source is already a persistent PETSc Vec
Examples: elevation (`self.hLocal`), cumulative erodep (`self.cumEDLocal`), erosion rate (`self.EbLocal`), flow accumulation (`self.FAL`), soil thickness (`self.Lsoil`).

**Recognise it by**: a `*Local` attribute exists, is destroyed in `destroy_DMPlex` (unstructuredmesh.py:746-779), and survives across timesteps.

Write directly: `f["myField"][:, 0] = self.myLocal.getArray()`. Use `.copy()` if you need to modify the array (apply a floor, clamp at sea cells, unit-convert) — otherwise the modification corrupts the underlying Vec.

### Case B: field is computed from existing fields
Examples: a unit conversion (FA → drainage area), a derived ratio (sed load / FA), a smoothed version of an existing field.

If the computation is **transient** (one expression, no further use), do it inline against a scratch Vec / numpy intermediate and write it. Reuse `self.tmp/tmpL/tmp1` — do NOT allocate a new Vec.

If the computation produces a value you need across **multiple methods** (not just inside the output writer), you need a new persistent Vec. Go to Step 2.

### Case C: source is already a numpy array
Examples: `self.upsub` (tectonic uplift), `self.rainVal` (precipitation), `self.localFlex` (flexure), `self.waterFilled` (water-filled topography).

**Recognise it by**: the attribute is set as `np.array(...)` or `np.zeros(...)`, not a `Vec().duplicate()`. The array is local (sized `self.lpoints` or `self.mpoints[locIDs]`).

Write directly: `f["myField"][:, 0] = self.myArray`. No `.copy()` needed unless you modify it. No `destroy_DMPlex` entry (numpy arrays are GC'd with the Model).

---

## Step 2: Allocate a scratch Vec (only if needed)

**Read AGENTS.md > Scratch vector contract first.** The existing scratch Vecs are reusable for transient computations:

| Vec | Allocated at | Local/Global |
|---|---|---|
| `self.tmp` | sedplex.py:30 | global |
| `self.tmpL` | sedplex.py:31 | local |
| `self.tmp1` | sedplex.py:32 | global |
| `self.Qs`, `self.QsL`, `self.nQs` | sedplex.py:33-35 | g, l, l |
| `self.h` | hillslope.py:39 | global |
| `self.hl` | hillslope.py:40 | local |
| `self.dh` | hillslope.py:41 | global |
| `self.newH` | SPL.py:41 | global (also held by SNES) |
| `self.stepED` | SPL.py:39 | global |
| `self.upsG`, `self.upsL` | outmesh.py:64-65 | global, local |

**Rule**: reuse `self.tmpL` for transient computations inside `_outputMesh`. The output writer runs only at output timesteps (not every simulation step), so by the time `_outputMesh` is called, the previous kernel's use of `tmpL` is irrelevant — you're free to overwrite it. After your use, the next caller will overwrite it again. This is exactly what AGENTS.md > Scratch vector contract says: *"Any kernel may overwrite them at any time. They are NOT persistent state."*

Allocate a **new named Vec only if** the field must survive across method calls — typically because another kernel writes to it during the timestep and `_outputMesh` later reads it. (None of the current output fields require this — see the existing patterns: every modification happens in the writer itself via `.copy()` on a Vec or a transient numpy variable.)

If you genuinely need a new persistent Vec (rare for outputs), allocate it in `WriteMesh.__init__` mirroring outmesh.py:64-66:

```python
# In WriteMesh.__init__
self.myFieldG = self.hGlobal.duplicate()    # global view
self.myFieldL = self.hLocal.duplicate()     # local view used for output
```

### CRITICAL: register the new Vec in destroy_DMPlex

Per AGENTS.md > High-risk modules: *"Adding a new persistent Vec elsewhere requires adding it to that destroy list or it leaks."*

Open `gospl/mesher/unstructuredmesh.py:739-779` and add your Vec to the explicit-list section near the other Local/Global pairs:

```python
def destroy_DMPlex(self):
    ...
    self.bL.destroy()
    self.bG.destroy()
    # NEW: persistent Vec for the myField output
    self.myFieldG.destroy()
    self.myFieldL.destroy()
    ...
```

If your Vec is feature-gated (e.g. only allocated when `self.iceOn`), gate the destroy too:

```python
if self.myFieldOn:
    self.myFieldG.destroy()
    self.myFieldL.destroy()
```

The lazy-solver loop at unstructuredmesh.py:789-800 is a SEPARATE mechanism for cached `self._ksp_*` / `self._snes_*` / `self._ts_*` objects. Forcing/output Vecs created at init go in the explicit list, not the loop.

---

## Step 3: Write to HDF5

Inside `_outputMesh` (outmesh.py:203-399), the HDF5 file is open as `f` in a `with h5py.File(h5file, "w") as f:` block running from line 251 to line 378. Add your dataset inside that block.

### Existing conventions — match them exactly

```python
f.create_dataset(
    "myField",                     # HDF5 dataset name (lowercase verbose)
    shape=(self.lpoints, 1),       # per-node, single component
    dtype="float32",               # 4 bytes — Paraview-friendly, halves disk
    compression="gzip",            # transparent gzip compression
)
f["myField"][:, 0] = <source>      # write the column-0 slice
```

**Conventions to match**:

| Convention | Value | Why |
|---|---|---|
| Shape | `(self.lpoints, 1)` for per-node scalars | matches existing reads in `readData` and the XMF Dimensions |
| Dtype | `float32` for output-only fields | matches existing reads; halves on-disk size |
| Dtype | `float64` if used in restart (compaction stratigraphy uses this) | full precision survives restart |
| Compression | `"gzip"` | adaptive, ~5-10× smaller for smooth fields |
| Chunking | (not specified) | h5py uses auto-chunking; do not override |
| Column slice | `[:, 0]` | shape is `(lpoints, 1)` because XMF Dimensions are `<N> 1` |

### The `.copy()` decision

Existing fields fall into three patterns:

```python
# Pattern 1 — write the Vec data directly (no modifications)
f["elev"][:, 0] = self.hLocal.getArray()                   # line 258
f["erodep"][:, 0] = self.cumEDLocal.getArray()             # line 265

# Pattern 2 — defensive .copy() (no modifications, but defensive)
data = self.EbLocal.getArray().copy()                      # line 272
f["EDrate"][:, 0] = data

# Pattern 3 — .copy() because we MUTATE the array (floor, clamp at sea)
data = self.FAL.getArray().copy()                          # line 299
data[data <= DISCHARGE_FLOOR] = DISCHARGE_FLOOR
if not self.fast:
    data[self.seaID] = 1.0
f["FA"][:, 0] = data

# Pattern 4 — write a numpy array directly (source is not a Vec)
f["flexIso"][:, 0] = self.localFlex                        # line 341
f["rain"][:, 0] = self.rainVal                             # line 375
```

Use Pattern 1 or 4 if you don't mutate. Use Pattern 3 if you apply a floor / sea-cell clamp / unit conversion before writing. Pattern 2 (defensive copy without modification) is leftover from an earlier convention — don't add new fields in this style.

### Conditional emission

Existing fields are gated by feature flags (line numbers from `_outputMesh`):

```python
if not self.fast:                # waterFill (L274), FA/fillFA sea clamp (L290, L301)
if self.iceOn:                   # iceFA, iceH (L305)
if self.flexOn:                  # flexIso (L334)
if self.cptSoil:                 # soilH (L342)
if self.upsub is not None:       # uplift (L360)
if self.rainVal is not None:     # rain (L368)
```

For a field that should always be written, no conditional. For a feature-gated field, follow these examples.

---

## Step 4: Update the XDMF manifest

`_save_DMPlex_XMF` (outmesh.py:495-684) writes one XMF file per step, listing all Paraview attributes. The function loops over `MPIsize` ranks (line 514) and writes one Block per rank. Inside each Block, each Attribute block follows the exact same shape. Add yours.

### Add an Attribute block

Inside the `for p in range(MPIsize):` loop, after the existing attribute blocks (or in feature-gated position to match your HDF5 conditional), insert:

```python
f.write('         <Attribute Type="Scalar" Center="Node" Name="MyField">\n')
f.write(
    '          <DataItem Format="HDF" NumberType="Float" Precision="4" '
)
f.write(
    'Dimensions="%d 1">%s:/myField</DataItem>\n' % (self.nodes[p], pfile)
)
f.write("         </Attribute>\n")
```

### XMF block conventions

| Field | Value | Why |
|---|---|---|
| `Type` | `"Scalar"` | one value per node. Use `"Vector"` only for 3-component fields (none exist today). |
| `Center` | `"Node"` | values at vertices, not at cell centers |
| `Name` | Paraview-facing name (e.g. `Z`, `ED`, `MyField`) | what shows up in Paraview's variable list |
| `Format` | `"HDF"` | external HDF5 backing |
| `NumberType` | `"Float"` | match the HDF5 dtype family |
| `Precision` | `"4"` | corresponds to `float32` in HDF5. Use `"8"` if you wrote `float64`. |
| `Dimensions` | `"%d 1" % self.nodes[p]` | MUST exactly match the HDF5 dataset's `shape` — see Common mistakes #3 |
| Data path | `%s:/myField` | `pfile` is the per-rank H5 file; the path is the HDF5 dataset name with leading slash |

### Feature-gating

If your HDF5 write is feature-gated, the XMF block must use the SAME gate. Otherwise Paraview will reference a non-existent dataset and silently corrupt the visualisation (or fail to load). See outmesh.py:611-619 for the `if self.flexOn:` pattern.

---

## Step 5: Register in the output loop

`_outputMesh` is the only place that calls `f.create_dataset` for per-step output. The function is one big sequential block (no helper functions per field), so your additions go between existing fields. Pick a logical neighbour:

- Floor-clamped flow-like fields → near `fillFA` / `FA` (L274-303)
- Ice fields → inside the `if self.iceOn:` block (L305-332)
- Tectonic / flexural fields → near `flexIso` / `uplift` (L334-341, L360-367)
- Forcing fields → near `rain` (L368-375)
- Pure new field → just before `if self.memclear:` at L376

Add the matching XMF block in `_save_DMPlex_XMF` in the SAME relative position so writers and readers stay in sync visually.

### Restart loop (optional)

`readData` (outmesh.py:401-493) restores the state Vecs from a previous run's HDF5 outputs. Only fields that participate in restarting the simulation appear there — `elev`, `erodep`, `EDrate`, `FA`, `fillFA`, `sedLoad`, optionally `sedLoadF` (dual lithology), `iceH`, `flexIso`, `soilH`.

- **If your field is purely a visualisation output**, do NOT add a read block in `readData`. Restart works fine without it.
- **If your field is state** (initialised at the start of a fresh run, expected to survive restarts), add a read block mirroring outmesh.py:434-441:

```python
# In readData, inside the `with h5py.File(h5file, "r") as hf:` block
if self.myFieldOn:
    if "/myField" in hf:                            # forward-compat: old runs
        self.myFieldL.setArray(np.array(hf["/myField"])[:, 0])
        self.dm.localToGlobal(self.myFieldL, self.myFieldG)
    else:
        self.myFieldL.set(0.)                       # default if absent
        self.dm.localToGlobal(self.myFieldL, self.myFieldG)
```

The `"/myField" in hf` guard makes restart files written by older code (without this field) still loadable — see the `stratK` precedent at outmesh.py:485-491.

---

## Step 6: Write a regression test (recommended)

Add to `tests/test_regression.py`:

```python
@pytest.mark.slow
def test_output_field_myField(minimal_model):
    """
    Protects: HOW_TO_ADD_OUTPUT.md > new myField field exists in HDF5.

    Silent failure prevented: a future refactor that drops the
    `_outputMesh` call to `f.create_dataset("myField", ...)` would
    silently produce an HDF5 file without the field. Paraview would
    then either skip the field (if the XMF block was also dropped)
    or fail to render (if XMF still references it).
    """
    import h5py
    import os

    model = minimal_model
    model.runProcesses()

    h5path = os.path.join(
        model.outputDir, "h5", f"gospl.0.p0.h5"
    )
    assert os.path.exists(h5path), f"output file missing: {h5path}"

    with h5py.File(h5path, "r") as hf:
        assert "myField" in hf, "myField dataset missing from HDF5 output"
        ds = hf["myField"]
        assert ds.shape == (model.lpoints, 1), \
            f"myField shape {ds.shape} != ({model.lpoints}, 1)"
        assert str(ds.dtype) == "float32", \
            f"myField dtype {ds.dtype} != float32"
```

Why `@pytest.mark.slow`: this needs a full `Model` instantiation and one `runProcesses()` call to populate the output directory. See `tests/conftest.py` for the `minimal_model` fixture.

The test runs at the fast tier when `pytest -m 'not slow'` is selected (i.e. it's skipped on PR CI), and at the full tier on nightly / pre-tag runs.

---

## Worked example: drainage basin area per node

Suppose we want a new output field `basinArea` = the area in m² that drains through each node, derived from `self.FAL` (flow accumulation, m³/yr) and a known rainfall rate. (For illustration; the actual physical conversion may need to involve `self.rainVal` rather than a constant.)

### Step 1: source

`self.FAL` is a persistent Local Vec (Case A). No new computation kernel needed — we'll do the unit conversion inline at output time on a `.copy()` of the array.

### Step 2: scratch Vec?

No new persistent Vec needed. The conversion is one numpy expression. Skip Step 2 entirely.

### Step 3: HDF5 write

In `_outputMesh`, add immediately after the existing `FA` block (around outmesh.py:303):

```python
f.create_dataset(
    "basinArea",
    shape=(self.lpoints, 1),
    dtype="float32",
    compression="gzip",
)
data = self.FAL.getArray().copy()
data[data < DISCHARGE_FLOOR] = DISCHARGE_FLOOR
# FA is m^3/yr; divide by representative rainfall (m/yr) → m^2 of basin
data = data / 1.0    # placeholder unit conversion; substitute real conversion
f["basinArea"][:, 0] = data
```

### Step 4: XMF manifest

In `_save_DMPlex_XMF`, add immediately after the existing `FA` Attribute block (around outmesh.py:581):

```python
f.write('         <Attribute Type="Scalar" Center="Node" Name="basinArea">\n')
f.write(
    '          <DataItem Format="HDF" NumberType="Float" Precision="4" '
)
f.write(
    'Dimensions="%d 1">%s:/basinArea</DataItem>\n' % (self.nodes[p], pfile)
)
f.write("         </Attribute>\n")
```

### Step 5: output loop registration

Already done — both blocks above sit inside the existing loops. No further wiring.

### Step 6: regression test

```python
@pytest.mark.slow
def test_output_basinArea(minimal_model):
    """Assert basinArea HDF5 dataset exists with correct shape."""
    import h5py, os
    model = minimal_model
    model.runProcesses()
    h5path = os.path.join(model.outputDir, "h5", "gospl.0.p0.h5")
    with h5py.File(h5path, "r") as hf:
        assert "basinArea" in hf
        assert hf["basinArea"].shape == (model.lpoints, 1)
        assert str(hf["basinArea"].dtype) == "float32"
        # Values should be strictly positive (we floored at DISCHARGE_FLOOR)
        from gospl.tools.constants import DISCHARGE_FLOOR
        import numpy as np
        assert (np.array(hf["basinArea"]) >= DISCHARGE_FLOOR).all()
```

That's it — no `destroy_DMPlex` change (Case A, no new Vec), no restart-loop change (visualisation only).

---

## Common mistakes

### 1. Using a scratch Vec and assuming it persists

```python
# WRONG — self.tmpL is scratch (AGENTS.md > Scratch vector contract).
# Another kernel will overwrite it before _outputMesh runs.
def _computeMyField(self):
    self.tmpL.setArray(<some expression>)        # later overwritten

def _outputMesh(self):
    ...
    f["myField"][:, 0] = self.tmpL.getArray()    # reads stale/wrong data
```

**Fix**: either compute inline inside `_outputMesh` (the recommended pattern — fields like `FA` do this), OR allocate a dedicated persistent Vec in `WriteMesh.__init__` and add it to `destroy_DMPlex`.

### 2. Forgetting destroy_DMPlex for a new persistent Vec

If Step 2 allocates a new `self.myFieldG = self.hGlobal.duplicate()` Vec, you MUST add `self.myFieldG.destroy()` to the explicit list in `unstructuredmesh.py:739-779`. The destroy list is hardcoded — there is no automatic scan. Omission leaks the PETSc object at simulation end, surfaces only as a warning under `petsc4py.PETSc.garbage_cleanup()` and as growing memory over many runs in HPC contexts.

The lazy-solver loop at outmesh.py:789-800 handles `self._ksp_*` / `self._snes_*` / `self._ts_*` only. Don't put Vecs in there.

### 3. XDMF dimension mismatch

The XMF Dimensions string MUST match the HDF5 dataset shape exactly. Paraview does not validate this — it silently reads the wrong number of bytes and renders garbage.

```python
# HDF5 dataset
f.create_dataset("myField", shape=(self.lpoints, 1), dtype="float32", ...)

# XMF block — Dimensions MUST be "<lpoints> 1", NOT "<lpoints>" alone
f.write(
    'Dimensions="%d 1">%s:/myField</DataItem>\n' % (self.nodes[p], pfile)
)
```

Other common dimension mismatches:
- Writing `shape=(self.lpoints, 3)` (vector) but using `Type="Scalar"` in XMF: silent garbage.
- Writing `dtype="float64"` but setting `Precision="4"` in XMF: silent half-precision.
- Per-rank Block referencing `self.nodes[p]` for the dataset size but writing `self.lpoints` to HDF5 — these must agree per partition (`self.nodes[p]` IS `lpoints` on rank `p`, set in outmesh.py:239).

If Paraview renders your field as garbled or wraps incorrectly across cells, suspect a dimension or dtype mismatch first.

---

## See also

- AGENTS.md > Scratch vector contract — when reuse is safe vs. when to allocate
- AGENTS.md > High-risk modules — destroy_DMPlex ownership of every persistent Vec
- [HOW_TO_ADD_FORCING.md](HOW_TO_ADD_FORCING.md) — sister runbook for the input side
- `tests/test_regression.py::test_mass_conservation` — closes the loop on the output path (reads `cumEDLocal` after a run)
