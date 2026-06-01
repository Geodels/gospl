# How to add a new forcing to goSPL

This runbook describes how to add a new time-varying forcing (e.g. precipitation, tectonic uplift, an erodibility multiplier) to goSPL. Following it mechanically reproduces the patterns already used by `raindata`, `sedfacdata`, `tedata`, and `tecdata`.

Read [`AGENTS.md`](../AGENTS.md) first — this document assumes you have read its **Forcing DataFrame layout contract** (column names are API, named `df.at[nb, col]` access), the **numpy ↔ PETSc boundary** rule, and the **High-risk modules** list.

---

## Overview

A *forcing* in goSPL is an external time-series of per-node values applied to the model at every timestep: precipitation rate, tectonic uplift, sediment-erodibility multiplier, elastic-thickness map, etc. Each forcing exists in **two zones** simultaneously per AGENTS.md "numpy ↔ PETSc boundary":

- **Numpy land** — a pandas `DataFrame` (`self.xxxdata`) built at parse time from the YAML input, plus the per-step numpy array of node values (`self.rainVal`, `self.sedFacMesh`, etc.).
- **PETSc land** — if the forcing feeds an implicit solve (like rain feeds `self.bL` for the flow-accumulation KSP), a halo-aware PETSc Vec.

A new forcing always touches the same six files. There are no exceptions.

---

## Files to touch (checklist)

```
gospl/tools/inputparser.py        — _defineXxx, _readXxx (YAML → DataFrame)
gospl/mesher/unstructuredmesh.py  — _updateXxx (per-step DataFrame → fields),
                                     applyForces (dispatcher), destroy_DMPlex
                                     (if you allocate any new Vec)
AGENTS.md                          — Forcing DataFrame layout contract table
tests/test_regression.py          — column-order regression guard
```

`unstructuredmesh.py` is **HIGH-RISK** per AGENTS.md > High-risk modules. Run the full regression suite (`pytest tests/`) after editing it.

---

## Step 1: Define the YAML schema

Decide on a top-level YAML key. Existing keys: `climate` (rain), `sedfactor`, `temap` (elastic thickness), `tectonics`. Pick a name that doesn't collide.

The minimum viable schema is a list of events, each event a dict with at least a `start` time:

```yaml
yourforcing:
  - start: 0.        # required, in years; events are sorted by start
    uniform: 1.0     # uniform branch — one of these is required
    # OR
    map: [path/to/file, fieldname]   # map branch — NPZ archive + key name
```

Read the YAML section in `inputparser.py` using the **outer-section try/except** pattern (this is one of the AGENTS.md > YAML parsing helpers "outer-section blocks" that you do NOT collapse with `_get_param`):

```python
def _readYourForcing(self):
    """Parse yourforcing conditions."""
    yfData = None
    try:
        yfDict = self.input["yourforcing"]
        yfSort = sorted(yfDict, key=itemgetter("start"))
        for k in range(len(yfSort)):
            yStart = None
            yUniform = None
            yMap = None
            # start is required — the print+raise pattern is intentional
            try:
                yStart = yfSort[k]["start"]
            except KeyError:
                print(
                    "For each yourforcing event a start time is required.",
                    flush=True,
                )
                raise ValueError(
                    "Yourforcing event {} has no parameter start".format(k)
                )
            yUniform = yfSort[k].get("uniform")
            yMap = yfSort[k].get("map")

            if yMap is not None:
                # ... open and validate the NPZ file, same shape as
                # _readErofactor at inputparser.py:887-1003 ...
                pass

            yfData = self._defineYourForcing(k, yStart, yMap, yUniform, yfData)

        # Prepend a row at tStart if the first event starts later
        if yfData["start"][0] > self.tStart:
            tmpYf = []
            tmpYf.insert(
                0, {"start": self.tStart, "yUni": <sentinel-default>,
                    "yMap": None, "yKey": None}
            )
            yfData = pd.concat(
                [pd.DataFrame(tmpYf, columns=["start", "yUni", "yMap", "yKey"]),
                 yfData],
                ignore_index=True,
            )
        self.yfdata = yfData.copy()
        self.yfdata.reset_index(drop=True, inplace=True)
        self.yfNb = len(self.yfdata)

    except KeyError:
        # Section missing → forcing is off. Set self.yfdata = None so
        # applyForces' dispatcher can gate the per-step update on it.
        self.yfdata = None
        self.yfNb = 0
```

Then add `self._readYourForcing()` to the call list inside `ReadYaml.__init__` at the bottom of the existing `self._readX()` block.

### Key rule from AGENTS.md > YAML parsing helpers
- Inside the function body, use `dict.get(key, default)` on a pre-extracted section dict (like `yfSort[k].get("uniform")`).
- Use `self._get_param("section", "key", default=X)` only when reading `self.input` directly without pre-extracting.
- Do NOT use bare `try/except KeyError` for default-on-miss outside the outer-section block and the "required key with raise" pattern shown above.

---

## Step 2: Build the forcing DataFrame

Define the inner builder method `_defineYourForcing` immediately above `_readYourForcing`. The reference is `_defineErofactor` at inputparser.py:832 — the minimum viable pattern.

```python
def _defineYourForcing(self, k, yStart, yMap, yUniform, yfData):
    """
    Define yourforcing conditions for event k.

    :arg k: event number (0-based)
    :arg yStart: event start time (years)
    :arg yMap: [npz_path, field_name] for the spatial-map branch, else None
    :arg yUniform: scalar value for the uniform branch, else None
    :arg yfData: accumulator DataFrame (None on the first call)
    :return: yfData with event k appended
    """
    if yMap is None and yUniform is None:
        print(
            "For each yourforcing event a value (uniform) or a map (map) "
            "is required.",
            flush=True,
        )
        raise ValueError(
            "Yourforcing event {} has no value (uniform) or map (map)."
            .format(k)
        )

    tmpRow = []
    if yMap is None:
        tmpRow.insert(
            0,
            {"start": yStart, "yUni": yUniform, "yMap": None, "yKey": None},
        )
    else:
        tmpRow.insert(
            0,
            {"start": yStart, "yUni": None,
             "yMap": yMap[0] + ".npz", "yKey": yMap[1]},
        )

    # IMPORTANT: column-list order here defines the DataFrame's column-name set,
    # which is the API contract. Consumers use df.at[nb, "yUni"] not iloc.
    columns = ["start", "yUni", "yMap", "yKey"]
    if k == 0:
        yfData = pd.DataFrame(tmpRow, columns=columns)
    else:
        yfData = pd.concat(
            [yfData, pd.DataFrame(tmpRow, columns=columns)],
            ignore_index=True,
        )
    return yfData
```

### Column-order rule (AGENTS.md > Forcing DataFrame layout contract)
- Pick column names with a short prefix that matches your forcing (rain → `r*`, sedfactor → `s*`, te → `t*`, tectonics → no prefix for `start`/`end`). For erodibility it'd be `e*`: `eUni`, `eMap`, `eKey`.
- **Use the SAME column list everywhere it appears** — the `pd.DataFrame(... columns=[...])` call at `k == 0`, the `pd.concat(... columns=[...])` call at `k > 0`, and the tStart-prepend block in `_readYourForcing`. A mismatch silently inserts NaN at runtime.
- Once consumers use `df.at[nb, "col_name"]` access, you can add new columns in any order without breaking them. Just keep the existing column names stable.

### The 2026-06 rUni/sUni bug — read this carefully
The original `_defineErofactor` had a typo: dict key `"rUni"` was passed into a DataFrame built with column `"sUni"`. The column ended up as NaN at runtime, and `_updateEroFactor` crashed on `np.load(None)`. The bug shape:

```python
# BUG — dict key does not match DataFrame column name
tmpRow.insert(0, {"start": sStart, "rUni": sUniform, ...})    # rUni ← WRONG
yfData = pd.DataFrame(tmpRow, columns=["start", "sUni", ...])  # sUni
```

The lesson: **the dict-key names in the row literal MUST exactly match the column names in `pd.DataFrame(... columns=[...])`**. The regression test in Step 6 will catch this if you make the same mistake.

---

## Step 3: Register in `applyForces`

`applyForces` is the per-step forcing dispatcher in `unstructuredmesh.py:578`. It runs at **init** (`model.py:181`) AND every timestep (`model.py:274`). Add your update call near the existing rain / sedfactor / tectonics calls. The order matters only if your forcing depends on another forcing's output for the same step — which is usually NOT the case.

```python
def applyForces(self):
    """
    Finds the different values for climatic, tectonic and sea-level forcing
    that will be applied at any given time interval during the simulation.
    """
    t0 = process_time()
    # Sea level ...
    # Climate information
    if self.oroOn:
        self.cptOrography()
    else:
        self._updateRain()

    # Erodibility factor information
    if self.sedfacdata is not None:
        self._updateEroFactor()

    # NEW: your forcing here. Gate on `self.yfdata is not None`
    # so missing-section users aren't penalised.
    if self.yfdata is not None:
        self._updateYourForcing()

    # ...rest of applyForces (timing print, applyTectonics, boundaries)
```

The init-vs-per-step distinction is **automatic**: `applyForces` is called at both. Your `_updateXxx` method must handle the `self.tNow == self.tStart` (init) case correctly — typically by treating it identically to a per-step call. The reference `_updateRain` does this implicitly: `self.rainNb == -1` (set in `_buildMesh` at unstructuredmesh.py:481-484) triggers the "load event 0" branch on the very first call.

---

## Step 4: Write the per-timestep update method

Define `_updateYourForcing` near `_updateRain` (unstructuredmesh.py:656) and `_updateEroFactor` (unstructuredmesh.py:706). The pattern is the **time-window advance** — fetch the next event when the current simulation time passes its `start`, otherwise reuse the cached event:

```python
def _updateYourForcing(self):
    """
    Finds current yourforcing values for the considered time interval.
    """
    nb = self.yfNb
    if nb < len(self.yfdata) - 1:
        # NOTE: <= for non-tectonics forcings. Tectonics uses
        # `< self.tNow + self.dt` instead (strict, looking one dt ahead) —
        # see tectonics.py:64 and AGENTS.md commentary on the
        # forcing time-window operator.
        if self.yfdata.at[nb + 1, "start"] <= self.tNow:
            nb += 1

    if nb > self.yfNb or nb == -1:
        if nb == -1:
            nb = 0
        self.yfNb = nb

        if pd.isnull(self.yfdata["yUni"][nb]):
            # Map branch: load the NPZ archive and extract the named field.
            loadData = np.load(self.yfdata.at[nb, "yMap"])
            yVal = loadData[self.yfdata.at[nb, "yKey"]]
            del loadData
        else:
            # Uniform branch: broadcast the scalar across all mesh nodes.
            yVal = np.full(self.mpoints, self.yfdata.at[nb, "yUni"])

        # Optional: clip to physical bounds, store globally for restart
        yVal[yVal < 0] = 0.0
        self.yfMesh = yVal

    # Slice to the local partition of this rank. Always do this even when
    # nb didn't change — keeps self.yfVal in sync with self.locIDs.
    self.yfVal = self.yfMesh[self.locIDs]
    return
```

### Named access (AGENTS.md > Forcing DataFrame layout contract)
- `self.yfdata.at[nb, "yUni"]` — single-scalar lookup by row & column name.
- Do NOT use `self.yfdata.iloc[nb, 1]` — position-based access makes column order load-bearing. The 30-site refactor in 2026-06 eliminated all of these; the regression test `test_forcing_column_order` exists specifically to catch reintroductions.

### Forcing time-window operator (intentional difference)
Three forcings use `<= self.tNow` (rain, sedfactor, te — see unstructuredmesh.py:668, 718; addprocess.py:153). One forcing uses `< self.tNow + self.dt` (tectonics — see tectonics.py:64). The difference is intentional:
- `<=` advances when the event's start time has been reached.
- `< tNow + dt` advances when the event's start falls within the upcoming step — applies the new tectonic field a step earlier, which the velocity-driven advection needs.

Use `<= self.tNow` unless you have a specific reason matching tectonics' need.

---

## Step 5: Allocate and register any new Vecs

If your forcing only feeds existing scratch arrays (like `self.bL` for rain), you don't need new Vecs.

If you need a **new persistent Vec** (one that survives across timesteps and is used in a PETSc solve), allocate it in the relevant mixin's `__init__` — typically `UnstMesh.__init__` via the `_buildMesh` chain, or in the consumer mixin if local to one solver. Use `.duplicate()` from an existing Vec to inherit the DMPlex layout:

```python
# Example: a persistent global Vec sized like hGlobal
self.yfFlux = self.hGlobal.duplicate()
self.yfFluxL = self.hLocal.duplicate()
```

### CRITICAL: register in `destroy_DMPlex`
Per AGENTS.md > High-risk modules: *"Adding a new persistent Vec elsewhere requires adding it to that destroy list or it leaks."*

Open `unstructuredmesh.py:739-799`. The function is a flat list of `obj.destroy()` calls plus a loop for cached solvers. Add yours to the explicit-list section near the other Vecs:

```python
def destroy_DMPlex(self):
    ...
    self.bL.destroy()
    self.bG.destroy()
    # NEW: persistent forcing Vecs you allocated
    self.yfFlux.destroy()
    self.yfFluxL.destroy()
    ...
```

If your Vec is only allocated under a feature flag (like `self.iceHL` is only created when `self.iceOn`), gate the destroy too:

```python
if self.yfOn:
    self.yfFlux.destroy()
    self.yfFluxL.destroy()
```

The 2026-06 cached-KSP/SNES/TS attributes are destroyed by the second loop (`for name in (...): ... obj.destroy()`). That's a SEPARATE mechanism — only use it for `self._ksp_*`, `self._snes_*`, `self._ts_*` and their helper Vecs created lazily inside solver methods. Forcing Vecs created at init go in the explicit list.

---

## Step 6: Write the regression test

`tests/test_regression.py::test_forcing_column_order` (around line 116) protects the column-name contract. **Add a block for your new forcing** following the existing four:

```python
# ---- yfdata: built via _readYourForcing (uniform-only path) ----
parser = _bare_parser()
parser.input = {"yourforcing": [{"start": 0.0, "uniform": 1.0}]}
parser._readYourForcing()
assert list(parser.yfdata.columns) == [
    "start", "yUni", "yMap", "yKey",
], (
    "yfdata column order drift. AGENTS.md says: "
    "start, yUni, yMap, yKey."
)
```

This block:
- Uses the `_bare_parser()` fixture (defined at the top of `test_regression.py`) which bypasses `ReadYaml.__init__` and only injects the attributes needed by the forcing parsers.
- Runs in the fast tier — no real mesh, no PETSc DMPlex. The test should complete in milliseconds.
- Catches three classes of regression: column-name typo (e.g. the rUni/sUni bug), column-order change, removal of a column the contract guarantees.

Also update **AGENTS.md > Forcing DataFrame layout contract** to add a row for your DataFrame:

```markdown
| `self.yfdata` | `_defineYourForcing` (inputparser.py) | `start, yUni, yMap, yKey` | unstructuredmesh.py |
```

---

## Worked example: spatially varying erodibility

A concrete walk-through. Suppose we want a new forcing `erodibility` that lets the user impose a per-node erodibility multiplier (similar to `sedfactor` but conceptually distinct — say, derived from a NetCDF dataset that the user has pre-converted to NPZ).

### YAML

```yaml
erodibility:
  - start: 0.
    uniform: 1.0
  - start: 50000.
    map: [data/geology_50ka, K_multiplier]
  - start: 100000.
    uniform: 0.5
```

### Step 1+2: parser

In `gospl/tools/inputparser.py`, add `_defineErodibility` (above) and `_readErodibility` (below). Call `self._readErodibility()` from `ReadYaml.__init__`.

```python
def _defineErodibility(self, k, eStart, eMap, eUniform, eData):
    if eMap is None and eUniform is None:
        print(
            "For each erodibility event a value (uniform) or a map (map) "
            "is required.",
            flush=True,
        )
        raise ValueError(
            "Erodibility event {} has no value (uniform) or map (map)."
            .format(k)
        )
    tmpRow = []
    if eMap is None:
        tmpRow.insert(
            0,
            {"start": eStart, "eUni": eUniform, "eMap": None, "eKey": None},
        )
    else:
        tmpRow.insert(
            0,
            {"start": eStart, "eUni": None,
             "eMap": eMap[0] + ".npz", "eKey": eMap[1]},
        )
    columns = ["start", "eUni", "eMap", "eKey"]
    if k == 0:
        eData = pd.DataFrame(tmpRow, columns=columns)
    else:
        eData = pd.concat(
            [eData, pd.DataFrame(tmpRow, columns=columns)],
            ignore_index=True,
        )
    return eData

def _readErodibility(self):
    """Parse erodibility forcing conditions."""
    eData = None
    try:
        eDict = self.input["erodibility"]
        eSort = sorted(eDict, key=itemgetter("start"))
        for k in range(len(eSort)):
            eStart = None
            eUniform = None
            eMap = None
            try:
                eStart = eSort[k]["start"]
            except KeyError:
                print(
                    "For each erodibility event a start time is required.",
                    flush=True,
                )
                raise ValueError(
                    "Erodibility event {} has no parameter start".format(k)
                )
            eUniform = eSort[k].get("uniform")
            eMap = eSort[k].get("map")
            # ... (file-exists check + NPZ-field-exists check,
            # mirroring _readErofactor at inputparser.py:914-944) ...
            eData = self._defineErodibility(k, eStart, eMap, eUniform, eData)

        if eData["start"][0] > self.tStart:
            tmpE = []
            tmpE.insert(
                0, {"start": self.tStart, "eUni": 1.0, "eMap": None,
                    "eKey": None}
            )
            eData = pd.concat(
                [pd.DataFrame(tmpE,
                              columns=["start", "eUni", "eMap", "eKey"]),
                 eData],
                ignore_index=True,
            )
        self.eData = eData.copy()
        self.eData.reset_index(drop=True, inplace=True)
        self.eNb = len(self.eData)
    except KeyError:
        self.eData = None
        self.eNb = 0
```

### Step 3+4: dispatcher and update method

In `gospl/mesher/unstructuredmesh.py`, near `_updateEroFactor`:

```python
def _updateErodibility(self):
    """Finds current erodibility multiplier values for the considered time interval."""
    nb = self.eNb
    if nb < len(self.eData) - 1:
        if self.eData.at[nb + 1, "start"] <= self.tNow:
            nb += 1

    if nb > self.eNb or nb == -1:
        if nb == -1:
            nb = 0
        self.eNb = nb

        if pd.isnull(self.eData["eUni"][nb]):
            loadData = np.load(self.eData.at[nb, "eMap"])
            eVal = loadData[self.eData.at[nb, "eKey"]]
            del loadData
        else:
            eVal = np.full(self.mpoints, self.eData.at[nb, "eUni"])

        # Floor at small positive to avoid divide-by-zero in downstream SPL.
        eVal[eVal < 0.01] = 0.01
        self.eroMesh = eVal

    self.eroFactor = self.eroMesh[self.locIDs]
    return
```

And in `applyForces`:

```python
# ... after the existing _updateEroFactor block ...
if self.eData is not None:
    self._updateErodibility()
```

The SPL kernels can then read `self.eroFactor` as a per-node multiplier (similar to how the existing `self.sedfacVal` is used in `eroder/SPL.py`).

### Step 5: no new Vecs

`self.eroFactor` is a per-rank numpy array (lpoints,) — no PETSc Vec, no entry in `destroy_DMPlex`. The `self.eroMesh` global array is GC'd with the Model. Skip Step 5.

### Step 6: regression test

In `tests/test_regression.py::test_forcing_column_order`, add:

```python
# ---- eData: built via _readErodibility (uniform-only path) ----
parser = _bare_parser()
parser.input = {"erodibility": [{"start": 0.0, "uniform": 1.0}]}
parser._readErodibility()
assert list(parser.eData.columns) == [
    "start", "eUni", "eMap", "eKey",
], (
    "eData column order drift. AGENTS.md says: "
    "start, eUni, eMap, eKey."
)
```

And in `AGENTS.md > Forcing DataFrame layout contract`:

```markdown
| `self.eData` | `_defineErodibility` (inputparser.py) | `start, eUni, eMap, eKey` | unstructuredmesh.py |
```

Run `pytest tests/` — all six existing tests plus your new column-order assertion should pass.

---

## Common mistakes

1. **Dict key ≠ DataFrame column name** (the 2026-06 rUni/sUni bug regressed). The row literal you pass into `pd.DataFrame(... columns=[...])` must use the SAME spelling for its keys as the column-list strings. A mismatch silently produces a NaN column; the consumer then crashes on `np.load(None)` or similar. The regression test in Step 6 catches this — but only if you remember to add the test.

2. **Using `df.iloc[nb, k]` instead of `df.at[nb, "col"]`** in the consumer. This makes column ORDER load-bearing (any future column-append silently shifts every consumer). 30 such sites were removed in 2026-06; the AGENTS.md "Forcing DataFrame layout contract" forbids reintroducing the pattern.

3. **Forgetting to add new persistent Vecs to `destroy_DMPlex`** (unstructuredmesh.py:739-799). Cited in AGENTS.md > High-risk modules: the destroy list is hardcoded, and forgetting your Vec leaks the PETSc object at simulation end. Scratch Vecs created locally inside a method and explicitly destroyed at end-of-method don't go in the list; persistent ones do. Cached lazy solvers (`self._ksp_*`, `self._snes_*`, `self._ts_*`) use the *separate* lazy-solver loop at the bottom of `destroy_DMPlex`, not the explicit list.

4. **Forgetting to gate the per-step update on `self.xxxdata is not None`** in `applyForces`. If the user omits your YAML section, the parser's outer-section except sets `self.xxxdata = None` (AGENTS.md > YAML parsing helpers > "outer-section blocks"). The dispatcher must check before calling the update method, otherwise users without your forcing get an `AttributeError` on every step.

5. **Using `<` instead of `<=` (or vice versa) in the time-window advance** in your `_updateXxx`. Tectonics uses `< self.tNow + self.dt` (strict, +dt offset); rain/sedfactor/te use `<= self.tNow` (non-strict, current step). They are not interchangeable — they shift the activation of a new event by up to one timestep. Use `<= self.tNow` unless your forcing needs to apply during the *next* step (the way tectonics does for advection).
