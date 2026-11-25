"""
Build the reference pages of the user guide from their data files.

Called by ``docs/conf.py`` at the start of every build. It reads, in
``docs/user_guide/reference/``:

- ``inputs.json``: every input key, by block (unit, default, meaning,
  typical range, dependence on the resolution);
- ``outputs.json``: every output field, by group (unit, meaning, typical
  range, dependence on the resolution, sign);
- ``modules.json``: the optional modules, what each requires or excludes;

and writes ``inputs.rst``, ``outputs_ref.rst`` and ``modules.rst`` next to
them (not tracked by git). Edit the JSON files, not the generated pages:
when a process gains a key or an output, add it there (see AGENTS.md,
checklist item 14).
"""

from __future__ import annotations

import json
import os
import re

HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "reference")

ROLES = """.. role:: refunit
   :class: gospl-ref-unit

.. role:: refdefault
   :class: gospl-ref-default

"""

# sections of the input page: the top-level blocks of each theme
THEMES = [
    ("Mesh, time and outputs", ["name", "domain", "time", "output"]),
    ("Climate and sea level", ["climate", "sea", "orography", "sst"]),
    ("Rivers and hillslopes", ["spl", "diffusion", "sedfactor"]),
    ("Soil", ["soil"]),
    ("Ice", ["ice"]),
    ("Groundwater, duricrust and geochemistry", ["groundwater"]),
    ("Stratigraphy and provenance", ["compaction", "strata", "provenance"]),
    ("Marine systems", ["gravity", "carbonate", "pelagic", "waves", "coastal_erosion",
                        "carbonate_supply", "karst"]),
    ("Tectonics and flexure", ["tectonics", "flexure", "temap"]),
]

FILTER = """.. raw:: html

   <input type="search" class="gospl-ref-filter"
          placeholder="Filter: type a key, a unit or a word (e.g. dt, m/yr, porosity)"
          aria-label="Filter the entries">

"""

SCRIPT = """.. raw:: html

   <script src="../../_static/modchart.js"></script>

"""


def _rst(text):
    """Plain text to safe inline RST: `x` becomes ``x``, special characters
    that would start markup are escaped."""
    if not text:
        return ""
    text = str(text).replace("\\", "\\\\")
    parts = re.split(r"(``[^`]+``|`[^`]+`)", text)
    out = []
    for p in parts:
        if p.startswith("``"):
            out.append(p)
        elif p.startswith("`") and p.endswith("`"):
            out.append("``%s``" % p[1:-1])
        else:
            p = re.sub(r"([*|])", r"\\\1", p)
            p = re.sub(r"(\w)_(?=\W|$)", r"\1\\_", p)
            out.append(p)
    return "".join(out)


def _anchor(path):
    """Stable anchor of an input block on the input reference page."""
    return "in-" + re.sub(r"[^a-z0-9]+", "-", path.lower()).strip("-")


def _title(text, char):
    return "%s\n%s\n\n" % (text, char * max(len(text), 3))


def _entry(name, unit, default, rows):
    head = "``%s``" % name
    if unit:
        head += " :refunit:`%s`" % _rst(unit).replace("`", "")
    if default not in (None, ""):
        d = str(default)
        head += " :refdefault:`%s`" % _rst(
            d if d.lower().startswith("required") else "default %s" % d).replace("`", "")
    L = [".. container:: gospl-ref-entry", "", "   .. container:: gospl-ref-head", "",
         "      %s" % head, ""]
    for label, text in rows:
        if not text or str(text).strip().lower() in ("", "none", "n/a", "-"):
            continue
        if label == "Meaning":
            L += ["   %s" % _rst(text), ""]
        else:
            L += ["   :gospl-ref-label:`%s` %s" % (label, _rst(text)), ""]
    return "\n".join(L) + "\n"


def _inputs(blocks):
    s = [".. _input_reference:\n\n", "=" * 34 + "\n", "Input parameters: units and ranges\n",
         "=" * 34 + "\n\n", ROLES,
         ".. role:: gospl-ref-label\n   :class: gospl-ref-label\n\n",
         "Every key of the YAML input file, block by block: its unit, its default, what\n"
         "it does, the values used in practice and how a sensible value changes with\n"
         "the mesh spacing and the time step. The keys themselves, with examples of\n"
         "each block, are explained in :ref:`inputfile` and the pages that follow it.\n"
         "Units use ``yr`` for years and ``m`` for metres; ``-`` marks a\n"
         "dimensionless value. Which blocks need which: :ref:`module_chooser`.\n\n",
         FILTER, SCRIPT,
         ".. contents:: On this page\n   :local:\n   :depth: 1\n\n\n"]
    groups, order = {}, []
    for b in blocks:
        top = re.split(r"[.\[]", b["path"])[0]
        theme = next((t for t, tops in THEMES if top in tops), "Other blocks")
        if theme not in groups:
            groups[theme] = []
            order.append(theme)
        groups[theme].append(b)
    order.sort(key=lambda t: [x for x, _ in THEMES].index(t) if t in dict(THEMES) else 99)
    for theme in order:
        s.append(_title(theme, "-"))
        for b in groups[theme]:
            s.append(".. _%s:\n\n" % _anchor(b["path"]))
            s.append(_title("``%s``" % b["path"], "^"))
            if b.get("summary"):
                s.append(_rst(b["summary"]) + "\n\n")
            for k in b.get("keys", []):
                s.append(_entry(k["key"], k.get("unit"), k.get("default"),
                                [("Meaning", k.get("meaning")), ("Range:", k.get("range")),
                                 ("Resolution:", k.get("resolution")), ("Note:", k.get("notes"))]))
                s.append("\n")
    return "".join(s)


def _outputs(groups):
    s = [".. _output_reference:\n\n", "=" * 34 + "\n", "Output fields: units and ranges\n",
         "=" * 34 + "\n\n", ROLES,
         ".. role:: gospl-ref-label\n   :class: gospl-ref-label\n\n",
         "Every field goSPL writes, by process: its unit, what it means, the values\n"
         "found in practice and how they depend on the mesh spacing and the time step.\n"
         "How to open and plot the outputs: :ref:`outputs`.\n\n",
         FILTER, SCRIPT,
         ".. contents:: On this page\n   :local:\n   :depth: 1\n\n\n"]
    for g in groups:
        s.append(_title(g["group"], "-"))
        if g.get("condition") and g["condition"].strip().lower() != "always":
            s.append("Written %s.\n\n" % _rst(g["condition"]).rstrip("."))
        for f in g.get("fields", []):
            where = {"surface": "", "strata": "stratigraphy file (per layer)",
                     "restart": "restart only"}.get(f.get("file", ""), f.get("file", ""))
            s.append(_entry(f["name"], f.get("unit"), "",
                            [("Meaning", f.get("meaning")), ("Sign:", f.get("sign")),
                             ("Range:", f.get("range")), ("Resolution:", f.get("resolution")),
                             ("File:", where)]))
            s.append("\n")
    return "".join(s)


def _strip_refs(obj):
    """Drop the code citations (``file.py:123``) kept in the data for the
    maintainers: the page is for users."""
    if isinstance(obj, dict):
        return {k: _strip_refs(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_strip_refs(v) for v in obj]
    if isinstance(obj, str):
        ref = r"(?:docs )?[\w/.]+\.(?:py|rst|md):\d+(?:-\d+)?"
        obj = re.sub(r"\s*\(\s*%s(?:[,;]\s*%s)*\s*(?::[^)]*)?\)" % (ref, ref), "", obj)
        obj = re.sub(r"[,;]?\s*%s(?::[^,;)]*)?" % ref, "", obj)
        obj = re.sub(r"\(\s*,\s*", "(", obj)
        obj = re.sub(r"\s*\(\s*\d+(?:\s*,\s*\d+)+\s*\)", "", obj)
        obj = re.sub(r" {2,}", " ", obj)
        obj = re.sub(r"\(\s*", "(", re.sub(r"\s*\)", ")", obj))
        return re.sub(r"\s+([.,;])", r"\1", obj).strip()
    return obj


# input-reference entries shown with each module of the chooser:
# (input path, keys or None for the whole block)
MOD_INPUTS = {
    "core": [("domain", ["npdata", "flowdir", "bc"]), ("time", ["start", "end", "dt", "tout"]),
             ("spl", ["K", "m", "d"])],
    "sea": [("sea", None)],
    "climate": [("climate[]", ["start", "uniform", "map"])],
    "evaporation": [("climate[]", ["evap_uniform", "evap_map"])],
    "orography": [("orography", None)],
    "sedfactor": [("sedfactor[]", None)],
    "spl_transport": [("spl", ["G"])],
    "nlspl": [("spl", ["n"])],
    "soil": [("soil", None)],
    "deposition": [("domain", ["nodep"])],
    "hillslope": [("diffusion", ["hillslopeKa", "hillslopeKm", "hillslopeSc", "hillslopeNb",
                                 "hillslopenl"])],
    "channelQ": [("diffusion", ["channelQ"])],
    "marine_deposition": [("domain", ["seadepo"]),
                          ("diffusion", ["nonlinKm", "clinSlp", "marineSolver", "picardSub"])],
    "gravity": [("gravity", None), ("gravity.flood", None)],
    "waves": [("waves", None), ("waves.climate[]", None)],
    "coastal_erosion": [("coastal_erosion", None)],
    "strata": [("time", ["strat"]), ("compaction", None)],
    "dual": [("strata", None), ("strata.coarse", None), ("strata.fine", None)],
    "lithologies": [("strata.lithologies[]", None)],
    "provenance": [("provenance", None)],
    "sst": [("sst", None), ("sst.events[]", None)],
    "carbonate": [("carbonate", None), ("carbonate.factories[]", None),
                  ("carbonate.factories[].response[]", None)],
    "pelagic": [("pelagic[]", None)],
    "carbonate_supply": [("carbonate_supply", None)],
    "karst": [("karst", None)],
    "groundwater": [("groundwater", None)],
    "duricrust": [("groundwater.duricrust", None)],
    "geochem": [("groundwater.geochem", None), ("groundwater.geochem.species[]", None)],
    "river_load": [("groundwater.geochem", ["river_load"])],
    "marine_coupling": [("groundwater.geochem", ["marine_coupling"])],
    "ice": [("ice", None), ("ice.abrasion", None), ("ice.till", None)],
    "flexure_flat": [("flexure", None)],
    "flexure_global": [("flexure", None)],
    "temap": [("temap[]", None)],
    "tectonics_vertical": [("tectonics[]", None)],
    "advection": [("domain", ["advect"]), ("tectonics[]", None)],
    "plate_interp": [("domain", ["advect"]), ("tectonics[]", None)],
    "interp_remap": [("domain", ["interp_remap"])],
    "cyclic": [("domain", ["bc"])],
}

# where to get a mesh, shown with the skeleton input file
MESH_LINKS = {
    "global_download": "https://www.dropbox.com/scl/fo/kior7a0kf9ylqxnnn7vk5/"
                       "AJfSqKiefxvb4fHik_uRfeQ?rlkey=mh9hpd1z0i7u8qwlgspdo2zdf&st=8c3i8fir&dl=0",
    "global_build": "examples/generated/continental_flux/build_inputs",
    "global_data": "examples/generated/global_mesh_data/index",
    "regional_build": "examples/generated/generate_2Dmesh/index",
}

# technical guide page (docname) explaining each module
MOD_TECH = {
    "core": "tech_guide/flow", "climate": "tech_guide/flow", "evaporation": "tech_guide/flow",
    "sedfactor": "tech_guide/ero", "spl_transport": "tech_guide/ero", "nlspl": "tech_guide/ero",
    "soil": "tech_guide/ero", "deposition": "tech_guide/dep", "hillslope": "tech_guide/hill",
    "channelQ": "tech_guide/hill", "marine_deposition": "tech_guide/hill",
    "gravity": "tech_guide/gravity", "waves": "tech_guide/carb", "coastal_erosion": "tech_guide/carb",
    "strata": "tech_guide/strat", "dual": "tech_guide/strat", "lithologies": "tech_guide/strat",
    "provenance": "tech_guide/provenance", "sst": "tech_guide/carb", "carbonate": "tech_guide/carb",
    "pelagic": "tech_guide/carb", "carbonate_supply": "tech_guide/carb", "karst": "tech_guide/carb",
    "groundwater": "tech_guide/groundwater", "duricrust": "tech_guide/groundwater",
    "geochem": "tech_guide/groundwater", "river_load": "tech_guide/groundwater",
    "marine_coupling": "tech_guide/groundwater", "ice": "tech_guide/ice",
    "flexure_flat": "tech_guide/tecto", "flexure_global": "tech_guide/tecto", "temap": "tech_guide/tecto",
    "tectonics_vertical": "tech_guide/tecto", "advection": "tech_guide/tecto",
    "plate_interp": "tech_guide/tecto", "interp_remap": "tech_guide/tecto", "cyclic": "tech_guide/tecto",
}

# minimal input of each module for the skeleton YAML of the chooser:
# (block path, key, value, input-reference path or None); a value starting
# with "\n" is a YAML list written as is under the key. "{end}", "{dt}",
# "{tout}", "{strat}" are filled from the resource form.
MOD_YAML = {
    "core": [("", "name", "'my goSPL run'", None),
             ("domain", "npdata", "['mesh', 'v', 'c', 'z']", "domain"),
             ("domain", "flowdir", "5", "domain"),
             ("domain", "bc", "'oooo'", "domain"),
             ("time", "start", "0.", "time"), ("time", "end", "{end}", "time"),
             ("time", "dt", "{dt}", "time"), ("time", "tout", "{tout}", "time"),
             ("output", "dir", "'output'", "output"), ("output", "makedir", "True", "output"),
             ("spl", "K", "1.e-6", "spl"), ("spl", "d", "0.", "spl"), ("spl", "m", "0.5", "spl")],
    "sea": [("sea", "position", "0.", "sea")],
    "climate": [("", "climate", "\n  - start: 0.\n    uniform: 1.   # rainfall (m/yr), or map: ['rain', 'r']", "climate[]")],
    "evaporation": [],
    "orography": [("orography", "wind_speed", "10.", "orography"), ("orography", "wind_dir", "270.", "orography"),
                  ("orography", "conv_time", "1000.", "orography"), ("orography", "fall_time", "1000.", "orography"),
                  ("orography", "precip_base", "1.", "orography"), ("orography", "precip_min", "0.01", "orography")],
    "sedfactor": [("", "sedfactor", "\n  - start: 0.\n    uniform: 1.   # erodibility factor, or map: ['efac', 'f']", "sedfactor[]")],
    "spl_transport": [("spl", "G", "1.", "spl")],
    "nlspl": [("spl", "n", "1.5", "spl")],
    "soil": [("soil", "soilK", "1.e-5", "soil"), ("soil", "maxProd", "5.e-5", "soil"),
             ("soil", "depthProd", "0.5", "soil"), ("soil", "uniform", "1.", "soil")],
    "deposition": [("domain", "nodep", "False", "domain")],
    "hillslope": [("diffusion", "hillslopeKa", "0.02", "diffusion"), ("diffusion", "hillslopeKm", "0.2", "diffusion")],
    "channelQ": [("diffusion", "channelQ", "1.e6", "diffusion")],
    "marine_deposition": [("domain", "seadepo", "True", "domain"), ("diffusion", "nonlinKm", "10.", "diffusion"),
                          ("diffusion", "clinSlp", "1.e-5", "diffusion")],
    "gravity": [("gravity", "fraction", "0.2", "gravity"), ("gravity", "runout", "5.e4", "gravity")],
    "waves": [("waves", "climate", "\n    - from: 270.   # bearing the waves come from (degrees)\n      weight: 1.", "waves.climate[]")],
    "coastal_erosion": [("coastal_erosion", "rate", "1.e-3", "coastal_erosion")],
    "strata": [("time", "strat", "{strat}", "time"), ("compaction", "phis", "0.49", "compaction"),
               ("compaction", "z0s", "3700.", "compaction")],
    "dual": [("strata", "dual", "True", "strata")],
    "lithologies": [("strata", "lithologies", "\n    - name: silt   # one more lithology\n      size: fine", "strata.lithologies[]")],
    "provenance": [("provenance", "classes", "2", "provenance"),
                   ("provenance", "source", "['provenance', 'rock']", "provenance")],
    "sst": [("sst", "events", "\n    - start: 0.\n      uniform: 25.   # sea-surface temperature (degC)", "sst.events[]")],
    "carbonate": [("carbonate", "factories", "\n    - name: euphotic\n      growth: 500.      # m/Myr\n      extinction: 0.1   # 1/m\n      saturation: 60.   # W/m2\n      framework: 0.5    # share built in place (reef)", "carbonate.factories[]")],
    "pelagic": [("", "pelagic", "\n  - start: 0.\n    carbonate: 10.    # m/Myr\n    lysocline: 3500.  # m\n    ccd: 4500.        # m", "pelagic[]")],
    "carbonate_supply": [("carbonate_supply", "species", "'carbonate'", "carbonate_supply"),
                         ("carbonate_supply", "yield", "1.", "carbonate_supply")],
    "karst": [("karst", "rate", "5.e-5", "karst")],
    "groundwater": [("groundwater", "Ksat", "3.65e4", "groundwater"), ("groundwater", "specific_yield", "0.1", "groundwater"),
                    ("groundwater", "aquifer_base", "50.", "groundwater"), ("groundwater", "infiltration", "0.3", "groundwater")],
    "duricrust": [("groundwater.duricrust", "form_rate", "1.e-4", "groundwater.duricrust"),
                  ("groundwater.duricrust", "max_thickness", "5.", "groundwater.duricrust"),
                  ("groundwater.duricrust", "fringe_depth", "3.", "groundwater.duricrust"),
                  ("groundwater.duricrust", "fringe_width", "2.", "groundwater.duricrust")],
    "geochem": [("groundwater.geochem", "species", "\n        - name: carbonate\n          weatherability: 1.\n          precip_rate: 1.e-6", "groundwater.geochem.species[]")],
    "river_load": [("groundwater.geochem", "river_load", "True", "groundwater.geochem")],
    "marine_coupling": [("groundwater.geochem", "marine_coupling", "True", "groundwater.geochem")],
    "ice": [("ice", "hela", "2000.", "ice"), ("ice", "hice", "2400.", "ice")],
    "flexure_flat": [("flexure", "method", "'fem'", "flexure"), ("flexure", "thick", "3.e4", "flexure"),
                     ("flexure", "rhoc", "2300.", "flexure"), ("flexure", "rhoa", "3300.", "flexure")],
    "flexure_global": [("flexure", "method", "'global'", "flexure"), ("flexure", "thick", "3.e4", "flexure"),
                       ("flexure", "rhoc", "2300.", "flexure"), ("flexure", "rhoa", "3300.", "flexure")],
    "temap": [("", "temap", "\n  - start: 0.\n    uniform: 3.e4   # elastic thickness (m), or map: ['te', 't']", "temap[]")],
    "tectonics_vertical": [("", "tectonics", "\n  - start: 0.\n    end: {end}\n    upsub: ['uplift', 'u']   # npz: uplift rate per vertex (m/yr)", "tectonics[]")],
    "advection": [("domain", "advect", "'iioe2'", "domain"),
                  ("", "tectonics", "\n  - start: 0.\n    end: {end}\n    hdisp: ['hdisp', 'h']   # npz: horizontal velocity per vertex (m/yr)", "tectonics[]")],
    "plate_interp": [("domain", "advect", "'interp'", "domain"),
                     ("", "tectonics", "\n  - start: 0.\n    end: {end}\n    hdisp: ['hdisp', 'h']   # npz: horizontal velocity per vertex (m/yr)", "tectonics[]")],
    "interp_remap": [("domain", "interp_remap", "1", "domain")],
    "cyclic": [("domain", "bc", "'cfcf'", "domain")],
}


def _yaml_entries(inputs):
    """MOD_YAML with the comment of each key (unit and meaning, from the
    input reference) and the anchor of its block."""
    by_path = {b["path"]: b for b in inputs}
    out = {}
    for mid, rows in MOD_YAML.items():
        ents = []
        for block, key, value, ref in rows:
            com = ""
            if ref and ref in by_path:
                for k in by_path[ref]["keys"]:
                    if k["key"].split()[0] == key:
                        mean = re.split(r"(?<=[a-z0-9)]{3})\. (?=[A-Z])", k.get("meaning", ""))[0].rstrip(".")
                        unit = k.get("unit", "")
                        com = (unit + "; " if unit and unit != "-" else "") + mean
                        break
            if len(com) > 110:
                com = com[:107].rsplit(" ", 1)[0] + "..."
            ents.append({"block": block, "key": key, "value": value, "comment": com,
                         "anchor": _anchor(ref) if ref else ""})
        out[mid] = ents
    return out


# measured extra time of each option, % of the base step (rows of the cost
# table, docs/user_guide/scaling/option_cost_table.rst, "added time")
COST_ROWS = {
    "transport-limited SPL": "spl_transport", "soil": "soil", "flexure": "flexure",
    "orography": "orography", "ice": "ice", "stratigraphy": "strata",
    "dual lithology": "dual", "provenance (3 classes)": "provenance",
    "groundwater": "groundwater", "+ duricrust": "duricrust", "+ geochemistry": "geochem",
    "gravity currents": "gravity", "pelagic rain": "pelagic", "carbonates": "carbonate",
    "+ waves": "waves", "+ coastal erosion": "coastal_erosion", "+ karst": "karst",
}

# where each module works: flat (regional) mesh, sphere (global) or both
MOD_DOMAIN = {"orography": "flat", "flexure_flat": "flat", "cyclic": "flat",
              "flexure_global": "sphere"}


# the same numbers (Nov 2026), for a build without the cost table
COST_FALLBACK = {"spl_transport": 33, "soil": 25, "orography": 2, "ice": 11, "strata": 2,
                 "dual": 5, "provenance": 23, "groundwater": 37, "duricrust": 0, "geochem": 26,
                 "gravity": 21, "pelagic": 40, "carbonate": 56, "waves": 7,
                 "coastal_erosion": 0, "karst": 1, "flexure_flat": 1, "flexure_global": 1}


def _has_cost_page():
    return os.path.exists(os.path.join(os.path.dirname(HERE), "cost.rst"))


def _costs():
    """Added time (% of the base step at 1 rank) per module, from the table."""
    p = os.path.join(os.path.dirname(HERE), "scaling", "option_cost_table.rst")
    out = {}
    if not os.path.exists(p):
        return {k: float(v) for k, v in COST_FALLBACK.items()}
    rows = re.findall(r"   \* - (.*)\n     - (.*)\n     - (.*)\n     - (.*)\n     - (.*)\n", open(p).read())
    for label, _, one, _, added in rows:
        label = label.strip()
        if label in COST_ROWS:
            try:
                out[COST_ROWS[label]] = float(added.strip() or 0)
            except ValueError:
                pass
    if "flexure" in out:
        out["flexure_flat"] = out["flexure_global"] = out.pop("flexure")
    return out


def _params(inputs):
    """Per module: the input keys it sets, with unit, default and range."""
    by_path = {b["path"]: b for b in inputs}
    out = {}
    for mid, spec in MOD_INPUTS.items():
        rows = []
        for path, keys in spec:
            b = by_path.get(path)
            if not b:
                continue
            for k in b["keys"]:
                name = k["key"].split()[0]
                if keys is None or name in keys:
                    rows.append({"key": (path + ": " if path not in ("domain", "time", "spl", "diffusion")
                                         or keys is None else path + ": ") + name,
                                 "unit": k.get("unit", ""), "default": str(k.get("default", "")),
                                 "range": k.get("range", "")})
        if rows:
            out[mid] = rows
    return out


def _examples():
    """Catalogue of the goSPL-examples inputs, limited to the examples that
    have a page in this build's gallery."""
    p = os.path.join(HERE, "examples.json")
    if not os.path.exists(p):
        return []
    gen = os.path.join(os.path.dirname(os.path.dirname(HERE)), "examples", "generated")
    recs = json.load(open(p))
    if os.path.isdir(gen):
        recs = [r for r in recs if os.path.isdir(os.path.join(gen, r["example"]))]
    return recs


def _modules(data):
    data = _strip_refs(data)
    costs = _costs()
    ip = os.path.join(HERE, "inputs.json")
    inputs = json.load(open(ip)) if os.path.exists(ip) else []
    params = _params(inputs)
    yaml_ents = _yaml_entries(inputs)
    ids = {m["id"] for m in data["modules"]}
    for m in data["modules"]:
        m["domain"] = MOD_DOMAIN.get(m["id"], "both")
        m["cost"] = costs.get(m["id"], 0.0)
        if m["id"] in params:
            m["params"] = params[m["id"]]
        m["tech"] = MOD_TECH.get(m["id"], "")
        if m["id"] in yaml_ents:
            m["yaml"] = yaml_ents[m["id"]]
    ex = _examples()
    for r in ex:
        r["modules"] = [i for i in r["modules"] if i in ids]
    data["examples"] = ex
    data["meshes"] = MESH_LINKS
    js = json.dumps(data, ensure_ascii=False).replace("</", "<\\/")
    return "".join([
        ".. _module_chooser:\n\n", "=" * 26 + "\n", "Which modules do I need?\n", "=" * 26 + "\n\n",
        "Choose a regional (flat mesh) or global (sphere) model and tick the\n"
        "processes you want to simulate. The graph shows what each one needs (red\n"
        "until it is selected, green when complete), the page lists the examples\n"
        "that come closest to your set-up with their mesh, resolution and time span,\n"
        "the main parameters of each module with their usual values, and an estimate\n"
        "of the time and memory your run will need. The units and ranges of every\n"
        "key are on :ref:`input_reference`%s.\n\n" % (
            "; the measured cost of each option on :ref:`option_cost`" if _has_cost_page() else ""),
        ".. raw:: html\n\n",
        "   <div class=\"gospl-modchart\"></div>\n",
        "   <script>window.GOSPL_MODULES = %s; window.GOSPL_DOCROOT = \"../../\";</script>\n" % js,
        "   <script src=\"../../_static/modchart.js\"></script>\n\n",
        ".. note::\n\n",
        "   The graph is built from the checks goSPL makes when it reads the input\n",
        "   file: a module is *required* when goSPL stops without it, ignores the\n",
        "   option, or would compute nothing useful. The resource estimate scales the\n",
        "   costs measured on a 47,672-node margin%s with the number\n" % (
            " (:ref:`option_cost`)" if _has_cost_page() else ""),
        "   of nodes and steps: read it as an order of magnitude. When in doubt, the\n",
        "   message goSPL prints at start-up is the reference.\n",
    ])


def build():
    io = [("inputs.json", "inputs.rst", _inputs), ("outputs.json", "outputs_ref.rst", _outputs),
          ("modules.json", "modules.rst", _modules)]
    for src, dst, fn in io:
        p = os.path.join(HERE, src)
        if not os.path.exists(p):
            continue
        open(os.path.join(HERE, dst), "w").write(fn(json.load(open(p))))


if __name__ == "__main__":
    build()
