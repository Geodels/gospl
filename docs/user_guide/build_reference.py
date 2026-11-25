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


def _costs():
    """Added time (% of the base step at 1 rank) per module, from the table."""
    p = os.path.join(os.path.dirname(HERE), "scaling", "option_cost_table.rst")
    out = {}
    if not os.path.exists(p):
        return out
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
    params = _params(json.load(open(ip))) if os.path.exists(ip) else {}
    ids = {m["id"] for m in data["modules"]}
    for m in data["modules"]:
        m["domain"] = MOD_DOMAIN.get(m["id"], "both")
        m["cost"] = costs.get(m["id"], 0.0)
        if m["id"] in params:
            m["params"] = params[m["id"]]
    ex = _examples()
    for r in ex:
        r["modules"] = [i for i in r["modules"] if i in ids]
    data["examples"] = ex
    js = json.dumps(data, ensure_ascii=False).replace("</", "<\\/")
    return "".join([
        ".. _module_chooser:\n\n", "=" * 26 + "\n", "Which modules do I need?\n", "=" * 26 + "\n\n",
        "Choose a regional (flat mesh) or global (sphere) model and tick the\n"
        "processes you want to simulate. The graph shows what each one needs (red\n"
        "until it is selected, green when complete), the page lists the examples\n"
        "that come closest to your set-up with their mesh, resolution and time span,\n"
        "the main parameters of each module with their usual values, and an estimate\n"
        "of the time and memory your run will need. The units and ranges of every\n"
        "key are on :ref:`input_reference`; the measured cost of each option on\n"
        ":ref:`option_cost`.\n\n",
        ".. raw:: html\n\n",
        "   <div class=\"gospl-modchart\"></div>\n",
        "   <script>window.GOSPL_MODULES = %s; window.GOSPL_DOCROOT = \"../../\";</script>\n" % js,
        "   <script src=\"../../_static/modchart.js\"></script>\n\n",
        ".. note::\n\n",
        "   The graph is built from the checks goSPL makes when it reads the input\n",
        "   file: a module is *required* when goSPL stops without it, ignores the\n",
        "   option, or would compute nothing useful. The resource estimate scales the\n",
        "   costs measured on a 47,672-node margin (:ref:`option_cost`) with the number\n",
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
