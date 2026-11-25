"""
Write ``examples.json``: what each goSPL-examples input sets up, for the
module chooser (``docs/user_guide/reference/modules.rst``).

For every input file of every example (``<group>/<example>/*.yml``) it
records the domain (regional flat mesh or global sphere), the mesh (nodes,
mean spacing, extent), the time span, step and stratigraphic interval, and
the modules the input switches on (ids of ``modules.json``). It reads the
meshes, which are not in the goSPL-examples repository for every example
(several are built by scripts), so run it on a checkout where the examples
have been set up, and commit the result:

    python docs/user_guide/reference/make_examples_catalog.py ~/Workspace/goSPL-examples
"""

import glob
import json
import os
import sys

import numpy as np
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
GROUPS = [("Local-examples", "regional"), ("Global-examples", "global")]


def _mesh(path):
    """Nodes, mean edge length (m), flat flag, extent (km) or radius (km)."""
    try:
        d = np.load(path)
    except Exception:
        return None
    v, c = np.asarray(d["v"], float), np.asarray(d["c"])
    if v.shape[1] == 2:
        v = np.column_stack([v, np.zeros(len(v))])
    flat = np.ptp(v[:, 2]) < 1.0e-6 * max(np.ptp(v[:, 0]), 1.0)
    cs = c[:: max(1, len(c) // 20000)]
    e = np.concatenate([np.linalg.norm(v[cs[:, i]] - v[cs[:, (i + 1) % 3]], axis=1) for i in range(3)])
    out = {"nodes": int(len(v)), "spacing_m": float(np.mean(e)), "flat": bool(flat)}
    if flat:
        out["extent_km"] = [round(float(np.ptp(v[:, 0])) / 1e3, 1), round(float(np.ptp(v[:, 1])) / 1e3, 1)]
    else:
        out["radius_km"] = round(float(np.linalg.norm(v, axis=1).mean()) / 1e3, 0)
    return out


def _modules(y, flat):
    """Module ids (modules.json) switched on by one input."""
    m = {"core"}
    dom, time, spl = y.get("domain") or {}, y.get("time") or {}, y.get("spl") or {}
    if "sea" in y:
        m.add("sea")
    clim = y.get("climate") or []
    if clim:
        m.add("climate")
        if any(("evap_uniform" in r or "evap_map" in r) for r in clim if isinstance(r, dict)):
            m.add("evaporation")
    for key, mid in (("orography", "orography"), ("sedfactor", "sedfactor"), ("soil", "soil"),
                     ("gravity", "gravity"), ("waves", "waves"), ("coastal_erosion", "coastal_erosion"),
                     ("provenance", "provenance"), ("sst", "sst"), ("carbonate", "carbonate"),
                     ("pelagic", "pelagic"), ("carbonate_supply", "carbonate_supply"),
                     ("karst", "karst"), ("ice", "ice"), ("temap", "temap")):
        if key in y and y[key] is not None:
            m.add(mid)
    if float(spl.get("G", 0) or 0) > 0:
        m.add("spl_transport")
    if float(spl.get("n", 1) or 1) != 1:
        m.add("nlspl")
    nodep = bool(dom.get("nodep", False))
    if not nodep:
        m.add("deposition")
        if dom.get("seadepo", True):
            m.add("marine_deposition")
    diff = y.get("diffusion") or {}
    if float(diff.get("hillslopeKa", 0) or 0) > 0 or float(diff.get("hillslopeKm", 0) or 0) > 0:
        m.add("hillslope")
    if diff.get("channelQ"):
        m.add("channelQ")
    if float(time.get("strat", 0) or 0) > 0:
        m.add("strata")
    st = y.get("strata") or {}
    if st.get("dual"):
        m.add("dual")
    if st.get("lithologies"):
        m.add("lithologies")
    gw = y.get("groundwater")
    if gw is not None:
        m.add("groundwater")
        if gw.get("duricrust") is not None:
            m.add("duricrust")
        gc = gw.get("geochem")
        if gc is not None:
            m.add("geochem")
            if gc.get("river_load"):
                m.add("river_load")
            if gc.get("marine_coupling"):
                m.add("marine_coupling")
    if "flexure" in y:
        meth = (y["flexure"] or {}).get("method")
        m.add("flexure_flat" if (meth == "fem" or (meth is None and flat)) else "flexure_global")
    tec = y.get("tectonics") or []
    if any(isinstance(r, dict) and "upsub" in r for r in tec):
        m.add("tectonics_vertical")
    if any(isinstance(r, dict) and "hdisp" in r for r in tec):
        m.add("plate_interp" if dom.get("advect") == "interp" else "advection")
    if dom.get("interp_remap"):
        m.add("interp_remap")
    if "c" in str(dom.get("bc", "")):
        m.add("cyclic")
    return sorted(m)


def main(src):
    out = []
    for group, kind in GROUPS:
        for edir in sorted(glob.glob(os.path.join(src, group, "*"))):
            if not os.path.isdir(edir):
                continue
            name = os.path.basename(edir)
            for yml in sorted(glob.glob(os.path.join(edir, "*.yml"))):
                try:
                    y = yaml.safe_load(open(yml))
                except Exception:
                    continue
                if not isinstance(y, dict) or "domain" not in y or "time" not in y:
                    continue
                npd = (y["domain"] or {}).get("npdata")
                mesh = _mesh(os.path.join(edir, npd[0] + ".npz")) if npd else None
                t = y["time"]
                start, end, dt = float(t.get("start", 0)), float(t.get("end", 0)), float(t.get("dt", 0))
                rec = {"example": name, "group": kind, "input": os.path.basename(yml),
                       "page": "examples/generated/%s/index" % name,
                       "duration_yr": end - start, "dt_yr": dt,
                       "steps": int(round((end - start) / dt)) if dt else None,
                       "strat_yr": float(t.get("strat", 0) or 0) or None,
                       "modules": _modules(y, mesh["flat"] if mesh else kind == "regional")}
                if mesh:
                    rec.update(mesh)
                out.append(rec)
    json.dump(out, open(os.path.join(HERE, "examples.json"), "w"), indent=1)
    print("%d inputs from %d examples" % (len(out), len({r["example"] for r in out})))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/Workspace/goSPL-examples"))
