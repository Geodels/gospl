#!/usr/bin/env python3
"""
Derive goSPL glacier-geometry maps (``hela`` / ``hice``) from a paleo-climate
temperature map, by inverting the atmospheric lapse rate.

The equilibrium-line altitude (ELA) is the elevation at which the air
temperature reaches a calibrated threshold ``T_ELA``. With a lapse rate
:math:`\\Gamma` (deg C per metre),

    ELA(x) = (T_ref(x) - T_ELA) / Gamma

where ``T_ref`` is the temperature reduced to sea level. If your temperature map
is given at the local surface elevation ``z(x)`` (the usual near-surface output
of a climate model, after remapping to the goSPL mesh), reduce it on the fly:

    T_ref(x) = T(x) + Gamma * z(x)   ==>   ELA(x) = z(x) + (T(x) - T_ELA) / Gamma

The ice-cap altitude (top of the accumulation ramp) is set a fixed band above
the ELA: ``hice = hela + band``.

The output is a single ``.npz`` with per-vertex ``hela`` and ``hice`` arrays, in
the convention goSPL expects for the ``ice`` maps. Run it once per paleo-climate
snapshot, then list the outputs in the ``ice.glaciers`` time series (each entry
``hela: [<out-without-.npz>, 'hela']`` etc.). See docs/user_guide/surfproc.rst
and docs/tech_guide/ice.rst.

NOTE
----
- The temperature and elevation maps must already be per-vertex on the **goSPL
  mesh** (same shape/order as the mesh ``npdata`` arrays). Remapping a climate
  model's native grid onto the mesh is a separate, upstream step.
- This derives the ELA *position* only. goSPL's ablation magnitude stays
  precipitation-scaled (the ELA-ramp SMB); it is not a degree-day / energy-balance
  melt model.

Running it — terminal
---------------------
One snapshot (prints a ready-to-paste ``ice.glaciers`` entry via ``--start``)::

    python scripts/ela_from_temperature.py \\
        --temperature climate/t2m_21ka.npz --t-key t2m \\
        --reference surface --elevation input/mesh.npz --z-key z \\
        --lapse 0.0065 --t-ela -2.0 --band 400 \\
        --out input/ela_21ka.npz --start -21000

Whole paleo-climate history (one map per snapshot) in a shell loop::

    for t in 21000 18000 15000 12000; do
        python scripts/ela_from_temperature.py \\
            --temperature climate/t2m_${t}ka.npz --t-key t2m \\
            --reference surface --elevation input/mesh.npz --z-key z \\
            --lapse 0.0065 --t-ela -2.0 --band 400 \\
            --out input/ela_${t}.npz --start -${t}
    done

``python scripts/ela_from_temperature.py --help`` lists every option.

Running it — Jupyter notebook
-----------------------------
Either shell out to the script with the ``!`` magic::

    !python scripts/ela_from_temperature.py \\
        --temperature climate/t2m_21ka.npz --t-key t2m \\
        --reference surface --elevation input/mesh.npz --z-key z \\
        --lapse 0.0065 --t-ela -2.0 --band 400 --out input/ela_21ka.npz

or import the conversion and loop over snapshots, building the ``glaciers``
list programmatically::

    import sys, numpy as np
    sys.path.append("scripts")              # make the helper importable
    from ela_from_temperature import derive_ela

    z = np.load("input/mesh.npz")["z"]      # mesh elevation (per vertex)
    snapshots = {-21000: "climate/t2m_21ka.npz",
                 -18000: "climate/t2m_18ka.npz"}

    glaciers = []
    for t, fpath in sorted(snapshots.items()):
        T = np.load(fpath)["t2m"]
        hela, hice = derive_ela(T, lapse=0.0065, t_ela=-2.0, band=400.0,
                                reference="surface", elevation=z)
        out = f"input/ela_{int(-t/1000)}ka.npz"
        np.savez(out, hela=hela, hice=hice)
        glaciers.append({"start": float(t),
                         "hela": [out[:-4], "hela"],
                         "hice": [out[:-4], "hice"]})

    glaciers   # paste under the YAML `ice:` block as `glaciers:`
"""

import argparse
import sys

import numpy as np


def _load(path, key, what):
    try:
        data = np.load(path)
    except IOError:
        sys.exit("error: cannot open %s map file: %s" % (what, path))
    if key not in data:
        sys.exit(
            "error: field '%s' missing from %s file %s (available: %s)"
            % (key, what, path, ", ".join(data.files))
        )
    return np.asarray(data[key], dtype=np.float64)


def derive_ela(temperature, lapse, t_ela, band, reference="surface", elevation=None):
    """
    Return ``(hela, hice)`` per-vertex arrays from a temperature map.

    :arg temperature: per-vertex temperature (deg C).
    :arg lapse: atmospheric lapse rate (deg C per metre, e.g. 0.0065).
    :arg t_ela: temperature at the equilibrium line (deg C).
    :arg band: hice - hela accumulation-band width (m).
    :arg reference: 'surface' (T at the local elevation; needs ``elevation``)
        or 'sealevel' (T already reduced to sea level).
    :arg elevation: per-vertex surface elevation (m); required for 'surface'.
    """
    if lapse <= 0.0:
        raise ValueError("lapse rate must be positive (deg C per metre).")
    if reference == "surface":
        if elevation is None:
            raise ValueError("reference='surface' requires an elevation map.")
        hela = elevation + (temperature - t_ela) / lapse
    elif reference == "sealevel":
        hela = (temperature - t_ela) / lapse
    else:
        raise ValueError("reference must be 'surface' or 'sealevel'.")
    hice = hela + band
    return hela, hice


def main(argv=None):
    p = argparse.ArgumentParser(
        description="Derive goSPL hela/hice maps from a temperature map.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--temperature", required=True, help="temperature .npz file")
    p.add_argument("--t-key", default="t2m", help="temperature field name in the npz")
    p.add_argument(
        "--reference",
        choices=["surface", "sealevel"],
        default="surface",
        help="elevation the temperature is referenced to",
    )
    p.add_argument("--elevation", help="elevation .npz file (for reference=surface)")
    p.add_argument("--z-key", default="z", help="elevation field name in the npz")
    p.add_argument(
        "--lapse", type=float, default=0.0065, help="lapse rate (deg C per metre)"
    )
    p.add_argument(
        "--t-ela", type=float, required=True, help="temperature at the ELA (deg C)"
    )
    p.add_argument(
        "--band", type=float, default=400.0, help="hice - hela band width (m)"
    )
    p.add_argument("--out", required=True, help="output .npz path")
    p.add_argument("--ela-key", default="hela", help="ELA field name in the output")
    p.add_argument("--hice-key", default="hice", help="ice-cap field name in the output")
    p.add_argument(
        "--start",
        type=float,
        default=None,
        help="if set, print a ready-to-paste ice.glaciers YAML entry for this time",
    )
    args = p.parse_args(argv)

    temperature = _load(args.temperature, args.t_key, "temperature")
    elevation = None
    if args.reference == "surface":
        if not args.elevation:
            sys.exit("error: --reference surface requires --elevation")
        elevation = _load(args.elevation, args.z_key, "elevation")
        if elevation.shape != temperature.shape:
            sys.exit(
                "error: elevation %s and temperature %s shapes differ"
                % (elevation.shape, temperature.shape)
            )

    if not np.isfinite(temperature).all():
        print(
            "warning: temperature map has non-finite values; the derived ELA will "
            "be non-finite there. Fill/mask them before running goSPL.",
            file=sys.stderr,
        )

    hela, hice = derive_ela(
        temperature,
        lapse=args.lapse,
        t_ela=args.t_ela,
        band=args.band,
        reference=args.reference,
        elevation=elevation,
    )

    np.savez(args.out, **{args.ela_key: hela, args.hice_key: hice})
    finite = np.isfinite(hela)
    print(
        "wrote %s (%d vertices): hela %.0f..%.0f m, hice = hela + %.0f m"
        % (
            args.out,
            hela.size,
            float(np.min(hela[finite])) if finite.any() else float("nan"),
            float(np.max(hela[finite])) if finite.any() else float("nan"),
            args.band,
        )
    )

    if args.start is not None:
        base = args.out[:-4] if args.out.endswith(".npz") else args.out
        print("\n# ice.glaciers entry:")
        print("  - start: %g" % args.start)
        print("    hela: ['%s', '%s']" % (base, args.ela_key))
        print("    hice: ['%s', '%s']" % (base, args.hice_key))

    return 0


if __name__ == "__main__":
    sys.exit(main())
