"""
Analytical benchmark — steady-state Dupuit-Boussinesq water table.

Validates the implicit groundwater head solve (`gospl.flow.gwplex._solveHead`,
DESIGN_WATERTABLE_DURICRUST.md Phase 2) against the exact 1-D Dupuit-Boussinesq
solution for an unconfined aquifer on a flat impermeable base, uniform recharge
`R`, conductivity `K`, draining to a fixed-head edge with a no-flow divide at the
far edge:

    s(x)^2 = s_d^2 + (R/K) (2 W x - x^2),    x in [0, W]

with `s = h - z_bed` the saturated thickness, `s_d` at the drain (x=0), and the
divide at x=W (ds/dx = 0). Derivation: steady Dupuit `d/dx(K s ds/dx) = -R`,
integrate with the no-flow condition at x=W.

Setup (a strip that ramps UP from a west drain, over a flat base at z_bed=0):
  - the surface ramps `z(x) = z_drain + slope*x` UP to the east — a monotone
    west-draining slope, so the flat-model drainage routes every cell to the open
    west edge and produces NO interior pits (a truly flat plateau spuriously ponds
    a few cells, which would pin the water table to the surface and corrupt the
    test);
  - the WEST edge is the fixed drain ('f' → head pinned to `z_drain` via seepage);
    N/E/S edges are walls ('w', no-flow), enforcing the Dupuit no-flow divide;
  - a per-vertex `aquifer_base = z` map holds the base flat at datum (`z_bed = 0`),
    so `h` IS the saturated thickness and the flat-base analytic parabola applies
    verbatim; the ramp is steep enough that the parabola (max ~80 m) stays below
    the surface everywhere (no interior seepage clip);
  - uniform recharge `R = infiltration * rain`.

The surface is held FIXED (erosion is not exercised — we call `updateGroundwater`
repeatedly, not `runProcesses`), so the head relaxes onto the steady Dupuit
profile and is compared column-by-column to the analytic parabola.
"""
import os
import numpy as np
import pytest

scipy = pytest.importorskip("scipy")
from scipy.spatial import Delaunay  # noqa: E402

from petsc4py import PETSc  # noqa: E402
from gospl.model import Model  # noqa: E402

# ---- Domain / aquifer parameters (SI) --------------------------------------
NX, NY = 121, 9          # strip: fine in x (flow direction), a few rows in y
DX = 100.0               # node spacing (m)  -> W = (NX-1)*DX = 12 km
W = (NX - 1) * DX
Z_DRAIN = 50.0           # west-edge drain surface (m) = pinned head = s_d
SLOPE = 0.02             # surface ramp up to the east (>0.0065 keeps wt below z)
KSAT = 3650.0            # hydraulic conductivity (m/yr) ~ 10 m/day
RAIN = 1.0               # uniform precip (m/yr)
INFIL = 0.1              # infiltration fraction -> R = 0.1 m/yr
R = INFIL * RAIN
N_RELAX = 250            # steady-state relaxation iterations
TOL_RMSE_PCT = 5.0       # RMSE of h vs analytic, as % of the head range


def _surface(x):
    """West-draining ramp: z = z_drain + slope*x (monotone, no interior pits)."""
    return Z_DRAIN + SLOPE * x


def _write_mesh(path):
    """Structured strip mesh: v=(x,y,0), surface ramps up from a west drain; a
    per-vertex `abase` = z map holds the aquifer base flat at datum (z_bed=0)."""
    xs = np.arange(NX) * DX
    ys = np.arange(NY) * DX
    X, Y = np.meshgrid(xs, ys)
    x = X.ravel()
    y = Y.ravel()
    v = np.column_stack([x, y, np.zeros_like(x)])
    cells = Delaunay(np.column_stack([x, y])).simplices.astype(np.int64)
    z = _surface(x).astype(np.float64)
    abase = z.copy()                          # aquifer_base = z -> z_bed = 0
    np.savez(path, v=v, c=cells, z=z, abase=abase)
    return x


def _write_yaml(path):
    with open(path, "w") as f:
        f.write(
            "name: dupuit groundwater benchmark\n"
            "domain:\n"
            "    npdata: ['dupuit_mesh','v','c','z']\n"
            "    flowdir: 2\n"
            "    bc: 'wwwf'\n"          # N,E,S wall; W = fixed drain
            "    seadepo: False\n"
            "    nodep: True\n"
            "time:\n"
            "    start: 0.\n    end: 10.\n    tout: 10.\n    dt: 10.\n"
            "spl:\n"
            "    K: 1.0e-20\n    d: 0.\n    m: 0.5\n"   # erosion effectively off
            "diffusion:\n    hillslopeKa: 0.\n    hillslopeKm: 0.\n"
            "groundwater:\n"
            f"    Ksat: {KSAT}\n"
            "    specific_yield: 0.1\n"
            "    aquifer_base: ['dupuit_mesh','abase']\n"   # z_bed = z - abase = 0
            f"    infiltration: {INFIL}\n"
            "    picard_its: 3\n"
            "    seepage_passes: 4\n"
            "sea:\n    position: -1000.\n"          # no marine anywhere
            "climate:\n  - start: 0.\n"
            f"    uniform: {RAIN}\n"
            "output:\n    dir: 'dupuit_out'\n"
        )


@pytest.mark.benchmark
@pytest.mark.slow
def test_dupuit_watertable(tmp_path, monkeypatch):
    """Steady goSPL water table vs the exact Dupuit parabola (RMSE < 5%)."""
    monkeypatch.chdir(tmp_path)
    _write_mesh(tmp_path / "dupuit_mesh.npz")
    _write_yaml(tmp_path / "dupuit.yml")

    m = Model("dupuit.yml", verbose=False, showlog=False)
    try:
        assert m.gwOn and m.flatModel, "expected a flat groundwater model"
        # One step to populate the flow fields (seaID / drainage / rainVal);
        # erosion is ~off (flat plateau, K~0) so the surface stays flat.
        m.tEnd = m.tNow + 0.5 * m.dt
        m.runProcesses()
        # Relax the head onto steady state on the fixed surface.
        for _ in range(N_RELAX):
            m.updateGroundwater()

        x = m.lcoords[:, 0]
        z = m.hLocal.getArray()
        h = m.headL.getArray()
        own = m.inIDs == 1
        # Interior plateau columns (exclude the pinned drain edge x=0).
        interior = own & (x > 0.5 * DX)

        # Analytic Dupuit parabola: s(x)^2 = s_d^2 + (R/K)(2 W x - x^2), s = h
        # (base at 0). s_d = z_drain (head pinned at the drain).
        s_analytic = np.sqrt(
            Z_DRAIN ** 2 + (R / KSAT) * (2.0 * W * x - x ** 2)
        )
        # Below-surface everywhere (the solve must not have clipped to seepage).
        assert (h[interior] <= z[interior] + 1.0e-6).all(), (
            "water table reached the surface"
        )
        assert (s_analytic[interior] <= z[interior]).all(), (
            "analytic parabola tops the ramp surface — retune SLOPE"
        )

        err = h[interior] - s_analytic[interior]
        rmse = float(np.sqrt(np.mean(err ** 2)))
        head_range = float(s_analytic[interior].max() - s_analytic[interior].min())
        rmse_pct = 100.0 * rmse / max(head_range, 1.0e-9)
        # R^2 of the fit.
        ss_res = float(np.sum(err ** 2))
        ss_tot = float(np.sum((h[interior] - h[interior].mean()) ** 2))
        r2 = 1.0 - ss_res / max(ss_tot, 1.0e-30)

        if PETSc.COMM_WORLD.getRank() == 0:
            print(
                f"\n[dupuit] W={W:.0f} m  R/K={R/KSAT:.2e}  head range={head_range:.2f} m"
                f"  RMSE={rmse:.3f} m ({rmse_pct:.2f}%)  R^2={r2:.4f}"
            )
        assert rmse_pct < TOL_RMSE_PCT, (
            f"steady water table off the Dupuit parabola: RMSE {rmse_pct:.2f}% "
            f"(> {TOL_RMSE_PCT}%)"
        )
        assert r2 > 0.98, f"poor fit to Dupuit parabola: R^2={r2:.4f}"
    finally:
        m.destroy()
