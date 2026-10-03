"""
Mesh-native orographic precipitation.

Protects: AGENTS.md > Fixed (wind_dir sign).

Split out of the former tests/test_regression.py. Run this group alone with
`pytest -m orography`; see tests/README.md for the marker list.
"""

from __future__ import annotations

import numpy as np
import pytest

# Skip the whole module (rather than erroring at collection) when the
# goSPL runtime stack is not installed.
pytest.importorskip(
    "gospl.tools.inputparser",
    reason="goSPL runtime deps (petsc4py / ruamel.yaml / scipy) not installed",
)

pytestmark = [pytest.mark.orography]


def test_orography_rain_shadow(oro_hill_model):
    """
    Protects: the mesh-native orographic precipitation solver (parallel, no
    regular grid / FFT). Smith-Barstad is recast as two steady advection-
    relaxation equations on the DMPlex (cloud water → hydrometeors → precip)
    with the terrain-forced uplift source Cw·(v·∇h). The fixture is an ocean →
    coastal mountain → ocean profile with wind from the west. Checks that:
      - the orographic rain is computed, finite and positive everywhere;
      - the uplift is driven by topography ABOVE sea level only: the windward
        ocean (submarine bathymetry) gets the background rate, no orographic
        enhancement (the elevation is clamped to sea level in the source);
      - the windward/lee RAIN SHADOW is reproduced: the coastal mountain is far
        wetter than the dried lee, and the upwind ocean is wetter than the
        (shadowed) downwind ocean;
      - the precipitation peak sits at/just downwind of the crest (advective
        shift), not on the far windward side.
    """
    model = oro_hill_model
    assert model.oroOn and model.flatModel

    model.cptOrography()

    own = model.inIDs == 1                 # owned nodes (exclude halo)
    x = model.lcoords[:, 0]
    rain = model.rainVal
    assert np.isfinite(rain[own]).all(), "orographic rain is not finite"
    assert (rain[own] > 0.0).all(), "orographic rain must stay positive"

    def region_mean(lo, hi):
        m = (x >= lo) & (x < hi) & own
        return rain[m].mean()

    # Background (non-orographic) rate the solver adds everywhere.
    background = model.oro_precip_base * 0.366 * model.rainfall_frequency

    windward_ocean = region_mean(1000.0, 4000.0)    # submarine, upwind
    mountain = region_mean(6500.0, 9000.0)          # coastal windward flank (land)
    lee_ocean = region_mean(12500.0, 15000.0)       # shadowed, downwind ocean

    # Sea-level clamp: no orographic forcing over the windward submarine
    # bathymetry — that region carries only the background rate.
    assert abs(windward_ocean - background) < 0.2 * background, (
        f"windward ocean got orographic rain ({windward_ocean:.3f}) instead of "
        f"background ({background:.3f}) — sea-level clamp not applied"
    )
    # ... and it is far drier than the coastal mountain it sits next to.
    assert mountain > 10.0 * windward_ocean, "no orographic uplift over the coast"

    # Rain shadow: the dried downwind ocean is much drier than the windward
    # mountain, and drier than the upwind ocean.
    assert mountain > 5.0 * lee_ocean, "no rain shadow over the lee"
    assert lee_ocean < windward_ocean, "lee ocean not shadowed vs upwind ocean"

    # Precip maximum is at/just downwind of the crest (advective shift).
    xpeak = x[own][np.argmax(rain[own])]
    assert 9000.0 <= xpeak <= 11500.0, f"precip peak misplaced at x={xpeak:.0f}"


def test_orography_wind_direction_convention(oro_hill_model):
    """
    Protects: `wind_dir` is the meteorological "comes FROM" bearing
    (0=N, 90=E, 180=S, 270=W) consistently in BOTH the East and North
    components (`addprocess._windVector`).

    Silent failure prevented: the North component was historically `+cos`
    (a "goes toward" sign) while the East component was `-sin` (a "comes from"
    sign). The mix was invisible on E/W winds (cos of 90/270 is 0) and on the
    rain-shadow fixture (wind_dir=270), so a northerly (wind_dir=0) blew toward
    the North instead of the South. A wind COMING FROM bearing d blows TOWARD
    d+180, i.e. velocity (u_east, v_north) = speed * (-sin d, -cos d).
    """
    model = oro_hill_model
    speed = model.wind_speed
    # (wind_dir, expected unit blow direction): comes-from -> blows opposite.
    cases = {
        0.0:   (0.0, -1.0),   # from N -> toward S (-y)
        90.0:  (-1.0, 0.0),   # from E -> toward W (-x)
        180.0: (0.0, 1.0),    # from S -> toward N (+y)
        270.0: (1.0, 0.0),    # from W -> toward E (+x)  (rain-shadow fixture)
    }
    for d, (ex, ey) in cases.items():
        model.wind_dir = d
        u, v = model._windVector()
        assert np.allclose([u, v], [ex * speed, ey * speed], atol=1.0e-9 * speed + 1.0e-12), (
            f"wind_dir={d} (comes-from) should blow toward "
            f"({ex:+.0f},{ey:+.0f})*speed; got ({u/speed:+.3f},{v/speed:+.3f})*speed"
        )
