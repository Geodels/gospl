"""
Tests for the scripts/ela_from_temperature.py preprocessing helper, which derives
goSPL hela/hice maps from a paleo-climate temperature map by lapse-rate inversion.
"""

import importlib.util
import os

import numpy as np
import pytest

_HELPER = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "scripts",
    "ela_from_temperature.py",
)


def _load_helper():
    if not os.path.exists(_HELPER):
        pytest.skip("ela_from_temperature.py helper not present")
    spec = importlib.util.spec_from_file_location("ela_from_temperature", _HELPER)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_derive_ela_surface_reference():
    """T given at the local surface: ELA = z + (T - T_ELA)/lapse, hice = ELA + band."""
    m = _load_helper()
    z = np.array([0.0, 1000.0, 2000.0, 3000.0])
    lapse, t_ela, band = 0.0065, -2.0, 400.0
    # A perfect lapse-rate field (15 C at sea level) -> a single ELA everywhere.
    t2m = 15.0 - lapse * z
    hela, hice = m.derive_ela(t2m, lapse, t_ela, band, reference="surface", elevation=z)

    expected = (15.0 - t_ela) / lapse                 # = 2615.38 m
    assert np.allclose(hela, expected)
    assert np.allclose(hice, hela + band)
    # By construction the air temperature at the ELA equals T_ELA.
    assert np.allclose(15.0 - lapse * hela, t_ela)


def test_derive_ela_sealevel_reference():
    """T already reduced to sea level: ELA = (T - T_ELA)/lapse, no elevation needed."""
    m = _load_helper()
    lapse, t_ela = 0.0065, -2.0
    t_sl = np.array([15.0, 10.0, 0.0])               # warmer -> higher ELA
    hela, hice = m.derive_ela(t_sl, lapse, t_ela, 300.0, reference="sealevel")
    assert np.allclose(hela, (t_sl - t_ela) / lapse)
    assert (np.diff(hela) < 0).all()                  # cooler sea-level T -> lower ELA
    assert np.allclose(hice - hela, 300.0)


def test_derive_ela_guards():
    """Invalid configurations raise rather than silently producing garbage."""
    m = _load_helper()
    with pytest.raises(ValueError):
        m.derive_ela(np.zeros(3), 0.0, -2.0, 400.0, reference="sealevel")  # lapse <= 0
    with pytest.raises(ValueError):
        m.derive_ela(np.zeros(3), 0.0065, -2.0, 400.0, reference="surface")  # no elev
    with pytest.raises(ValueError):
        m.derive_ela(np.zeros(3), 0.0065, -2.0, 400.0, reference="bogus")
