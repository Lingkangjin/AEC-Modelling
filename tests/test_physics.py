import numpy as np
import pytest

import aec_modelling as am
from aec_modelling import electrolyte, thermo


def test_thermoneutral_voltage_at_25C():
    assert thermo.thermoneutral_voltage(298.15) == pytest.approx(1.481, abs=2e-3)


def test_thermodynamic_matches_empirical_reversible_voltage():
    T = np.array([298.15, 323.15, 353.15])
    e_thermo = thermo.standard_reversible_voltage(T)
    e_emp = am.electrochemistry.empirical_reversible_voltage(T)
    np.testing.assert_allclose(e_thermo, e_emp, atol=5e-3)


def test_reversible_voltage_increases_with_pressure():
    v1 = am.reversible_voltage(353.15, 1.5, 0.3)
    v30 = am.reversible_voltage(353.15, 30, 0.3)
    assert v30 > v1


def test_electrolyte_properties_reasonable():
    # 30 wt% KOH at 80 degC: ~1.26 kg/L, ~6.7 mol/L, ~1.4 S/cm
    assert electrolyte.density(353.15, 0.3) == pytest.approx(1258, rel=0.02)
    assert electrolyte.molarity(353.15, 0.3) == pytest.approx(6.73, rel=0.02)
    assert 1.2 < electrolyte.conductivity(353.15, 0.3) < 1.6
    assert electrolyte.conductivity(353.15, 0.3) > electrolyte.conductivity(303.15, 0.3)


def test_voltage_breakdown_sums_and_is_monotonic():
    m = am.AECModel(P=16, w=0.25)
    j = np.linspace(10, 500, 50)
    df = m.overpotentials(j, 343.15)
    np.testing.assert_allclose(df.V_rev + df.V_act + df.V_ohm, df.V_cell)
    assert np.all(np.diff(df.V_cell) > 0)


@pytest.mark.parametrize("activation", ["tafel", "asinh"])
def test_inverse_current_density(activation):
    m = am.AECModel(P=7, w=0.3, activation=activation)
    j = np.array([20.0, 150.0, 400.0])
    V = m.voltage(j, 333.15)
    np.testing.assert_allclose(m.current_density(V, 333.15), j, rtol=1e-6)


def test_asinh_and_tafel_agree_at_high_current():
    t = am.AECModel(activation="tafel")
    a = am.AECModel(activation="asinh")
    assert t.voltage(300, 343.15) == pytest.approx(a.voltage(300, 343.15), abs=1e-6)
    assert a.activation_overpotential(1e-9, 343.15) >= 0
