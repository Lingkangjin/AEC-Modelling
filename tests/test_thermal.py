import numpy as np
import pytest

import aec_modelling as am


@pytest.fixture
def stack():
    cell = am.AECModel(delta_l=0.59, alpha=0.131, C1=-6646, C2=20.8, P=16, w=0.25)
    return am.Stack(cell, n_cells=50, cell_area=1.0e4)


def test_stack_scaling(stack):
    I = stack.current(300.0)
    assert stack.voltage(I, 353.15) == pytest.approx(50 * stack.cell.voltage(300.0, 353.15))
    # 1 A for 1 h in one cell -> 0.0376 g H2
    one = am.Stack(stack.cell, n_cells=1)
    assert one.hydrogen_production_kg_h(1.0) * 1000 == pytest.approx(0.0376, rel=1e-3)


def test_required_flow_gives_target_steady_state(stack):
    lm = am.LumpedThermalModel(stack, heat_capacity=1e6, thermal_resistance=0.05)
    I = stack.current(300.0)
    m = lm.required_lye_flow(I, T_out=353.15, T_in=343.15)
    assert m > 0
    assert lm.steady_state(I, 343.15, m) == pytest.approx(353.15, abs=1e-6)


def test_transient_converges_to_steady_state(stack):
    lm = am.LumpedThermalModel(stack, heat_capacity=1e6)
    I, m = stack.current(250.0), 2.0
    T_ss = lm.steady_state(I, 340.0, m)
    df = lm.simulate(np.linspace(0, 20000, 50), I=I, T_in=340.0, m_lye=m, T0=300.0)
    assert df["T"].iloc[-1] == pytest.approx(T_ss, abs=0.05)


def test_axial_model_energy_balance_and_current(stack):
    ax = am.AxialThermalModel(stack, n_segments=30, cp_lye=3200.0)
    I, m, T_in = stack.current(300.0), 3.0, 343.15
    r = ax.solve(I, T_in, m)
    assert r.profile["j"].mean() == pytest.approx(300.0, rel=1e-6)
    q = (r.profile["q"] * stack.cell_area / 30).sum() * stack.n_cells
    assert q == pytest.approx(m * 3200.0 * (r.T_out - T_in), rel=1e-9)
    # exothermic operation: lye heats up, downstream strips carry more current
    assert r.T_out > T_in
    assert r.profile["j"].iloc[-1] > r.profile["j"].iloc[0]
