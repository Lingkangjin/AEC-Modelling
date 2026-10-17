"""0D thermal model: lye-flow sizing, warm-up and a load step of an industrial stack.

Run:  python examples/04_lumped_thermal.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from aec_modelling import LumpedThermalModel, Stack, datasets, fit_polarization

OUT = Path(__file__).parent / "figures"
OUT.mkdir(exist_ok=True)

# 1) calibrated cell model
cell = fit_polarization(datasets.load("sakas")).model

# 2) stack: 100 cells of 2.6 m2 (illustrative MW-scale stack)
stack = Stack(cell, n_cells=100, cell_area=2.6e4)
I_nom = stack.current(300.0)  # A at 300 mA/cm2
print(
    f"Nominal: I = {I_nom:.0f} A, P_el = {stack.power(I_nom, 353.15) / 1e6:.2f} MW, "
    f"H2 = {stack.hydrogen_production_kg_h(I_nom):.1f} kg/h"
)

# 3) thermal model: heat capacity of stack + lye hold-up, heat losses to ambient
thermal = LumpedThermalModel(stack, heat_capacity=2.0e7, thermal_resistance=0.01, T_amb=298.15)

T_in, T_out = 343.15, 353.15
m_lye = thermal.required_lye_flow(I_nom, T_out=T_out, T_in=T_in)
print(f"Lye flow to keep 70 -> 80 degC: {m_lye:.2f} kg/s")
print(
    "Heat balance at nominal point [kW]:",
    {
        k: round(float(v) / 1000, 1)
        for k, v in thermal.heat_flows(T_out, I_nom, T_in, m_lye).items()
        if k != "dTdt"
    },
)

# 4) transient: cold start at 40 degC, full load, then 50 % load after 2 h
t = np.linspace(0, 4 * 3600, 400)
load = lambda t: I_nom if t < 2 * 3600 else 0.5 * I_nom  # noqa: E731
sim = thermal.simulate(t, I=load, T_in=T_in, m_lye=m_lye, T0=313.15)

fig, axes = plt.subplots(3, 1, figsize=(8, 9), sharex=True, layout="constrained")
h = sim["t"] / 3600
axes[0].plot(h, sim["T"] - 273.15, label="stack / lye outlet")
axes[0].plot(h, sim["T_in"] - 273.15, "--", label="lye inlet")
axes[0].set_ylabel("T [°C]")
axes[1].plot(h, sim["V_cell"], label="$V_{cell}$")
axes[1].plot(h, sim["efficiency"], label="HHV efficiency")
axes[1].set_ylabel("[V] / [-]")
for k in ("Q_gen", "Q_lye", "Q_loss"):
    axes[2].plot(h, sim[k] / 1000, label=k)
axes[2].set(ylabel="Heat [kW]", xlabel="time [h]")
for ax in axes:
    ax.legend(frameon=False)
    ax.grid(alpha=0.3)

fig.savefig(OUT / "04_lumped_thermal.png", dpi=150)
plt.show()
