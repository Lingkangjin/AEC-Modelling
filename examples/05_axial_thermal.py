"""1D thermal model: temperature rise along the cell and current maldistribution.

Lower lye flow -> larger temperature rise -> downstream regions draw more
current. The 0D model at the mean temperature is shown for comparison.

Run:  python examples/05_axial_thermal.py
"""

from pathlib import Path

import matplotlib.pyplot as plt

from aec_modelling import AxialThermalModel, Stack, datasets, fit_polarization

OUT = Path(__file__).parent / "figures"
OUT.mkdir(exist_ok=True)

cell = fit_polarization(datasets.load("sakas")).model
stack = Stack(cell, n_cells=100, cell_area=2.6e4)
model = AxialThermalModel(stack, n_segments=50)

I = stack.current(300.0)
T_in = 343.15

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), layout="constrained")
print(f"{'m_lye [kg/s]':>12} {'T_out [C]':>10} {'V_1D [V]':>9} {'V_0D(Tmean)':>12} {'maldist.':>9}")
for m_lye in (6.0, 8.0, 12.0, 16.0):
    r = model.solve(I, T_in, m_lye)
    v0d = float(stack.cell_voltage(I, r.T_mean))
    print(
        f"{m_lye:>12.1f} {r.T_out - 273.15:>10.2f} {r.V_cell:>9.4f} {v0d:>12.4f} "
        f"{100 * r.current_maldistribution:>8.1f}%"
    )
    axes[0].plot(r.profile["z"], r.profile["T"] - 273.15, label=f"{m_lye:g} kg/s")
    axes[1].plot(r.profile["z"], r.profile["j"], label=f"{m_lye:g} kg/s")

axes[0].set(xlabel="normalised position along lye flow [-]", ylabel="T [°C]")
axes[1].set(xlabel="normalised position along lye flow [-]", ylabel="j [mA/cm$^2$]")
for ax in axes:
    ax.legend(title="lye flow", frameon=False)
    ax.grid(alpha=0.3)

fig.savefig(OUT / "05_axial_thermal.png", dpi=150)
plt.show()
