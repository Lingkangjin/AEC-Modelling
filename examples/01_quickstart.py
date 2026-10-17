"""Quickstart: polarisation curves and voltage breakdown of a cell.

Run:  python examples/01_quickstart.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from aec_modelling import AECModel
from aec_modelling.plotting import plot_overpotentials

OUT = Path(__file__).parent / "figures"
OUT.mkdir(exist_ok=True)

# Parameters fitted on the Sakas et al. industrial stack (see 02_fit_sakas.py)
cell = AECModel(delta_l=0.593, alpha=0.131, C1=-6646, C2=20.8, P=16.0, w=0.25)

j = np.linspace(5, 500, 200)  # mA/cm^2

print(f"V_cell(300 mA/cm2, 70 degC) = {cell.voltage(300, 343.15):.3f} V")
print(f"i0(70 degC) = {cell.exchange_current_density(343.15):.3e} mA/cm2")
print(f"Tafel slope = {1000 * cell.tafel_slope(343.15):.0f} mV/dec")
print(f"Ea = {cell.activation_energy() / 1000:.1f} kJ/mol")
print(f"j at 1.9 V, 80 degC = {cell.current_density(1.9, 353.15):.1f} mA/cm2")
print()
print(cell.overpotentials([100, 300, 500], 353.15).round(4).to_string(index=False))

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), layout="constrained")

for T_c in (50, 60, 70, 80, 90):
    axes[0].plot(j, cell.voltage(j, T_c + 273.15), label=f"{T_c} °C")
axes[0].set(xlabel="j [mA/cm$^2$]", ylabel="$V_{cell}$ [V]", title="Effect of temperature")
axes[0].legend(frameon=False)
axes[0].grid(alpha=0.3)

plot_overpotentials(cell, j, 353.15, ax=axes[1])

fig.savefig(OUT / "01_quickstart.png", dpi=150)
print(f"\nFigure saved to {OUT / '01_quickstart.png'}")
plt.show()
