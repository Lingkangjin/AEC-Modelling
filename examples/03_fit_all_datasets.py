"""Fit every bundled literature dataset and compare the parameters.

Run:  python examples/03_fit_all_datasets.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from aec_modelling import datasets, fit_polarization
from aec_modelling.plotting import plot_fit

OUT = Path(__file__).parent / "figures"
OUT.mkdir(exist_ok=True)

print(datasets.list_datasets().drop(columns="description").to_string(), "\n")

names = ["ulleberg", "sakas", "sanchez", "degroot_pressure", "degroot_temperature"]
rows = []
fig, axes = plt.subplots(2, 3, figsize=(15, 9), layout="constrained")

for name, ax in zip(names, axes.flat):
    df = datasets.load(name)
    fit = fit_polarization(df)
    rows.append(
        {
            "dataset": name,
            **fit.params,
            "Ea [kJ/mol]": fit.model.activation_energy() / 1000,
            "R2": fit.r2,
            "RMSE [mV]": 1000 * fit.rmse,
        }
    )
    plot_fit(fit, df, ax=ax, title=f"{name}  (R$^2$={fit.r2:.3f})")

# Field data from the UNIVPM stack: noisy, 1-min logging incl. transients.
# KOH fraction was not logged, so it must be supplied.
field = datasets.load("univpm")
axes.flat[-1].scatter(field["j"], field["V"], c=field["T"] - 273.15, s=12, cmap="viridis")
axes.flat[-1].set(title="univpm field data (raw)", xlabel="j [mA/cm$^2$]", ylabel="$V_{cell}$ [V]")
axes.flat[-1].grid(alpha=0.3)

table = pd.DataFrame(rows).set_index("dataset")
print(table.round(4).to_string())
print("\nNote: degroot_pressure has one temperature only, so C1 and C2 are not identifiable.")

fig.savefig(OUT / "03_fit_all_datasets.png", dpi=150)
plt.show()
