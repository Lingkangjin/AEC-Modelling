"""Calibrate the four parameters on a few points and validate on the full dataset.

This is the workflow of the original ``Main.py``: draw 20 random points that
span at least two temperatures, fit ``delta_l, alpha, C1, C2`` on them and
check how well the model reproduces every measured point.

Run:  python examples/02_fit_sakas.py
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from aec_modelling import datasets, fit_polarization, sample_points
from aec_modelling.plotting import plot_fit, plot_overpotentials

OUT = Path(__file__).parent / "figures"
OUT.mkdir(exist_ok=True)

df = datasets.load("sakas")  # columns: j [mA/cm2], T [K], V [V], P [bar], w [-]
print(df.attrs["reference"], "| P =", df.attrs["P"], "bar | w =", df.attrs["w"])

calibration = sample_points(df, n=20, min_temperatures=2, min_per_temperature=2, random_state=42)

fit = fit_polarization(
    calibration,
    p0=(0.40, 0.13, -6330, 14),  # delta_l, alpha, C1, C2
    bounds=((0, 0.1, -8000, 0), (1, 0.4, -3000, 25)),
)
print(fit.summary())

print("\nValidation on all points, per temperature:")
print(fit.evaluate(df, by="T").round(4).to_string(index=False))

fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), layout="constrained")
plot_fit(fit, df, calibration=calibration, ax=axes[0], j=np.linspace(50, 400, 100))
axes[0].set_title("Sakas et al. (2022): calibration vs. validation")
plot_overpotentials(fit.model, np.linspace(50, 400, 100), 353.15, ax=axes[1])

fig.savefig(OUT / "02_fit_sakas.png", dpi=150)
plt.show()
