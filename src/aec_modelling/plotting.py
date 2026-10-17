"""Optional matplotlib helpers (install with ``pip install aec-modelling[plot]``)."""

from __future__ import annotations

import numpy as np


def _plt():
    try:
        import matplotlib.pyplot as plt
    except ImportError as e:  # pragma: no cover
        raise ImportError("Plotting requires matplotlib: pip install 'aec-modelling[plot]'") from e
    return plt


def plot_fit(fit, data, calibration=None, ax=None, j=None, title=None):
    """Measured points vs. fitted polarisation curves, one colour per temperature.

    Parameters
    ----------
    fit : FitResult
    data : DataFrame
        Full dataset (validation points, drawn as circles).
    calibration : DataFrame, optional
        Subset used for calibration (drawn as squares).
    """
    plt = _plt()
    ax = ax or plt.subplots(figsize=(6, 4.5))[1]
    j = np.linspace(max(1.0, data["j"].min() * 0.5), data["j"].max() * 1.1, 200) if j is None else j
    P = float(data["P"].median()) if "P" in data else None
    w = float(data["w"].median()) if "w" in data else None
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    for i, (T, grp) in enumerate(data.groupby("T")):
        c = colors[i % len(colors)]
        ax.scatter(grp["j"], grp["V"], color=c, s=25, alpha=0.7)
        ax.plot(j, fit.model.voltage(j, T, P, w), color=c, label=f"{T - 273.15:.1f} °C")
    if calibration is not None:
        ax.scatter(
            calibration["j"],
            calibration["V"],
            marker="s",
            s=60,
            facecolors="none",
            edgecolors="k",
            label="calibration points",
        )
    ax.set_xlabel("Current density $j$ [mA/cm$^2$]")
    ax.set_ylabel("Cell voltage $V_{cell}$ [V]")
    ax.set_title(title or f"R$^2$ = {fit.evaluate(data)['r2'].iloc[0]:.3f}")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, frameon=False)
    return ax


def plot_overpotentials(model, j, T, ax=None, P=None, w=None):
    """Stacked area plot of reversible, activation and ohmic contributions."""
    plt = _plt()
    ax = ax or plt.subplots(figsize=(6, 4.5))[1]
    df = model.overpotentials(j, T, P, w)
    ax.stackplot(
        df["j"],
        df["V_rev"],
        df["V_act"],
        df["V_ohm"],
        labels=["$V_{rev}$", "$V_{act}$", "$V_{ohm}$"],
        alpha=0.8,
    )
    ax.plot(df["j"], df["V_tn"], "k--", lw=1, label="$V_{tn}$")
    ax.set_xlabel("Current density $j$ [mA/cm$^2$]")
    ax.set_ylabel("Voltage [V]")
    ax.set_title(f"T = {T - 273.15:.0f} °C")
    ax.legend(loc="lower right", frameon=False)
    ax.grid(alpha=0.3)
    return ax
