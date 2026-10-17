"""Fit the four-parameter model to polarisation data.

The fitting works on a *long-format* table with one row per measured point:

======  ===========  ============================================
column  unit         meaning
======  ===========  ============================================
j       mA/cm^2      current density
T       K            temperature
V       V            measured cell voltage
P       bar (abs.)   system pressure  (optional if passed as kwarg)
w       -            KOH mass fraction (optional if passed as kwarg)
======  ===========  ============================================

All bundled datasets (:mod:`aec_modelling.datasets`) follow this layout.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

from .electrochemistry import PARAMETER_NAMES, AECModel, cell_voltage

#: Initial guess (delta_l [cm], alpha [-], C1 [K], C2 [-]) used in the paper.
DEFAULT_P0 = (0.40, 0.13, -6330.0, 14.0)
#: Default lower / upper bounds of the four parameters.
DEFAULT_BOUNDS = ((0.0, 0.05, -10000.0, -5.0), (2.0, 1.0, 0.0, 30.0))


# ----------------------------------------------------------------- metrics
def r2_score(y_true, y_pred) -> float:
    """Coefficient of determination R^2."""
    y_true, y_pred = np.asarray(y_true, float), np.asarray(y_pred, float)
    ss_res = np.sum((y_true - y_pred) ** 2)
    ss_tot = np.sum((y_true - y_true.mean()) ** 2)
    return float(1.0 - ss_res / ss_tot) if ss_tot > 0 else float("nan")


def rmse(y_true, y_pred) -> float:
    """Root-mean-square error (same unit as the data)."""
    return float(np.sqrt(np.mean((np.asarray(y_true, float) - np.asarray(y_pred, float)) ** 2)))


def mape(y_true, y_pred) -> float:
    """Mean absolute percentage error [%]."""
    y_true, y_pred = np.asarray(y_true, float), np.asarray(y_pred, float)
    return float(100.0 * np.mean(np.abs((y_true - y_pred) / y_true)))


def score(y_true, y_pred) -> dict:
    """R^2, RMSE [V] and MAPE [%] in one dict."""
    return {
        "r2": r2_score(y_true, y_pred),
        "rmse": rmse(y_true, y_pred),
        "mape": mape(y_true, y_pred),
    }


# --------------------------------------------------------------- data prep
def _columns(data: pd.DataFrame | None, j, T, V, P, w):
    if data is not None:
        missing = [c for c in ("j", "T", "V") if c not in data]
        if missing:
            raise KeyError(f"data is missing required column(s) {missing}")
        j, T, V = data["j"], data["T"], data["V"]
        if P is None and "P" in data:
            P = data["P"]
        if w is None and "w" in data:
            w = data["w"]
    if P is None or w is None:
        raise ValueError("Pressure P and KOH mass fraction w are required (as columns or kwargs)")
    j, T, V = (np.asarray(x, dtype=float) for x in (j, T, V))
    P = np.broadcast_to(np.asarray(P, dtype=float), j.shape)
    w = np.broadcast_to(np.asarray(w, dtype=float), j.shape)
    if np.isnan(P).any() or np.isnan(w).any():
        raise ValueError(
            "P or w contains NaN; pass the operating conditions explicitly, "
            "e.g. fit_polarization(df, w=0.30)"
        )
    return j, T, V, P, w


def sample_points(
    data: pd.DataFrame,
    n: int = 20,
    min_temperatures: int = 2,
    min_per_temperature: int = 2,
    random_state=None,
    max_tries: int = 1000,
):
    """Randomly draw ``n`` calibration points spanning several temperatures.

    This mirrors the procedure of the paper: the model is calibrated on a small
    random subset and validated on the full dataset, which tests its ability
    to extrapolate in temperature.
    """
    rng = np.random.default_rng(random_state)
    if n >= len(data):
        return data.copy()
    for _ in range(max_tries):
        sub = data.sample(n=n, random_state=int(rng.integers(2**31)))
        counts = sub.groupby("T").size()
        if len(counts) >= min_temperatures and counts.min() >= min_per_temperature:
            return sub
    raise RuntimeError(
        "Could not draw a sample satisfying the temperature constraints; "
        "increase n or relax min_temperatures / min_per_temperature."
    )


# ------------------------------------------------------------------ result
@dataclass
class FitResult:
    """Outcome of :func:`fit_polarization`."""

    model: AECModel
    covariance: np.ndarray
    n_points: int
    metrics: dict
    stderr: dict = field(init=False)

    def __post_init__(self):
        sd = np.sqrt(np.clip(np.diag(self.covariance), 0, None))
        self.stderr = dict(zip(PARAMETER_NAMES, sd))

    @property
    def params(self) -> dict:
        return self.model.params

    @property
    def r2(self) -> float:
        return self.metrics["r2"]

    @property
    def rmse(self) -> float:
        return self.metrics["rmse"]

    def predict(self, data: pd.DataFrame) -> np.ndarray:
        """Predicted cell voltage for every row of ``data`` (uses its P/w columns if present)."""
        P = data["P"] if "P" in data else self.model.P
        w = data["w"] if "w" in data else self.model.w
        return self.model.voltage(
            data["j"].to_numpy(), data["T"].to_numpy(), np.asarray(P), np.asarray(w)
        )

    def evaluate(self, data: pd.DataFrame, by: str | None = None) -> pd.DataFrame:
        """Goodness of fit on ``data``; ``by="T"`` returns one row per temperature."""
        pred = pd.Series(self.predict(data), index=data.index)
        if by is None:
            return pd.DataFrame([{"n": len(data), **score(data["V"], pred)}])
        rows = []
        for key, grp in data.groupby(by):
            rows.append({by: key, "n": len(grp), **score(grp["V"], pred.loc[grp.index])})
        return pd.DataFrame(rows)

    def summary(self) -> str:
        units = {"delta_l": "cm", "alpha": "-", "C1": "K", "C2": "-"}
        lines = ["Four-parameter AEC model fit", "-" * 44]
        for k in PARAMETER_NAMES:
            lines.append(
                f"{k:>8} = {self.params[k]:>12.5g} +/- {self.stderr[k]:<10.3g} [{units[k]}]"
            )
        lines += [
            "-" * 44,
            f"points = {self.n_points}   R2 = {self.r2:.4f}   RMSE = {1000 * self.rmse:.1f} mV",
            f"activation energy Ea = {self.model.activation_energy() / 1000:.1f} kJ/mol",
        ]
        return "\n".join(lines)

    def __str__(self):  # pragma: no cover - cosmetic
        return self.summary()


# --------------------------------------------------------------------- fit
def fit_polarization(
    data: pd.DataFrame | None = None,
    *,
    j=None,
    T=None,
    V=None,
    P=None,
    w=None,
    p0=DEFAULT_P0,
    bounds=DEFAULT_BOUNDS,
    activation: str = "tafel",
    **curve_fit_kwargs,
) -> FitResult:
    """Least-squares fit of ``delta_l, alpha, C1, C2`` to measured polarisation data.

    Pass either a DataFrame with columns ``j, T, V`` (+ optionally ``P, w``) or
    the arrays as keyword arguments. ``P`` and ``w`` can be scalars or arrays
    (e.g. a dataset measured at several pressures).

    Returns
    -------
    FitResult
        Holds the calibrated :class:`AECModel`, parameter standard errors and
        R^2/RMSE on the calibration points.
    """
    j, T, V, P, w = _columns(data, j, T, V, P, w)
    if np.any(j <= 0):
        raise ValueError("All current densities must be > 0 mA/cm^2 (drop open-circuit points).")

    def f(X, delta_l, alpha, C1, C2):
        jj, TT, PP, ww = X
        return cell_voltage(jj, TT, PP, ww, delta_l, alpha, C1, C2, activation=activation)

    curve_fit_kwargs.setdefault("maxfev", 20000)
    popt, pcov = curve_fit(f, (j, T, P, w), V, p0=p0, bounds=bounds, **curve_fit_kwargs)

    model = AECModel(
        *(float(p) for p in popt),
        P=float(np.median(P)),
        w=float(np.median(w)),
        activation=activation,
    )
    pred = f((j, T, P, w), *popt)
    return FitResult(model=model, covariance=pcov, n_points=len(j), metrics=score(V, pred))
