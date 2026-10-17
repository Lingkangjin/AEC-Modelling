"""Four-parameter semi-empirical polarisation model of an alkaline electrolysis cell.

The cell voltage is the sum of three contributions::

    V(j, T) = V_rev(T, P, w) + V_ohm(j, T, w) + V_act(j, T)

    V_rev = U0(T) + RT/(2F) * ln[(P - p_v)^1.5 / a_w]           reversible
    V_ohm = j * delta_l / sigma(T, w)                            ohmic
    V_act = RT/(2 alpha F) * ln(j / i0(T)),  ln i0 = C1/T + C2    activation (Tafel)

Only four parameters are fitted to polarisation data:

=========  =======  ==========================================================
Parameter  Unit     Meaning
=========  =======  ==========================================================
delta_l    cm       effective ionic path length (electrolyte + separator + bubbles)
alpha      -        lumped charge-transfer coefficient (anode + cathode)
C1         K        slope of the Arrhenius law of i0  (C1 = -Ea / R)
C2         -        intercept of the Arrhenius law of i0 (j in mA/cm^2)
=========  =======  ==========================================================

Units: current density ``j`` in mA/cm^2, temperature ``T`` in K, pressure ``P``
in bar (absolute), KOH mass fraction ``w`` in 0-1, voltages in V.

Reference: Jin L., Nogueira Nakashima R., Comodi G., Lund Frandsen H.,
Applied Thermal Engineering 262 (2025) 125154.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace

import numpy as np
import pandas as pd
from scipy.optimize import brentq

from . import electrolyte
from .constants import M_H2, Z_ELECTRONS, F, R
from .thermo import thermoneutral_voltage

PARAMETER_NAMES = ("delta_l", "alpha", "C1", "C2")
ACTIVATION_MODES = ("tafel", "asinh")


def empirical_reversible_voltage(T):
    """Reversible voltage at standard pressure, empirical fit [V] (LeRoy et al.)."""
    T = np.asarray(T, dtype=float)
    return 1.5184 - 1.5421e-3 * T + 9.523e-5 * T * np.log(T) + 9.84e-8 * T**2


def reversible_voltage(T, P, w):
    """Reversible cell voltage including pressure and electrolyte effects [V].

    Parameters
    ----------
    T : temperature [K]
    P : system pressure [bar, absolute]
    w : KOH mass fraction [-]
    """
    T = np.asarray(T, dtype=float)
    p_v = electrolyte.vapour_pressure(T, w)
    a_w = electrolyte.water_activity(T, w)
    p_gas = np.maximum(np.asarray(P, dtype=float) - p_v, 1e-6)
    return empirical_reversible_voltage(T) + R * T / (Z_ELECTRONS * F) * np.log(p_gas**1.5 / a_w)


def cell_voltage(j, T, P, w, delta_l, alpha, C1, C2, activation="tafel"):
    """Functional form of the four-parameter model (used by the fitting routines).

    All arguments broadcast against each other; returns the cell voltage [V].
    """
    rev, ohm, act = _contributions(j, T, P, w, delta_l, alpha, C1, C2, activation)
    return rev + ohm + act


def _contributions(j, T, P, w, delta_l, alpha, C1, C2, activation):
    j = np.asarray(j, dtype=float)
    T = np.asarray(T, dtype=float)
    rev = reversible_voltage(T, P, w)
    ohm = j * delta_l / (electrolyte.conductivity(T, w) * 1000.0)  # mA -> A
    b = R * T / (Z_ELECTRONS * alpha * F)
    ln_i0 = C1 / T + C2
    if activation == "tafel":
        act = b * (np.log(j) - ln_i0)
    elif activation == "asinh":
        act = b * np.arcsinh(j / (2.0 * np.exp(ln_i0)))
    else:
        raise ValueError(f"activation must be one of {ACTIVATION_MODES}, got {activation!r}")
    return rev, ohm, act


@dataclass
class AECModel:
    """Four-parameter alkaline electrolysis cell model.

    Parameters
    ----------
    delta_l, alpha, C1, C2
        The four fitted parameters (see module docstring).
    P : float
        Default system pressure [bar, absolute].
    w : float
        Default KOH mass fraction [-].
    activation : {"tafel", "asinh"}
        ``"tafel"`` reproduces the published model. ``"asinh"`` (Butler-Volmer
        symmetric form) is identical at practical current densities but stays
        positive and finite as ``j -> 0``.

    Examples
    --------
    >>> m = AECModel(delta_l=0.40, alpha=0.13, C1=-6330, C2=14, P=16, w=0.25)
    >>> float(m.voltage(200, 343.15))  # doctest: +SKIP
    1.93
    """

    delta_l: float = 0.40
    alpha: float = 0.13
    C1: float = -6330.0
    C2: float = 14.0
    P: float = 1.01325
    w: float = 0.30
    activation: str = "tafel"

    def __post_init__(self):
        if self.activation not in ACTIVATION_MODES:
            raise ValueError(f"activation must be one of {ACTIVATION_MODES}")

    # ------------------------------------------------------------------ params
    @property
    def params(self) -> dict:
        """The four fitted parameters as a dict."""
        return {k: getattr(self, k) for k in PARAMETER_NAMES}

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> AECModel:
        return cls(**d)

    def with_conditions(self, *, P=None, w=None) -> AECModel:
        """Return a copy with different default operating conditions."""
        return replace(self, P=self.P if P is None else P, w=self.w if w is None else w)

    def _pw(self, P, w):
        return (self.P if P is None else P), (self.w if w is None else w)

    # ---------------------------------------------------------------- voltages
    def exchange_current_density(self, T):
        """Lumped exchange current density i0(T) [mA/cm^2]."""
        return np.exp(self.C1 / np.asarray(T, dtype=float) + self.C2)

    def tafel_slope(self, T):
        """Tafel slope per decade [V/dec]."""
        return np.log(10) * R * np.asarray(T, dtype=float) / (Z_ELECTRONS * self.alpha * F)

    def activation_energy(self):
        """Apparent activation energy of i0 [J/mol] (``-C1 * R``)."""
        return -self.C1 * R

    def reversible_voltage(self, T, P=None, w=None):
        P, w = self._pw(P, w)
        return reversible_voltage(T, P, w)

    def ohmic_overpotential(self, j, T, w=None):
        w = self.w if w is None else w
        return np.asarray(j, dtype=float) * self.delta_l / (electrolyte.conductivity(T, w) * 1000.0)

    def activation_overpotential(self, j, T):
        return _contributions(j, T, self.P, self.w, *self.params.values(), self.activation)[2]

    def voltage(self, j, T, P=None, w=None):
        """Cell voltage [V] at current density ``j`` [mA/cm^2] and temperature ``T`` [K]."""
        P, w = self._pw(P, w)
        return cell_voltage(j, T, P, w, *self.params.values(), activation=self.activation)

    def overpotentials(self, j, T, P=None, w=None) -> pd.DataFrame:
        """Breakdown of the cell voltage into its contributions."""
        P, w = self._pw(P, w)
        j, T = np.broadcast_arrays(np.asarray(j, dtype=float), np.asarray(T, dtype=float))
        rev, ohm, act = _contributions(j, T, P, w, *self.params.values(), self.activation)
        rev = np.broadcast_to(rev, j.shape)
        return pd.DataFrame(
            {
                "j": j.ravel(),
                "T": T.ravel(),
                "V_rev": rev.ravel(),
                "V_act": act.ravel(),
                "V_ohm": ohm.ravel(),
                "V_cell": (rev + act + ohm).ravel(),
                "V_tn": np.broadcast_to(thermoneutral_voltage(T), j.shape).ravel(),
            }
        )

    def polarization_curve(self, j, temperatures, P=None, w=None) -> pd.DataFrame:
        """Polarisation curves for one or more temperatures (long-format DataFrame)."""
        frames = [self.overpotentials(j, T, P, w) for T in np.atleast_1d(temperatures)]
        return pd.concat(frames, ignore_index=True)

    def current_density(self, V, T, P=None, w=None, j_max=1e5):
        """Invert the model: current density [mA/cm^2] that gives voltage ``V``."""
        P, w = self._pw(P, w)
        V, T, P, w = np.broadcast_arrays(*(np.asarray(x, dtype=float) for x in (V, T, P, w)))
        out = np.empty(V.shape)
        for idx in np.ndindex(V.shape):
            out[idx] = self._invert(V[idx], T[idx], P[idx], w[idx], j_max)
        return out if out.ndim else float(out)

    def _invert(self, V, T, P, w, j_max):
        j_min = 1e-9
        f = lambda j: self.voltage(j, T, P, w) - V  # noqa: E731
        if f(j_min) >= 0:
            return 0.0
        if f(j_max) <= 0:
            return np.inf
        return brentq(f, j_min, j_max, xtol=1e-10, rtol=1e-12)

    # ------------------------------------------------------------- performance
    def voltage_efficiency(self, j, T, P=None, w=None):
        """HHV voltage efficiency ``V_tn / V_cell`` [-]."""
        return thermoneutral_voltage(T) / self.voltage(j, T, P, w)

    def heat_generation(self, j, T, P=None, w=None):
        """Net heat released per unit area ``(V - V_tn) j`` [W/cm^2].

        Positive values mean the cell heats up (exothermic operation).
        """
        return (self.voltage(j, T, P, w) - thermoneutral_voltage(T)) * np.asarray(j) / 1000.0

    @staticmethod
    def hydrogen_flux(j, faradaic_efficiency=1.0):
        """Hydrogen production per unit area [mol/(s cm^2)]."""
        return faradaic_efficiency * np.asarray(j, dtype=float) / 1000.0 / (Z_ELECTRONS * F)

    @classmethod
    def hydrogen_mass_flux(cls, j, faradaic_efficiency=1.0):
        """Hydrogen production per unit area [kg/(s cm^2)]."""
        return cls.hydrogen_flux(j, faradaic_efficiency) * M_H2 / 1000.0
