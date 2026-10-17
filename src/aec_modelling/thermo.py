"""Thermodynamics of the water-splitting reaction  H2O(l) -> H2(g) + 1/2 O2(g).

Species properties use the NIST Shomate equations (valid for 298-1000 K for
the gases and 298-500 K for liquid water). Enthalpies include the standard
enthalpy of formation, so ``delta_h(298.15)`` is the higher heating value of
hydrogen (~285.8 kJ/mol).

All functions accept scalars or numpy arrays of temperature in kelvin.
"""

from __future__ import annotations

import numpy as np

from .constants import Z_ELECTRONS, F

# Shomate coefficients A..H (NIST Chemistry WebBook)
_SHOMATE = {
    "H2O(l)": (-203.6060, 1523.290, -3196.413, 2474.455, 3.855326, -256.5478, -488.7163, -285.8304),
    "H2": (33.066178, -11.363417, 11.432816, -2.772874, -0.158558, -9.980797, 172.707974, 0.0),
    "O2": (31.32234, -20.23531, 57.86644, -36.50624, -0.007374, -8.903471, 246.7945, 0.0),
}
SPECIES = tuple(_SHOMATE)


def species_properties(species: str, T):
    """Return ``(cp, h, s)`` of a species at temperature ``T`` [K].

    Units: cp [J/(mol K)], h [J/mol] (formation-referenced), s [J/(mol K)].
    ``species`` is one of ``"H2O(l)"``, ``"H2"``, ``"O2"``.
    """
    try:
        A, B, C, D, E, F_, G, _H = _SHOMATE[species]
    except KeyError:
        raise ValueError(f"Unknown species {species!r}; choose from {SPECIES}") from None
    t = np.asarray(T, dtype=float) / 1000.0
    cp = A + B * t + C * t**2 + D * t**3 + E / t**2
    # H - H298 + dHf  (the "- H" term of the Shomate form is omitted on purpose)
    h = (A * t + B * t**2 / 2 + C * t**3 / 3 + D * t**4 / 4 - E / t + F_) * 1000.0
    s = A * np.log(t) + B * t + C * t**2 / 2 + D * t**3 / 3 - E / (2 * t**2) + G
    return cp, h, s


def _reaction(T, index):
    w = species_properties("H2O(l)", T)[index]
    h2 = species_properties("H2", T)[index]
    o2 = species_properties("O2", T)[index]
    return h2 + 0.5 * o2 - w


def delta_h(T):
    """Reaction enthalpy of liquid-water splitting [J/mol]."""
    return _reaction(T, 1)


def delta_s(T):
    """Reaction entropy of liquid-water splitting [J/(mol K)]."""
    return _reaction(T, 2)


def delta_g(T):
    """Gibbs free energy of liquid-water splitting at 1 bar [J/mol]."""
    return delta_h(T) - np.asarray(T, dtype=float) * delta_s(T)


def thermoneutral_voltage(T):
    """Thermoneutral voltage ``dH / (zF)`` [V] (~1.48 V at 25 degC).

    A cell operated above this voltage produces net heat; below it, the
    reaction absorbs heat from the surroundings.
    """
    return delta_h(T) / (Z_ELECTRONS * F)


def standard_reversible_voltage(T):
    """Thermodynamic reversible voltage at 1 bar, ``dG / (zF)`` [V]."""
    return delta_g(T) / (Z_ELECTRONS * F)


class ThermalProperties:
    """Object wrapper kept for backwards compatibility with the original scripts.

    >>> ThermalProperties(353.15).delta_h()  # doctest: +SKIP
    """

    def __init__(self, T_kel):
        self.T_kel = T_kel

    def water(self):
        return species_properties("H2O(l)", self.T_kel)

    def hydrogen(self):
        return species_properties("H2", self.T_kel)

    def oxygen(self):
        return species_properties("O2", self.T_kel)

    def delta_h(self):
        return delta_h(self.T_kel)

    def delta_s(self):
        return delta_s(self.T_kel)
