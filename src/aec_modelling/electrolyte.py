"""Properties of aqueous KOH electrolyte.

All functions are vectorised: ``T`` is temperature in kelvin and ``w`` is the
KOH mass fraction (0-1, e.g. 0.30 for 30 wt%).

References
----------
Gilliam RJ et al., Int J Hydrogen Energy 32 (2007) 359-364.
    Density, molarity and ionic conductivity correlations.
Balej J., Int J Hydrogen Energy 10 (1985) 233-243 / LeRoy RL et al. (1980).
    Vapour pressure and water activity of KOH solutions.
"""

from __future__ import annotations

import numpy as np

from .constants import KELVIN, M_KOH


def density(T, w):
    """Electrolyte density [kg/m^3] (Gilliam et al. 2007)."""
    t = np.asarray(T, dtype=float) - KELVIN
    rho_water = -0.0032 * t**2 - 0.1108 * t + 1001.7
    return rho_water * np.exp(0.86 * np.asarray(w, dtype=float))


def molarity(T, w):
    """KOH molar concentration [mol/L]."""
    return np.asarray(w, dtype=float) * density(T, w) / M_KOH


def molality(w):
    """KOH molality [mol per kg of water]."""
    w = np.asarray(w, dtype=float)
    return 1000.0 * w / (M_KOH * (1.0 - w))


def conductivity(T, w):
    """Specific ionic conductivity of the electrolyte [S/cm] (Gilliam et al. 2007).

    Valid roughly for 0-12 mol/L and 0-100 degC.
    """
    T = np.asarray(T, dtype=float)
    M = molarity(T, w)
    return (
        -2.041 * M
        - 0.0028 * M**2
        + 0.005332 * M * T
        + 207.2 * M / T
        + 0.001043 * M**3
        - 0.0000003 * M**2 * T**2
    )


def vapour_pressure_water(T):
    """Saturation pressure of pure water [bar]."""
    T = np.asarray(T, dtype=float)
    return np.exp(81.6179 - 7699.68 / T - 10.9 * np.log(T) + 9.5891e-3 * T)


def vapour_pressure(T, w):
    """Water vapour pressure above the KOH solution [bar]."""
    m = molality(w)
    a = -0.0151 * m - 1.6788e-3 * m**2 + 2.2588e-5 * m**3
    b = 1 - 1.2062e-3 * m + 5.6024e-4 * m**2 - 7.8228e-6 * m**3
    return np.exp(2.302 * a + b * np.log(vapour_pressure_water(T)))


def water_activity(T, w):
    """Activity of water in the KOH solution [-]."""
    T = np.asarray(T, dtype=float)
    m = molality(w)
    return np.exp(-0.05192 * m + 0.003302 * m**2 + (3.177 * m - 2.131 * m**2) / T)


def specific_heat(T, w):
    """Approximate specific heat capacity of the lye [J/(kg K)].

    Linear blend between water (4180 J/kg/K) and the value of 30 wt% KOH
    (~3100 J/kg/K) commonly used in system studies. Adequate for energy
    balances; replace with a dedicated correlation for detailed design.
    """
    w = np.asarray(w, dtype=float)
    return (4180.0 - 3600.0 * w) * np.ones_like(np.asarray(T, dtype=float))
