"""aec-modelling: thermo-electrochemical modelling of alkaline water electrolysers.

Main entry points
-----------------
- :class:`AECModel`              four-parameter polarisation model of a cell
- :func:`fit_polarization`       calibrate the four parameters on measured data
- :class:`Stack`                 stack = cell model + number of cells + area
- :class:`LumpedThermalModel`    0D dynamic / steady-state stack energy balance
- :class:`AxialThermalModel`     1D steady-state temperature & current distribution
- :mod:`aec_modelling.datasets`  literature polarisation curves ready for fitting
- :mod:`aec_modelling.electrolyte`, :mod:`aec_modelling.thermo`  property functions
"""

from . import datasets, electrolyte, thermo
from .electrochemistry import AECModel, cell_voltage, reversible_voltage
from .fitting import FitResult, fit_polarization, r2_score, rmse, sample_points
from .thermal import AxialThermalModel, LumpedThermalModel, Stack
from .thermo import thermoneutral_voltage

__version__ = "0.1.0"

__all__ = [
    "AECModel",
    "AxialThermalModel",
    "FitResult",
    "LumpedThermalModel",
    "Stack",
    "cell_voltage",
    "datasets",
    "electrolyte",
    "fit_polarization",
    "r2_score",
    "reversible_voltage",
    "rmse",
    "sample_points",
    "thermo",
    "thermoneutral_voltage",
    "__version__",
]
