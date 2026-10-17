"""Thermal models of an alkaline electrolysis stack coupled to the electrochemical model.

Heat is released whenever the cell voltage exceeds the thermoneutral voltage::

    Q_gen = n_cells * (V_cell(j, T) - V_tn(T)) * I                         [W]

Two levels of detail are provided:

* :class:`LumpedThermalModel` (0D) -- the stack is a single, well-mixed thermal
  mass. Lye leaves at the stack temperature. Supports transient simulation,
  steady state and lye-flow sizing::

      C_th dT/dt = Q_gen - (T - T_amb)/R_th - m_lye c_p (T - T_in)

* :class:`AxialThermalModel` (1D) -- steady-state, adiabatic cell discretised
  along the lye flow direction. Each segment sees the same (equipotential) cell
  voltage, so warmer segments downstream draw more current. Captures the
  temperature rise across the stack and the resulting current maldistribution.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from . import electrolyte
from .constants import M_H2, Z_ELECTRONS, F
from .electrochemistry import AECModel
from .thermo import thermoneutral_voltage


def _as_function(x):
    return x if callable(x) else (lambda t, _x=x: _x)


@dataclass
class Stack:
    """A stack of identical cells in series.

    Parameters
    ----------
    cell : AECModel
        Calibrated cell model (its ``P`` and ``w`` are the operating conditions).
    n_cells : int
        Number of cells in series.
    cell_area : float
        Active area per cell [cm^2].
    faradaic_efficiency : float
        Fraction of current producing hydrogen [-].
    """

    cell: AECModel
    n_cells: int = 1
    cell_area: float = 100.0
    faradaic_efficiency: float = 1.0

    def current_density(self, I):
        """Stack current [A] -> current density [mA/cm^2]."""
        return np.asarray(I, dtype=float) * 1000.0 / self.cell_area

    def current(self, j):
        """Current density [mA/cm^2] -> stack current [A]."""
        return np.asarray(j, dtype=float) * self.cell_area / 1000.0

    def cell_voltage(self, I, T):
        return self.cell.voltage(self.current_density(I), T)

    def voltage(self, I, T):
        """Stack voltage [V]."""
        return self.n_cells * self.cell_voltage(I, T)

    def power(self, I, T):
        """Electrical power consumption [W]."""
        return self.voltage(I, T) * np.asarray(I, dtype=float)

    def heat_generation(self, I, T):
        """Net heat released by the stack [W] (negative = endothermic)."""
        return self.n_cells * (self.cell_voltage(I, T) - thermoneutral_voltage(T)) * np.asarray(I)

    def hydrogen_production(self, I):
        """Hydrogen production [mol/s]."""
        return (
            self.faradaic_efficiency * self.n_cells * np.asarray(I, dtype=float) / (Z_ELECTRONS * F)
        )

    def hydrogen_production_kg_h(self, I):
        """Hydrogen production [kg/h]."""
        return self.hydrogen_production(I) * M_H2 / 1000.0 * 3600.0

    def efficiency(self, I, T):
        """HHV efficiency of the stack (voltage x Faradaic efficiency) [-]."""
        return self.faradaic_efficiency * thermoneutral_voltage(T) / self.cell_voltage(I, T)


# ---------------------------------------------------------------------- 0D
@dataclass
class LumpedThermalModel:
    """0D (well-mixed) energy balance of the stack.

    Parameters
    ----------
    stack : Stack
    heat_capacity : float
        Lumped heat capacity of stack + lye hold-up [J/K].
    thermal_resistance : float
        Resistance between stack and ambient [K/W]; ``inf`` = adiabatic.
    T_amb : float
        Ambient temperature [K].
    cp_lye : float, optional
        Lye specific heat [J/(kg K)]; defaults to :func:`electrolyte.specific_heat`.
    """

    stack: Stack
    heat_capacity: float = 5.0e5
    thermal_resistance: float = np.inf
    T_amb: float = 298.15
    cp_lye: float | None = None

    def _cp(self, T):
        return (
            self.cp_lye
            if self.cp_lye is not None
            else electrolyte.specific_heat(T, self.stack.cell.w)
        )

    def heat_flows(self, T, I, T_in, m_lye):
        """Return the terms of the energy balance [W] at stack temperature ``T``."""
        q_gen = self.stack.heat_generation(I, T)
        q_loss = (T - self.T_amb) / self.thermal_resistance
        q_lye = m_lye * self._cp(T) * (T - T_in)
        return {
            "Q_gen": q_gen,
            "Q_loss": q_loss,
            "Q_lye": q_lye,
            "dTdt": (q_gen - q_loss - q_lye) / self.heat_capacity,
        }

    def steady_state(self, I, T_in, m_lye, T_bounds=(250.0, 500.0)):
        """Steady-state stack (= lye outlet) temperature [K]."""

        def balance(T):
            hf = self.heat_flows(T, I, T_in, m_lye)
            return hf["Q_gen"] - hf["Q_loss"] - hf["Q_lye"]

        lo, hi = T_bounds
        if balance(lo) * balance(hi) > 0:
            raise ValueError(
                "No steady state within T_bounds; the stack needs more cooling "
                "(increase m_lye or decrease the current)."
            )
        return brentq(balance, lo, hi, xtol=1e-8)

    def required_lye_flow(self, I, T_out, T_in):
        """Lye mass flow [kg/s] that keeps the stack at ``T_out`` with inlet ``T_in``."""
        if T_out <= T_in:
            raise ValueError("T_out must be greater than T_in")
        hf = self.heat_flows(T_out, I, T_in, 0.0)
        return (hf["Q_gen"] - hf["Q_loss"]) / (self._cp(T_out) * (T_out - T_in))

    def simulate(self, t, I, T_in, m_lye, T0=None, **ivp_kwargs) -> pd.DataFrame:
        """Transient simulation.

        Parameters
        ----------
        t : array_like
            Output times [s] (first value = start time).
        I, T_in, m_lye : float or callable ``f(t)``
            Stack current [A], lye inlet temperature [K], lye flow [kg/s].
        T0 : float, optional
            Initial stack temperature [K] (default: ``T_in(t0)``).

        Returns
        -------
        DataFrame with time series of temperature, voltage, heat flows and H2 production.
        """
        t = np.asarray(t, dtype=float)
        I_f, Tin_f, m_f = _as_function(I), _as_function(T_in), _as_function(m_lye)
        T0 = Tin_f(t[0]) if T0 is None else T0

        def rhs(tt, y):
            return [self.heat_flows(y[0], I_f(tt), Tin_f(tt), m_f(tt))["dTdt"]]

        ivp_kwargs.setdefault("max_step", max((t[-1] - t[0]) / 500.0, 1e-3))
        ivp_kwargs.setdefault("method", "LSODA")
        sol = solve_ivp(rhs, (t[0], t[-1]), [T0], t_eval=t, **ivp_kwargs)
        if not sol.success:
            raise RuntimeError(sol.message)

        T = sol.y[0]
        I_t = np.array([I_f(tt) for tt in t], dtype=float)
        Tin_t = np.array([Tin_f(tt) for tt in t], dtype=float)
        m_t = np.array([m_f(tt) for tt in t], dtype=float)
        hf = self.heat_flows(T, I_t, Tin_t, m_t)
        return pd.DataFrame(
            {
                "t": t,
                "T": T,
                "T_in": Tin_t,
                "I": I_t,
                "m_lye": m_t,
                "j": self.stack.current_density(I_t),
                "V_cell": self.stack.cell_voltage(I_t, T),
                "P_el": self.stack.power(I_t, T),
                "Q_gen": hf["Q_gen"],
                "Q_loss": hf["Q_loss"] * np.ones_like(T),
                "Q_lye": hf["Q_lye"],
                "H2_kg_h": self.stack.hydrogen_production_kg_h(I_t),
                "efficiency": self.stack.efficiency(I_t, T),
            }
        )


# ---------------------------------------------------------------------- 1D
@dataclass
class AxialResult:
    """Solution of :meth:`AxialThermalModel.solve`."""

    V_cell: float
    T_in: float
    T_out: float
    j_mean: float
    profile: pd.DataFrame

    @property
    def T_mean(self) -> float:
        return float(self.profile["T"].mean())

    @property
    def current_maldistribution(self) -> float:
        """(j_max - j_min) / j_mean [-]."""
        j = self.profile["j"]
        return float((j.max() - j.min()) / self.j_mean)


@dataclass
class AxialThermalModel:
    """1D steady-state model along the lye flow direction (adiabatic cells).

    Each cell is split into ``n_segments`` strips perpendicular to the flow.
    All strips share the cell voltage; the local current density follows from
    the local temperature, and the heat released in a strip warms the lye that
    flows to the next one.

    Parameters
    ----------
    stack : Stack
    n_segments : int
        Number of strips per cell (discretisation is first-order upwind).
    cp_lye : float, optional
        Lye specific heat [J/(kg K)].
    """

    stack: Stack
    n_segments: int = 40
    cp_lye: float | None = None

    def _march(self, V, T_in, m_cell):
        cell, n = self.stack.cell, self.n_segments
        a_seg = self.stack.cell_area / n
        T = np.empty(n + 1)
        j = np.empty(n)
        T[0] = T_in
        for k in range(n):
            j[k] = cell._invert(V, T[k], cell.P, cell.w, 1e5)
            q = (V - thermoneutral_voltage(T[k])) * j[k] * a_seg / 1000.0  # W
            cp = self.cp_lye if self.cp_lye is not None else electrolyte.specific_heat(T[k], cell.w)
            T[k + 1] = T[k] + q / (m_cell * cp)
        return T, j

    def solve(self, I, T_in, m_lye) -> AxialResult:
        """Solve for the cell voltage and the temperature / current profiles.

        Parameters
        ----------
        I : float
            Stack current [A].
        T_in : float
            Lye inlet temperature [K].
        m_lye : float
            Total lye flow to the stack [kg/s] (split equally between cells).
        """
        cell = self.stack.cell
        m_cell = m_lye / self.stack.n_cells
        j_mean = float(self.stack.current_density(I))

        def mismatch(V):
            _, j = self._march(V, T_in, m_cell)
            return j.mean() - j_mean

        lo = float(cell.reversible_voltage(T_in)) - 0.3
        hi = float(cell.voltage(j_mean, T_in))
        while mismatch(hi) < 0:
            hi += 0.1
            if hi > 10:
                raise RuntimeError("Could not bracket the cell voltage")
        V = brentq(mismatch, lo, hi, xtol=1e-9)

        T, j = self._march(V, T_in, m_cell)
        z = (np.arange(self.n_segments) + 0.5) / self.n_segments
        T_mid = 0.5 * (T[:-1] + T[1:])
        profile = pd.DataFrame(
            {
                "z": z,
                "T": T_mid,
                "j": j,
                "V_tn": thermoneutral_voltage(T[:-1]),
                "q": (V - thermoneutral_voltage(T[:-1])) * j / 1000.0,  # W/cm^2
            }
        )
        return AxialResult(V_cell=V, T_in=T_in, T_out=float(T[-1]), j_mean=j_mean, profile=profile)
