# Model description

This page documents the equations implemented in `aec_modelling`. Symbols and
units follow the code: current density **j in mA/cm²**, temperature **T in K**,
pressure **P in bar (absolute)**, KOH mass fraction **w in 0–1**.

## 1. Electrochemical model (`aec_modelling.electrochemistry`)

The cell voltage is split into three contributions:

$$
V_{cell}(j,T) = V_{rev}(T,P,w) + V_{act}(j,T) + V_{ohm}(j,T,w)
$$

### Reversible voltage

$$
V_{rev} = U_0(T) + \frac{RT}{2F}\ln\frac{(P - p_v)^{1.5}}{a_{w}}
$$

$$
U_0(T) = 1.5184 - 1.5421\cdot10^{-3}\,T + 9.523\cdot10^{-5}\,T\ln T + 9.84\cdot10^{-8}\,T^2
$$

with $p_v(T,m)$ the water vapour pressure above the lye and $a_w(T,m)$ the water
activity, both functions of the KOH molality $m$ (Balej 1985; LeRoy et al. 1980).
`thermo.standard_reversible_voltage` gives the purely thermodynamic
$\Delta G/2F$ from NIST Shomate data. At 1 bar it agrees with $U_0$ to within 1 mV.

### Activation overpotential (lumped anode + cathode)

$$
V_{act} = \frac{RT}{2\alpha F}\ln\frac{j}{i_0(T)}, \qquad \ln i_0 = \frac{C_1}{T} + C_2
$$

The Arrhenius form gives an apparent activation energy $E_a = -C_1 R$.
`AECModel(activation="asinh")` uses
$\frac{RT}{2\alpha F}\,\mathrm{asinh}\!\left(\frac{j}{2 i_0}\right)$ instead.
It matches the Tafel form at practical current densities and stays finite as $j \to 0$.

### Ohmic overpotential

$$
V_{ohm} = \frac{\delta_l}{\sigma_l(T,w)}\,\frac{j}{1000}
$$

$\sigma_l$ is the KOH conductivity from Gilliam et al. (2007), in S/cm.
$\delta_l$ is an *effective* ionic path length. It lumps the electrode gap, the
separator (via its MacMullin number) and the bubble void fraction. This is why
the fitted values of 0.4–0.6 cm are larger than the geometric gap.

### Fitted parameters

| Parameter | Unit | Typical range | Role |
|-----------|------|---------------|------|
| `delta_l` | cm   | 0.3 – 0.8     | slope of the linear (ohmic) region |
| `alpha`   | –    | 0.1 – 0.35    | Tafel slope $2.303\,RT/(2\alpha F)$ |
| `C1`      | K    | −7000 – −3000 | temperature dependence of $i_0$ |
| `C2`      | –    | 10 – 21       | magnitude of $i_0$ |

With data at a single temperature, `C1` and `C2` are perfectly correlated and
only $\ln i_0(T)$ is identifiable. Use at least two temperatures; the paper
shows that ~20 points spread over two temperatures are enough to extrapolate.

## 2. Electrolyte properties (`aec_modelling.electrolyte`)

| Function | Correlation |
|----------|-------------|
| `density(T, w)` | $\rho = (-0.0032\,t^2 - 0.1108\,t + 1001.7)\,e^{0.86w}$ [kg/m³], $t$ in °C |
| `molarity(T, w)` | $c = w\rho / M_{KOH}$ [mol/L] |
| `molality(w)` | $m = 1000\,w / (M_{KOH}(1-w))$ [mol/kg H₂O] |
| `conductivity(T, w)` | Gilliam et al. (2007) [S/cm] |
| `vapour_pressure(T, w)` | Balej (1985) [bar] |
| `water_activity(T, w)` | Balej (1985) [–] |
| `specific_heat(T, w)` | linear 4180 → 3100 J/(kg K) between 0 and 30 wt% |

## 3. Thermodynamics (`aec_modelling.thermo`)

The NIST Shomate equations for H₂O(l), H₂ and O₂ give $\Delta H(T)$, $\Delta S(T)$ and $\Delta G(T)$.
The **thermoneutral voltage** is

$$
V_{tn}(T) = \frac{\Delta H(T)}{2F} \approx 1.48\ \text{V}
$$

## 4. Thermal models (`aec_modelling.thermal`)

The net heat released by a stack of $N$ cells carrying current $I$ is

$$
\dot Q_{gen} = N\,(V_{cell} - V_{tn})\,I
$$

### 0D lumped model, `LumpedThermalModel`

The stack is treated as one well-mixed thermal mass, and the lye leaves at the stack temperature:

$$
C_{th}\frac{dT}{dt} = \dot Q_{gen}(T,I) - \frac{T - T_{amb}}{R_{th}} - \dot m_{lye}\,c_p\,(T - T_{in})
$$

* `simulate(t, I, T_in, m_lye)`: transient response. Each input can be a constant or a function of time.
* `steady_state(I, T_in, m_lye)`: the stack temperature at which the right-hand side is zero.
* `required_lye_flow(I, T_out, T_in)`: the lye flow needed to hold a target outlet temperature.

### 1D axial model, `AxialThermalModel`

This is a steady-state, adiabatic model of each cell, discretised into $n$ strips along the lye flow.
All strips share the cell voltage $V$ because the electrodes are equipotential.
For strip $k$:

$$
V = V_{cell}(j_k, T_k), \qquad
\dot m_{cell}\,c_p\,(T_{k+1}-T_k) = (V - V_{tn}(T_k))\,j_k A_k
$$

The solver finds $V$ such that $\frac{1}{n}\sum_k j_k = \bar j$. The temperature
rises along the cell, so downstream strips draw more current. The result
reports this as `current_maldistribution` $=(j_{max}-j_{min})/\bar j$. This is the
mechanism behind the paper's conclusion that a small lye temperature difference
(76–80 °C) improves efficiency.

## 5. Units cheat-sheet

| Quantity | Unit in API |
|----------|-------------|
| current density `j` | mA/cm² |
| current `I` | A |
| cell area | cm² |
| temperature | K |
| pressure | bar (absolute) |
| KOH fraction `w` | mass fraction 0–1 |
| lye flow `m_lye` | kg/s (whole stack) |
| heat flows | W |
| heat capacity | J/K |
| thermal resistance | K/W |
