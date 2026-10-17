# aec-modelling: four-parameter thermo-electrochemical model of alkaline electrolysers

[![tests](https://github.com/Lingkangjin/AEC-Modelling/actions/workflows/tests.yml/badge.svg)](https://github.com/Lingkangjin/AEC-Modelling/actions/workflows/tests.yml)
[![PyPI](https://img.shields.io/pypi/v/aec-modelling.svg)](https://pypi.org/project/aec-modelling/)
[![Python](https://img.shields.io/pypi/pyversions/aec-modelling.svg)](https://pypi.org/project/aec-modelling/)
[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.applthermaleng.2024.125154-blue)](https://doi.org/10.1016/j.applthermaleng.2024.125154)

A lightweight Python package for **alkaline water electrolysis (AWE / AEC)**. It provides:

- an **electrochemical model** that describes the polarisation curve with only **four fitted parameters**;
- **fitting tools** to calibrate those parameters on your own measurements, from as few as ~20 points;
- **thermal models** of the stack: a 0D transient energy balance and a 1D along-the-channel model;
- **literature datasets** ready to fit (Ulleberg, Sakas, Sánchez, de Groot) and field data from UNIVPM.

It is intended for **system-level studies** such as power-to-hydrogen plants, temperature
control, lye-cooling sizing and techno-economic analysis. In these studies you need a model
that is fast, robust and extrapolates sensibly in temperature.

The package was developed during a visiting PhD period (2022) at DTU Energy under the
supervision of Henrik Lund Frandsen.

<p float="left">
  <img src="https://raw.githubusercontent.com/Lingkangjin/AEC-Modelling/main/docs/images/logounivpm.png" alt="UNIVPM" width="150"/>
  <img src="https://raw.githubusercontent.com/Lingkangjin/AEC-Modelling/main/docs/images/logodtu.png" alt="DTU" width="100"/>
</p>

---

## Contents

- [Installation](#installation)
- [Quick start](#quick-start)
- [The model in one minute](#the-model-in-one-minute)
- [Fitting your own data](#fitting-your-own-data)
- [Thermal modelling](#thermal-modelling)
- [Bundled datasets](#bundled-datasets)
- [Examples](#examples)
- [Repository structure](#repository-structure)
- [Development](#development)
- [Citation](#citation)
- [References](#references)

---

## Installation

```bash
pip install aec-modelling            # core: numpy, scipy, pandas
pip install "aec-modelling[plot]"    # + matplotlib plotting helpers
```

From source, for development:

```bash
git clone https://github.com/Lingkangjin/AEC-Modelling.git
cd AEC-Modelling
uv sync --extra dev          # or: python -m pip install -e ".[dev]"
```

Python ≥ 3.9 is required.

## Quick start

```python
import numpy as np
from aec_modelling import AECModel, datasets, fit_polarization

# 1. Fit the four parameters on a literature dataset (or your own DataFrame)
df = datasets.load("sakas")          # columns: j [mA/cm2], T [K], V [V], P [bar], w [-]
fit = fit_polarization(df)
print(fit.summary())

# 2. Use the calibrated model
cell = fit.model
cell.voltage(300, 353.15)                    # V at 300 mA/cm2, 80 °C
cell.overpotentials(np.linspace(10, 500, 50), 353.15)   # DataFrame: V_rev, V_act, V_ohm, V_cell, V_tn
cell.current_density(1.9, 353.15)            # inverse: j at 1.9 V
cell.voltage_efficiency(300, 353.15)         # V_tn / V_cell
```

```text
Four-parameter AEC model fit
--------------------------------------------
 delta_l =      0.59344 +/- 0.0364     [cm]
   alpha =      0.13121 +/- 0.00803    [-]
      C1 =      -6646.3 +/- 397        [K]
      C2 =       20.816 +/- 1          [-]
--------------------------------------------
points = 40   R2 = 0.9982   RMSE = 3.6 mV
activation energy Ea = 55.3 kJ/mol
```

## The model in one minute

```text
V_cell = V_rev + V_act + V_ohm

V_rev = U0(T) + RT/(2F) · ln[ (P − p_v)^1.5 / a_w ]     reversible
V_act = RT/(2·alpha·F) · ( ln j − C1/T − C2 )            activation (Tafel)
V_ohm = delta_l · j / sigma(T, w)                        ohmic
```

| Parameter | Unit | Physical meaning |
|-----------|------|------------------|
| `delta_l` | cm   | effective ionic path length (gap + separator + bubbles) |
| `alpha`   | –    | lumped charge-transfer coefficient, which sets the Tafel slope |
| `C1`      | K    | Arrhenius slope of the exchange current density (activation energy Ea = −C1·R) |
| `C2`      | –    | Arrhenius intercept of ln i0 |

The reversible voltage, electrolyte conductivity, vapour pressure and water activity come
from literature correlations that depend on temperature, pressure and KOH concentration.
Only the four parameters above are fitted. The full equations are in
**[docs/model.md](https://github.com/Lingkangjin/AEC-Modelling/blob/main/docs/model.md)**.

<p float="left">
  <img src="https://raw.githubusercontent.com/Lingkangjin/AEC-Modelling/main/docs/images/01_quickstart.png" alt="Polarisation curves and voltage breakdown" width="700"/>
</p>

## Fitting your own data

Put your measurements in a DataFrame with one row per point:

| column | unit | |
|--------|------|-|
| `j` | mA/cm² | current density (> 0) |
| `T` | K | cell / stack temperature |
| `V` | V | **per-cell** voltage |
| `P` | bar abs. | optional column, or pass `P=...` |
| `w` | – | KOH mass fraction, optional column, or pass `w=...` |

```python
import pandas as pd
from aec_modelling import fit_polarization, sample_points

df = pd.read_csv("my_test_bench.csv")          # j, T, V
fit = fit_polarization(df, P=10.0, w=0.30)

# Calibrate on 20 random points spanning 2+ temperatures and validate on the rest
calib = sample_points(df, n=20, min_temperatures=2, random_state=0)
fit = fit_polarization(calib, P=10.0, w=0.30)
print(fit.evaluate(df, by="T"))                # R2 / RMSE / MAPE per temperature

# Custom initial guess and bounds: (delta_l, alpha, C1, C2)
fit = fit_polarization(df, P=10, w=0.3, p0=(0.4, 0.13, -6330, 14),
                       bounds=((0, 0.05, -1e4, -5), (2, 1, 0, 30)))
```

**Tips**

- Use data at **two or more temperatures**. With a single temperature, `C1` and `C2` cannot be separated.
- For stack data, divide the stack voltage by the number of cells.
- `activation="asinh"` avoids negative activation overpotential when the data include very low current densities.

<p float="left">
  <img src="https://raw.githubusercontent.com/Lingkangjin/AEC-Modelling/main/docs/images/02_fit_sakas.png" alt="Calibration on 20 points, validation on all" width="700"/>
</p>

## Thermal modelling

A cell releases heat when its voltage is above the thermoneutral voltage V_tn ≈ 1.48 V:
Q_gen = N_cells · (V_cell − V_tn) · I. The package couples this heat source to two stack
energy-balance models.

```python
from aec_modelling import Stack, LumpedThermalModel, AxialThermalModel

stack = Stack(fit.model, n_cells=100, cell_area=2.6e4)     # area in cm2
I = stack.current(300)                                      # A at 300 mA/cm2

# 0D: well-mixed stack with lye cooling and heat losses
th = LumpedThermalModel(stack, heat_capacity=2e7, thermal_resistance=0.01)
m_lye = th.required_lye_flow(I, T_out=353.15, T_in=343.15)   # kg/s
T_ss  = th.steady_state(I, T_in=343.15, m_lye=m_lye)
t     = np.linspace(0, 4 * 3600, 400)                       # s
sim   = th.simulate(t, I=lambda t: I if t < 7200 else I/2, T_in=343.15, m_lye=m_lye, T0=313.15)

# 1D: temperature rise along the cell and the resulting current maldistribution
res = AxialThermalModel(stack, n_segments=50).solve(I, T_in=343.15, m_lye=8.0)
res.V_cell, res.T_out, res.current_maldistribution, res.profile    # profile: z, T, j, q
```

| 0D: warm-up and load step | 1D: profiles vs. lye flow |
|---|---|
| <img src="https://raw.githubusercontent.com/Lingkangjin/AEC-Modelling/main/docs/images/04_lumped_thermal.png" width="380"/> | <img src="https://raw.githubusercontent.com/Lingkangjin/AEC-Modelling/main/docs/images/05_axial_thermal.png" width="420"/> |

## Bundled datasets

```python
from aec_modelling import datasets
datasets.list_datasets()        # overview table
df = datasets.load("ulleberg")  # tidy DataFrame; metadata in df.attrs (reference, doi, notes)
```

| name | source | T range | P | KOH | points | R² (all points) |
|------|--------|---------|---|-----|--------|-----------------|
| `ulleberg` | Ulleberg 2003, PHOEBUS 26 kW | 30–80 °C | 7 bar | 30 wt% | 86 | 0.997 |
| `sakas` | Sakas et al. 2022, industrial 3 MW | 59–70 °C | 16 bar | 25 wt% | 40 | 0.998 |
| `sanchez` | Sánchez et al. 2018, 15 kW stack | 55–75 °C | 7 bar | 35 wt% | 24 | 0.998 |
| `degroot_pressure` | de Groot et al. 2022 | 50 °C | 10/30 bar | 30 wt% | 13 | 0.999 |
| `degroot_temperature` | de Groot et al. 2022 | 50–70 °C | 30 bar | 30 wt% | 24 | 0.998 |
| `univpm` | UNIVPM field data (2023), 1-min log | 32–56 °C | 1–5.5 bar | not logged | 184 | n/a (raw) |

> **Note:** pressure and KOH fraction that are not part of the raw files are the values used
> in the original analysis scripts. Check them against the source paper before using the
> fitted parameters quantitatively. The `univpm` set is raw field data that includes
> transients. It is provided for exploration, and fitting it requires `w=...`.

<p float="left">
  <img src="https://raw.githubusercontent.com/Lingkangjin/AEC-Modelling/main/docs/images/03_fit_all_datasets.png" alt="Fits on all literature datasets" width="800"/>
</p>

## Examples

| Script | What it shows |
|--------|---------------|
| [`examples/01_quickstart.py`](https://github.com/Lingkangjin/AEC-Modelling/blob/main/examples/01_quickstart.py) | polarisation curves at several temperatures, voltage breakdown, inverse model |
| [`examples/02_fit_sakas.py`](https://github.com/Lingkangjin/AEC-Modelling/blob/main/examples/02_fit_sakas.py) | calibration on 20 random points and validation on the full dataset (the original `Main.py` workflow) |
| [`examples/03_fit_all_datasets.py`](https://github.com/Lingkangjin/AEC-Modelling/blob/main/examples/03_fit_all_datasets.py) | fits every literature dataset and compares parameters and activation energies |
| [`examples/04_lumped_thermal.py`](https://github.com/Lingkangjin/AEC-Modelling/blob/main/examples/04_lumped_thermal.py) | 0D stack: lye-flow sizing, cold start, load step, heat balance |
| [`examples/05_axial_thermal.py`](https://github.com/Lingkangjin/AEC-Modelling/blob/main/examples/05_axial_thermal.py) | 1D cell: temperature and current profiles for different lye flows, 1D vs 0D |

```bash
python examples/02_fit_sakas.py     # figures are saved to examples/figures/
```

## Repository structure

```
AEC-Modelling/
├── src/aec_modelling/        # the installable package
│   ├── electrochemistry.py   # AECModel, reversible voltage, cell_voltage()
│   ├── fitting.py            # fit_polarization, FitResult, sample_points, metrics
│   ├── thermal.py            # Stack, LumpedThermalModel (0D), AxialThermalModel (1D)
│   ├── electrolyte.py        # KOH density, molarity, conductivity, vapour pressure, ...
│   ├── thermo.py             # Shomate data, ΔH/ΔS/ΔG, thermoneutral voltage
│   ├── datasets.py           # dataset registry and loaders
│   ├── plotting.py           # optional matplotlib helpers
│   └── data/*.csv            # tidy datasets shipped with the package
├── examples/                 # runnable scripts
├── tests/                    # pytest suite
├── data/raw/                 # original raw data files
├── data/legacy_fitting_parameters/  # Kibria/Miles fits of the earlier detailed model
├── scripts/prepare_datasets.py      # raw -> tidy dataset conversion
└── docs/                     # model description and figures
```

## Development

```bash
uv sync --extra dev
uv run pytest                 # tests
uv run ruff check .           # lint
uv run ruff format .          # format
uv build                      # sdist + wheel in dist/
```

**Adding a dataset:** put the raw file in `data/raw/`, add a `DatasetInfo` entry in
`src/aec_modelling/datasets.py` and a builder function in `scripts/prepare_datasets.py`,
then run `python scripts/prepare_datasets.py`.

**Releasing to PyPI:** bump `__version__` in `src/aec_modelling/__init__.py`, update
`CHANGELOG.md`, and create a GitHub release. The `publish` workflow then uploads the release
through PyPI trusted publishing. Configure the trusted publisher on pypi.org once before the
first release.

Contributions are welcome through issues and pull requests.

## Citation

If you use this package, please cite:

[1] Lingkang Jin, Rafael Nogueira Nakashima, Gabriele Comodi, Henrik Lund Frandsen. *Alkaline electrolysis for green hydrogen production: A novel, simple model for thermo-electrochemical coupled system analysis.* Applied Thermal Engineering 262 (2025) 125154. https://doi.org/10.1016/j.applthermaleng.2024.125154

[2] Lingkang Jin, Rafael Nogueira Nakashima, Henrik Lund Frandsen, Gabriele Comodi. *Alkaline Electrolysis for Green Hydrogen Production: Techno-Economic Analysis of Temperature Influence and Control.* ECOS 2023, 908–919. https://doi.org/10.52202/069564-0082

```bibtex
@article{JIN2025125154,
  title   = {Alkaline electrolysis for green hydrogen production: A novel, simple model for thermo-electrochemical coupled system analysis},
  journal = {Applied Thermal Engineering},
  volume  = {262},
  pages   = {125154},
  year    = {2025},
  issn    = {1359-4311},
  doi     = {10.1016/j.applthermaleng.2024.125154},
  author  = {Lingkang Jin and Rafael {Nogueira Nakashima} and Gabriele Comodi and Henrik Lund Frandsen}
}

@inproceedings{Jin2023AlkalineControl,
  title     = {{Alkaline Electrolysis for Green Hydrogen Production: Techno-Economic Analysis of Temperature Influence and Control}},
  booktitle = {36th International Conference on Efficiency, Cost, Optimization, Simulation and Environmental Impact of Energy Systems (ECOS 2023)},
  author    = {Jin, Lingkang and Nogueira Nakashima, Rafael and Lund Frandsen, Henrik and Comodi, Gabriele},
  year      = {2023},
  pages     = {908--919},
  address   = {Las Palmas de Gran Canaria, Spain},
  isbn      = {978-1-7138-7492-8},
  doi       = {10.52202/069564-0082}
}
```

GitHub also reads [`CITATION.cff`](https://github.com/Lingkangjin/AEC-Modelling/blob/main/CITATION.cff) and shows a "Cite this repository" button.

## License

Released under the [MIT License](https://github.com/Lingkangjin/AEC-Modelling/blob/main/LICENSE).

## Code contributors

- Lingkang Jin (lingkang32@gmail.com / l.jin@tue.nl)
- Rafael Nogueira Nakashima (rafnn@dtu.dk)

## References

### Datasets
1. Ulleberg Ø. Modeling of advanced alkaline electrolyzers: a system simulation approach. Int J Hydrogen Energy 2003;28:21–33. https://doi.org/10.1016/S0360-3199(02)00033-2
2. Sakas G, Ibáñez-Rioja A, Ruuskanen V, Kosonen A, Ahola J, Bergmann O. Dynamic energy and mass balance model for an industrial alkaline water electrolyzer plant process. Int J Hydrogen Energy 2022;47:4328–45. https://doi.org/10.1016/J.IJHYDENE.2021.11.126
3. Sánchez M, Amores E, Rodríguez L, Clemente-Jul C. Semi-empirical model and experimental validation for the performance evaluation of a 15 kW alkaline water electrolyzer. Int J Hydrogen Energy 2018;43:20332–45. https://doi.org/10.1016/j.ijhydene.2018.09.029
4. de Groot MT, Kraakman J, Garcia Barros RL. Optimal operating parameters for advanced alkaline water electrolysis. Int J Hydrogen Energy 2022. https://doi.org/10.1016/j.ijhydene.2022.08.075

### Reversible voltage
- Domain E, Conductivity SE. A Comprehensive Survey of Alkaline Electrolyzer Modeling: Electrical Domain and Specific Electrolyte Conductivity 2022.
- Kibria MF, Mridha MS. Electrochemical studies of the nickel electrode for the oxygen evolution reaction. 1996;21.

### Electrode kinetics
- Miles MH, Kissel G, Lu PWT, Srinivasan S. Effect of temperature on electrode kinetic parameters for hydrogen and oxygen evolution reactions on nickel electrodes in alkaline solutions. J Electrochem Soc 1976;123:332–6. https://doi.org/10.1149/1.2132820
- Kibria MF, Mridha MS, Khan AH. Electrochemical studies of a nickel electrode for the hydrogen evolution reaction. Int J Hydrogen Energy 1995;20:435–40. https://doi.org/10.1016/0360-3199(94)00073-9
- Huot J-Y. Hydrogen evolution and interface phenomena on a nickel cathode in 30 w/o KOH. J Electrochem Soc 1989;136:1933–9. https://doi.org/10.1149/1.2097088

### Electrolyte conductivity, density and concentration
- Gilliam RJ, Graydon JW, Kirk DW, Thorpe SJ. A review of specific conductivities of potassium hydroxide solutions for various concentrations and temperatures. Int J Hydrogen Energy 2007;32:359–64. https://doi.org/10.1016/j.ijhydene.2006.10.062

### Separator, gas solubility and diffusion
- Vermeiren P, Adriansens W, Moreels JP, Leysen R. Evaluation of the Zirfon® separator for use in alkaline water electrolysis and Ni-H2 batteries. Int J Hydrogen Energy 1998;23:321–4. https://doi.org/10.1016/s0360-3199(97)00069-4
- Tromans D. Temperature and pressure dependent solubility of oxygen in water: a thermodynamic analysis. Hydrometallurgy 1998;48:327–42. https://doi.org/10.1016/s0304-386x(98)00007-3
- Tham MK, Walker RD, Gubbins KE. Diffusion of oxygen and hydrogen in aqueous potassium hydroxide solutions. J Phys Chem 1970;74:1747–51. https://doi.org/10.1021/j100703a015

### Bubbles
- Vogt H, Balzer RJ. The bubble coverage of gas-evolving electrodes in stagnant electrolytes. Electrochim Acta 2005;50:2073–9. https://doi.org/10.1016/j.electacta.2004.09.025
- Hammoudi M, Henao C, Agbossou K, Dubé Y, Doumbia ML. New multi-physics approach for modelling and design of alkaline electrolyzers. Int J Hydrogen Energy 2012;37:13895–913. https://doi.org/10.1016/j.ijhydene.2012.07.015
- Mandin P, Derhoumi Z, Roustan H, Rolf W. Bubble over-potential during two-phase alkaline water electrolysis. Electrochim Acta 2014;128:248–58. https://doi.org/10.1016/j.electacta.2013.11.068
- Henao C, Agbossou K, Hammoudi M, Dubé Y, Cardenas A. Simulation tool based on a physics model and an electrical analogy for an alkaline electrolyser. J Power Sources 2014;250:58–67. https://doi.org/10.1016/j.jpowsour.2013.10.086
