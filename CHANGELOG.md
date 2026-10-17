# Changelog

## 0.1.0 (unreleased)

The research scripts were restructured into the installable package `aec-modelling`.

### Added
- `AECModel`: a vectorised four-parameter cell model. It provides a voltage breakdown,
  the inverse (V → j), efficiency, heat generation and H₂ production. There is an
  optional `asinh` activation form.
- `fit_polarization` / `FitResult`: fitting on long-format data, with per-point P and w,
  standard errors, R²/RMSE/MAPE, evaluation per temperature, and `sample_points` for
  calibration on random subsets.
- `Stack`, `LumpedThermalModel` (0D transient / steady state / lye sizing) and
  `AxialThermalModel` (1D temperature and current distribution).
- `datasets`: tidy literature datasets (Ulleberg, Sakas, Sánchez, de Groot) and the
  UNIVPM field data, with metadata.
- Optional plotting helpers, tests, examples, CI, and a model description in `docs/model.md`.

### Changed (compared with the original `src/func.py` and `Main.py`)
- CODATA constants: R = 8.314462618 instead of 8.31, F = 96485.33 instead of 96500.
- °C conversion uses 273.15 instead of 273.
- KOH molality for the vapour pressure and water activity is now per kg of **water**:
  `1000 w / (M_KOH (1 - w))`. The old code divided by kg of solution. The reversible
  voltage changes by a few mV. Refit the parameters rather than reusing old values.
- Shomate enthalpies are no longer rounded.
- The original scripts are still available in git history (commit `e086dec`).
