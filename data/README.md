# Data

| Folder | Content |
|--------|---------|
| `raw/` | Original digitised and measured files as collected during the project (CSV/XLSX). |
| `legacy_fitting_parameters/` | Kibria / Miles exchange-current and Tafel-limit fits from the earlier, detailed anode/cathode model. Each file is available as both `.pkl` and `.csv`. The four-parameter model does not use them. |

The package ships tidy copies of the polarisation datasets in
`src/aec_modelling/data/`, with columns `j, T, V, P, w`. To regenerate them
after editing `raw/` or the metadata in `aec_modelling/datasets.py`, run:

```bash
python scripts/prepare_datasets.py
```

### Raw files not (yet) packaged

| File | Content |
|------|---------|
| `COMSOL.csv`, `Experiment.csv` | COMSOL model and experimental curve from a bachelor thesis (single temperature) |
| `roy.csv`, `ursua.csv` | Single-temperature literature polarisation curves |
| `Vermeiren 30_wt.csv` | Zirfon separator area resistance vs. temperature (Vermeiren et al. 1998) |
| `Temperature variation overvoltage.csv` | Overvoltage vs. temperature |
| `Sakas *.csv` | Individual Sakas curves (already merged in `SAKAS.csv`) |
