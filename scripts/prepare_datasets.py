"""Convert the raw files in ``data/raw`` into the tidy CSVs shipped with the package.

Every output file has the columns ``j`` [mA/cm^2], ``T`` [K], ``V`` [V],
``P`` [bar abs.] and ``w`` [-]. Operating conditions that are not part of the
raw files (pressure, KOH fraction) come from ``aec_modelling.datasets.DATASETS``,
so that the metadata lives in a single place.

Run from the repository root::

    python scripts/prepare_datasets.py
"""

from pathlib import Path

import numpy as np
import pandas as pd

from aec_modelling.constants import KELVIN
from aec_modelling.datasets import DATASETS

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
OUT = ROOT / "src" / "aec_modelling" / "data"


def read(name):
    return pd.read_csv(RAW / name, encoding="latin1")


def ulleberg():
    d = read("Ullberg exp.csv")
    return pd.DataFrame({"j": d.iloc[:, 0], "T": d.iloc[:, 2] + KELVIN, "V": d.iloc[:, 1]})


def sakas():
    d = read("SAKAS.csv")
    return pd.DataFrame({"j": d["J [mA/cm^2]"], "T": d["T [K]"], "V": d["V [V]"]})


def sanchez():
    d = read("Sanchez_data.csv")
    return pd.DataFrame({"j": d["I (A/cm^2)"] * 1000, "T": d["T"] + KELVIN, "V": d["V_cell"]})


def degroot_pressure():
    d = read("De Groot.csv")
    return pd.DataFrame(
        {"j": d.iloc[:, 0], "T": d.iloc[:, 2] + KELVIN, "V": d.iloc[:, 1], "P": d.iloc[:, 3]}
    )


def degroot_temperature():
    d = read("De_Groot_2.csv")
    return pd.DataFrame(
        {"j": d.iloc[:, 0], "T": d.iloc[:, 2] + KELVIN, "V": d.iloc[:, 1], "P": d.iloc[:, 3]}
    )


def univpm():
    d = read("UNIVPM_data.csv")
    d = d[d["J from faraday $[mA/cm^2]$"] > 0]
    return pd.DataFrame(
        {
            "j": d["J from faraday $[mA/cm^2]$"],
            "T": d["AIM_Stack_1_temp"] + KELVIN,
            "V": d["Cell voltage [V]"],
            "P": d["Pres. cat [bar]"] + 1.01325,  # logged as gauge pressure
            "time": d["Timeindex"],
        }
    )


BUILDERS = {
    f.__name__: f for f in (ulleberg, sakas, sanchez, degroot_pressure, degroot_temperature, univpm)
}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    for name, build in BUILDERS.items():
        info = DATASETS[name]
        df = build()
        if "P" not in df:
            df["P"] = info.P
        df["w"] = np.nan if info.w is None else info.w
        cols = ["j", "T", "V", "P", "w"] + [c for c in df if c not in ("j", "T", "V", "P", "w")]
        df = df[cols].sort_values(["T", "j"]).reset_index(drop=True)
        df.to_csv(OUT / info.filename, index=False, float_format="%.6g")
        print(f"{name:<20} {len(df):>4} rows -> {info.filename}")


if __name__ == "__main__":
    main()
