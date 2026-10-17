"""Polarisation datasets bundled with the package.

Every dataset is returned as a tidy :class:`pandas.DataFrame` with columns

``j`` [mA/cm^2], ``T`` [K], ``V`` [V, per cell], ``P`` [bar abs.], ``w`` [-]

ready to be passed to :func:`aec_modelling.fit_polarization`. Metadata
(reference, DOI, operating conditions) are attached in ``df.attrs`` and
available through :func:`info`.

>>> from aec_modelling import datasets
>>> datasets.list_datasets()           # doctest: +SKIP
>>> df = datasets.load("sakas")        # doctest: +SKIP
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from importlib import resources

import pandas as pd


@dataclass(frozen=True)
class DatasetInfo:
    name: str
    filename: str
    description: str
    reference: str
    doi: str | None
    P: float | None  # system pressure [bar abs.]; None if given per row
    w: float | None  # KOH mass fraction; None if unknown
    notes: str = ""


DATASETS: dict[str, DatasetInfo] = {
    d.name: d
    for d in (
        DatasetInfo(
            name="ulleberg",
            filename="ulleberg.csv",
            description="PHOEBUS 26 kW pressurised alkaline electrolyser, 30-80 degC.",
            reference="Ulleberg O. Int J Hydrogen Energy 2003;28:21-33.",
            doi="10.1016/S0360-3199(02)00033-2",
            P=7.0,
            w=0.30,
        ),
        DatasetInfo(
            name="sakas",
            filename="sakas.csv",
            description="Industrial 3 MW alkaline electrolyser plant, 59-70 degC.",
            reference="Sakas G. et al. Int J Hydrogen Energy 2022;47:4328-45.",
            doi="10.1016/J.IJHYDENE.2021.11.126",
            P=16.0,
            w=0.25,
        ),
        DatasetInfo(
            name="sanchez",
            filename="sanchez.csv",
            description="15 kW alkaline stack (12 cells), 55-75 degC; per-cell voltage.",
            reference="Sanchez M. et al. Int J Hydrogen Energy 2018;43:20332-45.",
            doi="10.1016/j.ijhydene.2018.09.029",
            P=7.0,
            w=0.35,
        ),
        DatasetInfo(
            name="degroot_pressure",
            filename="degroot_pressure.csv",
            description="Advanced AWE single cell at 50 degC, 10 and 30 bar.",
            reference="de Groot MT, Kraakman J, Garcia Barros RL. Int J Hydrogen Energy 2022.",
            doi="10.1016/j.ijhydene.2022.08.075",
            P=None,
            w=0.30,
        ),
        DatasetInfo(
            name="degroot_temperature",
            filename="degroot_temperature.csv",
            description="Advanced AWE single cell at 30 bar, 50-70 degC.",
            reference="de Groot MT, Kraakman J, Garcia Barros RL. Int J Hydrogen Energy 2022.",
            doi="10.1016/j.ijhydene.2022.08.075",
            P=30.0,
            w=0.30,
        ),
        DatasetInfo(
            name="univpm",
            filename="univpm.csv",
            description="Field data of a commercial AWE stack logged at UNIVPM (June 2023), "
            "1-minute resolution; pressure and stack temperature per row.",
            reference="Jin L. et al. Appl Therm Eng 2025;262:125154 (new experimental data).",
            doi="10.1016/j.applthermaleng.2024.125154",
            P=None,
            w=None,
            notes="KOH mass fraction was not logged: pass w=... when fitting. "
            "Current density derived from the measured H2 flow via Faraday's law.",
        ),
    )
}

_COMMON_NOTE = (
    "Pressure / KOH fraction not stated in the raw file are the values used in the "
    "original analysis scripts; verify against the source before quantitative use."
)


def list_datasets() -> pd.DataFrame:
    """Overview table of the bundled datasets."""
    rows = []
    for d in DATASETS.values():
        df = load(d.name)
        rows.append(
            {
                "name": d.name,
                "points": len(df),
                "T range [K]": f"{df['T'].min():.1f}-{df['T'].max():.1f}",
                "j max [mA/cm2]": round(df["j"].max(), 1),
                "P [bar]": d.P if d.P is not None else "per row",
                "w [-]": d.w if d.w is not None else "unknown",
                "description": d.description,
            }
        )
    return pd.DataFrame(rows).set_index("name")


def info(name: str) -> DatasetInfo:
    try:
        return DATASETS[name]
    except KeyError:
        raise KeyError(f"Unknown dataset {name!r}. Available: {list(DATASETS)}") from None


def load(name: str) -> pd.DataFrame:
    """Load a bundled dataset as a tidy DataFrame (metadata in ``df.attrs``)."""
    meta = info(name)
    with resources.files(__package__).joinpath("data", meta.filename).open("r") as fh:
        df = pd.read_csv(fh)
    df.attrs.update(asdict(meta))
    df.attrs["notes"] = (meta.notes + " " + _COMMON_NOTE).strip()
    return df
