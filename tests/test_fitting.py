import numpy as np
import pytest

import aec_modelling as am
from aec_modelling import datasets

LITERATURE = ["ulleberg", "sakas", "sanchez", "degroot_pressure", "degroot_temperature"]


def test_fit_recovers_synthetic_parameters():
    true = am.AECModel(delta_l=0.5, alpha=0.15, C1=-5500, C2=16, P=10, w=0.3)
    j = np.tile(np.linspace(50, 600, 15), 3)
    T = np.repeat([323.15, 343.15, 363.15], 15)
    V = true.voltage(j, T)
    fit = am.fit_polarization(j=j, T=T, V=V, P=10, w=0.3)
    for k, v in true.params.items():
        assert fit.params[k] == pytest.approx(v, rel=1e-3)
    assert fit.r2 > 0.99999


@pytest.mark.parametrize("name", LITERATURE)
def test_literature_datasets_fit_well(name):
    df = datasets.load(name)
    assert set("jTVPw") <= set(df.columns)
    fit = am.fit_polarization(df)
    assert fit.r2 > 0.98
    assert fit.rmse < 0.015


def test_extrapolation_from_small_sample():
    df = datasets.load("ulleberg")
    calib = am.sample_points(df, n=15, random_state=0)
    assert calib["T"].nunique() >= 2
    fit = am.fit_polarization(calib)
    assert fit.evaluate(df)["r2"].iloc[0] > 0.95
    per_T = fit.evaluate(df, by="T")
    assert len(per_T) == df["T"].nunique()


def test_missing_operating_conditions_raise():
    df = datasets.load("univpm")
    with pytest.raises(ValueError, match="NaN"):
        am.fit_polarization(df)


def test_dataset_metadata():
    df = datasets.load("sakas")
    assert df.attrs["doi"].startswith("10.")
    assert len(datasets.list_datasets()) == len(datasets.DATASETS)
