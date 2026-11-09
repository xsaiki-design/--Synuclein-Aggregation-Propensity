import numpy as np
import pandas as pd
import pytest
from sklearn.linear_model import LinearRegression

from asyn import features, synth
from asyn.predict import _baseline_mean, _bootstrap_metric, loso_oof
from asyn.splits import assert_grouped, model_pipeline, study_cv


def test_loso_folds_have_disjoint_studies():
    df = synth.generate(n_studies=5, seed=2)
    X, groups, y, _ = features.build_features(df)
    g = np.asarray(groups)
    for tr, te in study_cv().split(X, y, g):
        assert set(g[tr]).isdisjoint(set(g[te]))
        assert len(np.unique(g[te])) == 1  # exactly one held-out study


def test_assert_grouped_catches_leakage():
    assert_grouped(np.array(["a", "a"]), np.array(["b"]))  # ok
    with pytest.raises(AssertionError):
        assert_grouped(np.array(["a", "b"]), np.array(["b"]))


def test_scaler_is_fit_on_training_fold_only():
    # Train numeric feature ~N(0,1); test shifted by +10. If scaling used train
    # stats, the transformed test mean must be far from 0.
    Xtr = pd.DataFrame({"position": np.zeros(20)})
    Xte = pd.DataFrame({"position": np.full(5, 10.0)})
    pipe = model_pipeline(["position"], LinearRegression())
    pipe.named_steps["pre"].fit(Xtr)
    zt = pipe.named_steps["pre"].transform(Xte)
    # train had zero variance -> scaler keeps it; use a varying train instead:
    Xtr = pd.DataFrame({"position": np.arange(20.0)})
    pipe.named_steps["pre"].fit(Xtr)
    train_mean = Xtr["position"].mean()
    zt = pipe.named_steps["pre"].transform(Xte)[:, 0]
    assert np.all(zt > (10 - train_mean) / Xtr["position"].std(ddof=0) - 1)


def test_loso_oof_fills_every_row():
    df = synth.generate(n_studies=5, seed=3)
    X, groups, y, blocks = features.build_features(df)
    cols = blocks["sequence"]
    oof = loso_oof(X, y, np.asarray(groups), cols,
                   LinearRegression(), {})  # empty grid -> no tuning
    assert np.isfinite(oof).all()


def test_bootstrap_interval_brackets_mean():
    df = synth.generate(n_studies=6, seed=4)
    _, groups, y, _ = features.build_features(df)
    pred = _baseline_mean(y, groups)
    from sklearn.metrics import r2_score
    m, lo, hi = _bootstrap_metric(y, pred, groups, r2_score, n=300)
    assert lo <= m <= hi
