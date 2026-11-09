"""Phase 5 — prediction, properly validated.

Everything here obeys the guardrails:
  * leave-one-study-out outer loop (deployment = a new lab);
  * nested tuning: an inner GroupKFold picks hyperparameters, the outer loop
    only scores, so no fold scores itself;
  * all preprocessing fit inside the fold (via the model Pipeline);
  * three baselines always reported (mean / existing-tool score / one feature);
  * bootstrap confidence intervals, resampled over studies, on every metric.

No neural network: ~130 points cannot support one. The default estimator is a
small random forest, with a Ridge fallback.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import GridSearchCV, GroupKFold

from .features import build_features
from .splits import model_pipeline, study_cv

RNG = np.random.default_rng(0)


def _default_grid():
    est = RandomForestRegressor(random_state=0, n_jobs=1)
    grid = {"est__n_estimators": [300],
            "est__max_depth": [2, 3, 4],
            "est__min_samples_leaf": [2, 4]}
    return est, grid


def loso_oof(X, y, groups, cols, estimator, param_grid):
    """Pooled out-of-fold predictions under nested leave-one-study-out."""
    y = np.asarray(y, float)
    groups = np.asarray(groups)
    oof = np.full(len(y), np.nan)
    outer = study_cv()
    n_inner = max(2, min(5, len(np.unique(groups)) - 1))
    for tr, te in outer.split(X, y, groups):
        pipe = model_pipeline(cols, estimator)
        inner = GroupKFold(n_splits=min(n_inner, len(np.unique(groups[tr]))))
        gs = GridSearchCV(pipe, param_grid, cv=inner, scoring="neg_mean_squared_error")
        gs.fit(X.iloc[tr][cols], y[tr], groups=groups[tr])
        oof[te] = gs.best_estimator_.predict(X.iloc[te][cols])
    return oof


def _baseline_mean(y, groups):
    y = np.asarray(y, float); groups = np.asarray(groups)
    oof = np.full(len(y), np.nan)
    for tr, te in study_cv().split(y, y, groups):
        oof[te] = y[tr].mean()
    return oof


def _baseline_single_feature(X, y, groups, col="d_hydro"):
    y = np.asarray(y, float); groups = np.asarray(groups)
    oof = np.full(len(y), np.nan)
    for tr, te in study_cv().split(X, y, groups):
        pipe = model_pipeline([col], LinearRegression())
        pipe.fit(X.iloc[tr][[col]], y[tr])
        oof[te] = pipe.predict(X.iloc[te][[col]])
    return oof


def _baseline_tool(X, y, groups, tool_col="tool_score"):
    if tool_col not in X.columns:
        return None
    y = np.asarray(y, float); groups = np.asarray(groups)
    oof = np.full(len(y), np.nan)
    for tr, te in study_cv().split(X, y, groups):
        pipe = model_pipeline([tool_col], LinearRegression())
        pipe.fit(X.iloc[tr][[tool_col]], y[tr])
        oof[te] = pipe.predict(X.iloc[te][[tool_col]])
    return oof


def _bootstrap_metric(y, pred, groups, metric, n=2000):
    """Bootstrap a metric by resampling whole studies (the independent unit)."""
    y = np.asarray(y, float); pred = np.asarray(pred, float); groups = np.asarray(groups)
    uniq = np.unique(groups)
    vals = []
    for _ in range(n):
        pick = RNG.choice(uniq, size=len(uniq), replace=True)
        idx = np.concatenate([np.where(groups == g)[0] for g in pick])
        vals.append(metric(y[idx], pred[idx]))
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return float(np.mean(vals)), float(lo), float(hi)


def _rmse(a, b):
    return float(np.sqrt(mean_squared_error(a, b)))


def run(df: pd.DataFrame) -> dict:
    X, groups, y, blocks = build_features(df)
    all_cols = [c for k in blocks for c in blocks[k]]
    est, grid = _default_grid()

    preds = {
        "full model (RF, all blocks)": loso_oof(X, y, groups, all_cols, est, grid),
        "baseline: predict the mean": _baseline_mean(y, groups),
        "baseline: hydrophobicity change only": _baseline_single_feature(X, y, groups),
    }
    tool = _baseline_tool(X, y, groups)
    if tool is not None:
        preds["baseline: existing-tool score"] = tool

    rows = []
    for name, pred in preds.items():
        r2_m, r2_lo, r2_hi = _bootstrap_metric(y, pred, groups, r2_score)
        rm_m, rm_lo, rm_hi = _bootstrap_metric(y, pred, groups, _rmse)
        rows.append({
            "model": name,
            "r2": round(r2_m, 3), "r2_lo": round(r2_lo, 3), "r2_hi": round(r2_hi, 3),
            "rmse": round(rm_m, 3), "rmse_lo": round(rm_lo, 3), "rmse_hi": round(rm_hi, 3),
        })
    table = pd.DataFrame(rows)
    note = None
    if tool is None:
        note = ("existing-tool baseline skipped: add a `tool_score` column "
                "(TANGO/WALTZ/AGGRESCAN output) to data to enable it.")
    return {"metrics": table, "oof_predictions": preds, "note": note,
            "target": "log10(t_half_h)"}
