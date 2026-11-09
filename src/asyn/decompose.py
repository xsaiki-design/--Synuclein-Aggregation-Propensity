"""Phase 4 — the variance decomposition (the project's main result).

Two complementary views of "how much is sequence vs. conditions":

  1. Nested models, scored by leave-one-study-out CV: conditions-only,
     sequence-only, sequence+structure, everything. The CV score is the honest
     "variance explained out of sample" for each block set.

  2. A mixed-effects null model with study as a random intercept. Its ICC is
     the share of total variance that sits at the lab/study level before any
     sequence feature is used — the clean statement of how much the conditions
     each lab measured under drive the reported rate.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import cross_val_score

from .features import build_features
from .splits import model_pipeline, study_cv

BLOCK_SETS = {
    "conditions only": ["conditions"],
    "sequence only": ["sequence"],
    "sequence + structure": ["sequence", "structure"],
    "everything": ["sequence", "structure", "conditions"],
}


def nested_cv_r2(df: pd.DataFrame) -> pd.DataFrame:
    X, groups, y, blocks = build_features(df)
    rows = []
    for name, keys in BLOCK_SETS.items():
        cols = [c for k in keys for c in blocks[k]]
        pipe = model_pipeline(cols, LinearRegression())
        scores = cross_val_score(
            pipe, X[cols], y, groups=groups, cv=study_cv(),
            scoring="r2",
        )
        rows.append({
            "model": name,
            "n_features": len(cols),
            "cv_r2_mean": float(np.mean(scores)),
            "cv_r2_std": float(np.std(scores)),
        })
    return pd.DataFrame(rows)


def study_icc(df: pd.DataFrame) -> dict:
    """Share of variance in log10(t_half) absorbed by a study random intercept."""
    import statsmodels.formula.api as smf

    d = df.copy()
    d["y"] = np.log10(d["t_half_h"].astype(float))
    md = smf.mixedlm("y ~ 1", d, groups=d["study"])
    res = md.fit(reml=True, method="lbfgs")
    group_var = float(res.cov_re.iloc[0, 0])
    resid_var = float(res.scale)
    icc = group_var / (group_var + resid_var) if (group_var + resid_var) else float("nan")
    return {
        "study_variance": group_var,
        "residual_variance": resid_var,
        "icc_study": icc,  # fraction of variance at the study/condition level
        "n_studies": int(d["study"].nunique()),
        "n_obs": int(len(d)),
    }


def run(df: pd.DataFrame) -> dict:
    table = nested_cv_r2(df)
    icc = study_icc(df)
    return {"nested_cv_r2": table, "study_icc": icc}
