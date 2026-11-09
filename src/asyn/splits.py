"""Study-grouped splitting and the fold-local preprocessing pipeline.

Two guardrails live here:
  * split by study, never by row (LeaveOneGroupOut on the study key);
  * every transform (impute, scale, encode) is fit inside the training fold
    only, by making it part of the sklearn Pipeline rather than a global step.
"""
from __future__ import annotations

import numpy as np
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.model_selection import LeaveOneGroupOut
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

CATEGORICAL = ["agitation", "buffer"]  # raw string condition columns


def study_cv() -> LeaveOneGroupOut:
    return LeaveOneGroupOut()


def assert_grouped(train_groups, test_groups) -> None:
    """Fail loudly if any study straddles the train/test boundary."""
    overlap = set(np.unique(train_groups)) & set(np.unique(test_groups))
    if overlap:
        raise AssertionError(f"study leakage across the split: {sorted(overlap)}")


def preprocessor(feature_cols):
    """ColumnTransformer: impute+scale numerics, one-hot the categoricals.

    Fitted per-fold because it is wrapped in the model Pipeline.
    """
    cat = [c for c in CATEGORICAL if c in feature_cols]
    num = [c for c in feature_cols if c not in cat]
    return ColumnTransformer(
        transformers=[
            ("num", Pipeline([
                ("impute", SimpleImputer(strategy="median")),
                ("scale", StandardScaler()),
            ]), num),
            ("cat", Pipeline([
                ("impute", SimpleImputer(strategy="most_frequent")),
                ("onehot", OneHotEncoder(handle_unknown="ignore")),
            ]), cat),
        ],
        remainder="drop",
    )


def model_pipeline(feature_cols, estimator) -> Pipeline:
    return Pipeline([("pre", preprocessor(feature_cols)), ("est", estimator)])
