"""The real sequence-vs-conditions analysis, on the data this SI supports.

For the 6 familial mutants, each condition table gives an aggregation score per
(mutant × condition category) with a support count. The question "how much is
sequence vs. conditions" becomes, at this level:

  how much of the variation in these scores is explained by WHICH MUTANT
  (sequence identity) versus WHICH CONDITION CATEGORY (buffer/salt/shaking/
  concentration)?

For each condition dimension we fit support-weighted linear models and compare
variance explained (R²): mutant identity only, condition category only, and
both. We also report, in the paper's own framing, the between-mutant spread of
scores versus the within-mutant spread across conditions — if the within-mutant
condition spread rivals the between-mutant spread, conditions matter as much as
the mutation, which is the project's thesis.

Scores are averages of a {-1,0,+1} label, so these are effect-level statements,
not continuous-kinetics claims.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

DATA = Path("data/conditions_familial.csv")


def load(path: Path = DATA) -> pd.DataFrame:
    # "None" is a real category (no added salt, no shaking) and is in pandas'
    # default NA list, so it must not be parsed as missing.
    df = pd.read_csv(path, keep_default_na=False, na_values=[""])
    return df[df["support"] > 0].copy()  # a cell with n=0 carries no measurement


def _wls_r2(df: pd.DataFrame, formula: str) -> float:
    m = smf.wls(formula, data=df, weights=df["support"]).fit()
    return float(m.rsquared)


def variance_explained(df: pd.DataFrame) -> pd.DataFrame:
    """Per condition table: weighted R² for mutant / condition / both.

    There is deliberately no pooled row. The four tables are four partitions of
    the SAME experiments (each one has a buffer, a salt, a shaking speed and a
    concentration), so stacking them counts each experiment up to four times.
    """
    rows = []
    for dim, g in df.groupby("dimension"):
        if g["category"].nunique() < 2 or g["mutant"].nunique() < 2:
            continue
        r_mut = _wls_r2(g, "score ~ C(mutant)")
        r_cat = _wls_r2(g, "score ~ C(category)")
        r_both = _wls_r2(g, "score ~ C(mutant) + C(category)")
        rows.append({
            "dimension": dim,
            "n_cells": len(g),
            "r2_sequence_mutant": round(r_mut, 3),
            "r2_conditions_category": round(r_cat, 3),
            "r2_both": round(r_both, 3),
        })
    return pd.DataFrame(rows)


def spread_comparison(df: pd.DataFrame) -> pd.DataFrame:
    """Between-mutant vs within-mutant (across-condition) spread, per dimension.

    Support-weighted. within/between > ~1 means conditions move the score as
    much as the mutation does.
    """
    rows = []
    for dim, g in df.groupby("dimension"):
        mut_means = g.groupby("mutant").apply(
            lambda x: np.average(x["score"], weights=x["support"]), include_groups=False)
        between = _wstd(mut_means.values, np.ones(len(mut_means)))
        within = []
        wts = []
        for _, x in g.groupby("mutant"):
            if len(x) > 1:
                mu = np.average(x["score"], weights=x["support"])
                within.append(_wstd(x["score"].values, x["support"].values))
                wts.append(x["support"].sum())
        within_sd = float(np.average(within, weights=wts)) if within else np.nan
        rows.append({
            "dimension": dim,
            "between_mutant_sd": round(float(between), 3),
            "within_mutant_sd": round(within_sd, 3),
            "within_over_between": round(within_sd / between, 2) if between else np.nan,
        })
    return pd.DataFrame(rows)


def _wstd(x, w):
    x = np.asarray(x, float); w = np.asarray(w, float)
    mu = np.average(x, weights=w)
    return float(np.sqrt(np.average((x - mu) ** 2, weights=w)))


def a30p_buffer_example(df: pd.DataFrame) -> dict:
    """Reproduce the paper's headline: A30P scores -0.75 (PBS) vs -0.31 (Tris)."""
    a = df[(df.mutant == "A30P") & (df.dimension == "buffer")].set_index("category")
    return {c: float(a.loc[c, "score"]) for c in a.index}


def run() -> dict:
    df = load()
    return {
        "variance_explained": variance_explained(df),
        "spread": spread_comparison(df),
        "a30p_buffer": a30p_buffer_example(df),
        "n_measurements": int(df["support"].sum()),
    }


if __name__ == "__main__":
    out = run()
    print("Variance explained (support-weighted R²):")
    print(out["variance_explained"].to_string(index=False))
    print("\nBetween-mutant vs within-mutant spread:")
    print(out["spread"].to_string(index=False))
    print(f"\nA30P buffer scores: {out['a30p_buffer']}")
    print(f"Total measurements behind these cells: {out['n_measurements']}")


# --- sequence signal, from Tables S3–S6 -------------------------------------
TYPES_CSV = Path("data/type_scores.csv")
# The 16 substitution transitions (the "Cntrl" rows are same-class controls).
TRANSITIONS = [t for t in (
    "+ -> -", "- -> +", "nP -> P", "P -> nP", "Sm -> Lg", "Lg -> Sm",
    "nB -> B", "B -> nB", "nβ -> β", "β -> nβ",
    "Pro -> X", "X -> Pro", "Gly -> X", "X -> Gly", "PG -> X", "X -> PG")]
DOMAINS = ["N-terminus", "NAC", "C-terminus"]


def load_types(path: Path = TYPES_CSV) -> pd.DataFrame:
    return pd.read_csv(path)


def domain_means(df: pd.DataFrame) -> pd.DataFrame:
    """Support-weighted mean aggregation score per sequence region."""
    rows = []
    for scope, g in df.groupby("scope"):
        g = g[g["average"].notna() & (g["support"] > 0)]
        rows.append({"scope": scope,
                     "weighted_mean_score": round(
                         float(np.average(g["average"], weights=g["support"])), 3),
                     "measurements": int(g["support"].sum())})
    return pd.DataFrame(rows).sort_values("weighted_mean_score")


def type_by_domain(df: pd.DataFrame) -> pd.DataFrame:
    """16 substitution transitions x 3 regions, as a score matrix."""
    d = df[df["type"].isin(TRANSITIONS) & df["scope"].isin(DOMAINS)]
    return d.pivot(index="type", columns="scope", values="average").reindex(
        TRANSITIONS)[DOMAINS]


def permutation_test(df: pd.DataFrame, n_perm: int = 1000, seed: int = 0) -> pd.DataFrame:
    """Compare each factor's R² with the R² it reaches on shuffled labels.

    Raw R² is not comparable between the two factors here: mutant has 5 degrees
    of freedom against a category's 1-2, on only 14-19 cells, so the bigger
    factor explains more variance mechanically. Permuting the labels gives each
    factor its own chance baseline, which is the fair comparison.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for dim, g in df.groupby("dimension"):
        obs_m = _wls_r2(g, "score ~ C(mutant)")
        obs_c = _wls_r2(g, "score ~ C(category)")
        null_m, null_c = [], []
        for _ in range(n_perm):
            h = g.copy()
            h["mutant"] = rng.permutation(h["mutant"].values)
            h["category"] = rng.permutation(h["category"].values)
            null_m.append(_wls_r2(h, "score ~ C(mutant)"))
            null_c.append(_wls_r2(h, "score ~ C(category)"))
        null_m, null_c = np.array(null_m), np.array(null_c)
        rows.append({
            "dimension": dim,
            "n_cells": len(g),
            "r2_mutant": round(obs_m, 3),
            "chance_mutant": round(float(null_m.mean()), 3),
            "p_mutant": round(float((np.sum(null_m >= obs_m) + 1) / (n_perm + 1)), 4),
            "r2_category": round(obs_c, 3),
            "chance_category": round(float(null_c.mean()), 3),
            "p_category": round(float((np.sum(null_c >= obs_c) + 1) / (n_perm + 1)), 4),
        })
    return pd.DataFrame(rows)
