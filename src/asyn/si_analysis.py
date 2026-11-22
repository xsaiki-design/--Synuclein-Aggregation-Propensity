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


def domain_means(df: pd.DataFrame, transitions_only: bool = False) -> pd.DataFrame:
    """Support-weighted mean aggregation score per sequence region.

    By default this includes the ten same-class "Cntrl" rows (e.g. nP -> nP)
    alongside the 16 class-changing transitions. That choice matters: on
    transitions alone the N-terminus and C-terminus are indistinguishable
    (-0.19 vs -0.20), whereas including the controls separates them
    (-0.10 vs -0.35). The NAC being the most affected region is robust either
    way (-0.57 with controls, -0.70 without) and across support thresholds
    1 to 10; the N-terminus/C-terminus ordering is not.
    """
    if transitions_only:
        df = df[df["type"].isin(TRANSITIONS)]
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


def strongest_by_region(df: pd.DataFrame, min_support: int = 5) -> pd.DataFrame:
    """Strongest substitution within each sequence region.

    The whole-protein rows of Table S3 aggregate across regions, so a figure
    like "X -> Pro, -0.88, n=33" is not a regional result: those 33 split 16 in
    the N-terminus and 17 in the NAC, with none in the C-terminus. Rank within
    each region instead, and apply a support floor so a single unanimous
    measurement cannot be called the strongest effect.

    SELECTION EFFECT. This returns the maximum of 5-7 ranked candidates per
    region, and the maximum of a small sample is biased away from zero. On this
    data the N-terminus and C-terminus winners both land on exactly +/-1.000
    with sd = 0 at n = 5 and n = 7. Treat those as order statistics, not effect
    estimates; only the NAC row (n = 17) has the support to be quoted. Use
    `region_ranking` to see every candidate.
    """
    d = df[df["type"].isin(TRANSITIONS) & df["average"].notna()
           & (df["support"] >= min_support) & df["scope"].isin(DOMAINS)].copy()
    d["abs_effect"] = d["average"].abs()
    d = d.sort_values(["abs_effect", "support"], ascending=[False, False])
    top = d.groupby("scope", sort=False).head(1)
    return (top.set_index("scope").reindex(DOMAINS).reset_index()
              [["scope", "type", "average", "sd", "support"]])


def region_ranking(df: pd.DataFrame, scope: str) -> pd.DataFrame:
    """Every substitution with data in one region, ranked by absolute effect."""
    d = df[(df["scope"] == scope) & df["type"].isin(TRANSITIONS)
           & df["average"].notna() & (df["support"] > 0)].copy()
    d["abs_effect"] = d["average"].abs()
    return (d.sort_values(["abs_effect", "support"], ascending=[False, False])
             [["type", "average", "sd", "support"]].reset_index(drop=True))


def variance_explained_bounds(df: pd.DataFrame | None = None) -> pd.DataFrame:
    """Bracket the mutant-identity R² at the level of individual experiments.

    `variance_explained` and `permutation_test` fit cell *means*: each of the
    14-19 points is already an average of up to 16 experiments, so disagreement
    inside a cell is averaged away before the variance is computed. That makes
    the reported R² an upper bound on the experiment-level figure.

    S7-S10 give no standard deviations, so within-cell variance is unknown. It
    is bounded though: for labels in {-1,0,+1} with cell mean m it cannot exceed
    1 - m^2. Putting it at that maximum gives the lower bound. The true value
    lies between the two columns and cannot be resolved from this source.
    """
    d = load() if df is None else df
    rows = []
    for dim, g in d.groupby("dimension"):
        n = g["support"].values.astype(float)
        m = g["score"].values
        grand = np.average(m, weights=n)
        mut = g.groupby("mutant").apply(
            lambda x: np.average(x["score"], weights=x["support"]),
            include_groups=False)
        mm = g["mutant"].map(mut).values
        ss_model = float(np.sum(n * (mm - grand) ** 2))
        ss_between = float(np.sum(n * (m - grand) ** 2))
        ss_within_max = float(np.sum(n * (1 - m ** 2)))
        rows.append({
            "dimension": dim,
            "cells": len(g),
            "experiments": int(n.sum()),
            "r2_cell_means": round(ss_model / ss_between, 3),
            "r2_experiment_lower": round(ss_model / (ss_between + ss_within_max), 3),
        })
    return pd.DataFrame(rows)
