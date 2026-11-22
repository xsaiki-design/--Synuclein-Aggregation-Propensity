"""Where does the published literature agree about alpha-synuclein mutations,
and where does it contradict itself?

The project's premise was that the literature disagrees with itself — A30P has
been reported to slow aggregation, to accelerate it, and to have no effect —
and that the disagreement is the finding. This measures it from data already
extracted.

Pancoe's Tables S3-S6 report, for each (region x mutation type), the mean
aggregation score, its standard deviation, and the number of underlying
experiments. Scores are means of per-experiment labels in {-1, 0, +1}, so
those three numbers determine the vote exactly:

    E[x]   = (n+ - n-) / n = m
    E[x^2] = (n+ + n-) / n = s^2 + m^2        (x^2 is 1 unless the label is 0)

    => n+ = n(s^2 + m^2 + m)/2,  n- = n(s^2 + m^2 - m)/2,  n0 = n - n+ - n-

Each row reads back as "k experiments reported acceleration, j reported
reduction". That is algebra, not modelling, and the recovered counts must be
non-negative near-integers — a further check on the extraction.

PUBLICATION BIAS. These are counts of *published* experiments. A result that
matches a strong prior expectation is easier to publish than a null or
contradictory one, so unanimity in this data is consistent with a real effect,
with selective publication, or with both. Nothing here can separate them.
Report these as counts of published reports, never as counts of experiments
performed.

An earlier version of this module scored dispersion as sd / sqrt(1 - mean^2).
That ratio was removed: it labelled X -> Pro in the NAC "mixed" when the vote
is 16-0-1, because the attainable spread shrinks at extreme means. The vote
counts carry the same information without that failure mode.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .si_analysis import DOMAINS, TRANSITIONS, load_types

# Scores are means of {-1,0,+1} labels, so sd cannot exceed sqrt(1-mean^2).
_TOL = 0.02  # rounding slack: the SI prints 3 decimals


def _max_sd(mean) -> np.ndarray:
    return np.sqrt(np.clip(1.0 - np.asarray(mean, float) ** 2, 0.0, None))


def validate(df: pd.DataFrame) -> pd.DataFrame:
    """Rows whose SD exceeds what a {-1,0,+1} mean allows.

    An independent check on the extraction: any violation means a mean or SD
    was decoded wrongly, since no real pair can breach this bound.
    """
    d = df[df["average"].notna() & (df["support"] > 0)].copy()
    d["max_sd"] = _max_sd(d["average"])
    return d[d["sd"] > d["max_sd"] + _TOL]


def vote_split(df: pd.DataFrame | None = None, min_support: int = 5) -> pd.DataFrame:
    """Recover how many published experiments reported each direction.

    `split` is the minority share of the directional votes: 0 means every
    report agreed on the sign, 0.5 means an even split.
    """
    d = (load_types() if df is None else df).copy()
    d = d[d["average"].notna() & (d["support"] >= min_support)
          & d["type"].isin(TRANSITIONS)].copy()
    n, m, s = d["support"].astype(float), d["average"], d["sd"]
    nz = s ** 2 + m ** 2                      # fraction with a non-zero score
    d["n_accel"] = n * (nz + m) / 2
    d["n_reduce"] = n * (nz - m) / 2
    d["n_noeffect"] = n - d["n_accel"] - d["n_reduce"]
    d["integrality_error"] = np.max(np.abs(
        np.vstack([d[c] - d[c].round() for c in ("n_accel", "n_reduce", "n_noeffect")])
    ), axis=0).round(3)
    for c in ("n_accel", "n_reduce", "n_noeffect"):
        d[c] = d[c].round().astype(int).clip(lower=0)
    directional = d["n_accel"] + d["n_reduce"]
    d["split"] = np.where(directional > 0,
                          np.minimum(d["n_accel"], d["n_reduce"]) / directional, 0.0)
    d["split"] = d["split"].round(3)
    return (d[["scope", "type", "support", "average", "n_accel", "n_noeffect",
               "n_reduce", "split", "integrality_error"]]
            .sort_values(["split", "support"], ascending=[False, False])
            .reset_index(drop=True))


def by_region(df: pd.DataFrame | None = None, min_support: int = 5) -> pd.DataFrame:
    """Is the published record more self-consistent about some regions?

    Only the three regions are returned. The whole-protein rows of Table S3 are
    the sum of the three regional tables — verified exact for all 26 types — so
    listing them alongside would invite double counting, the same nesting error
    as pooling the four condition tables.
    """
    v = vote_split(df, min_support)
    rows = []
    for scope in DOMAINS:                     # regions only; see docstring
        g = v[v.scope == scope]
        if g.empty:
            continue
        rows.append({
            "scope": scope,
            "types": len(g),
            "experiments": int(g.support.sum()),
            "mean_split": round(float(np.average(g.split, weights=g.support)), 3),
            "unanimous": int((g.split == 0).sum()),
            "deadlocked": int((g.split >= 0.4).sum()),
        })
    return pd.DataFrame(rows)


def headline_robustness(df: pd.DataFrame | None = None) -> pd.DataFrame:
    """Do this project's headline findings rest on an agreed-upon record?"""
    v = vote_split(df, min_support=5)
    want = [("NAC", "X -> Pro"), ("N-terminus", "X -> Pro"),
            ("whole_protein", "X -> Pro"), ("NAC", "X -> PG"),
            ("whole_protein", "nB -> B"), ("N-terminus", "nP -> P")]
    rows = []
    for scope, typ in want:
        g = v[(v.scope == scope) & (v["type"] == typ)]
        if g.empty:
            continue
        r = g.iloc[0]
        rows.append({
            "finding": f"{typ} in {scope}",
            "accel": int(r.n_accel), "none": int(r.n_noeffect),
            "reduce": int(r.n_reduce), "split": r.split,
            "verdict": ("unanimous in direction" if r.split == 0
                        else "deadlocked" if r.split >= 0.4 else "majority"),
        })
    return pd.DataFrame(rows)


def run() -> dict:
    t = load_types()
    return {"violations": validate(t), "votes": vote_split(t),
            "by_region": by_region(t), "headlines": headline_robustness(t)}


if __name__ == "__main__":
    r = run()
    print(f"extraction sanity — rows breaching the sd bound: {len(r['violations'])}\n")
    print("MOST SPLIT (published reports pointing opposite ways):")
    print(r["votes"].head(6).to_string(index=False))
    print("\nUNANIMOUS IN DIRECTION:")
    print(r["votes"][r["votes"].split == 0].head(6).to_string(index=False))
    print("\nBY REGION (regions only — whole_protein is their sum):")
    print(r["by_region"].to_string(index=False))
    print("\nHEADLINES (counts are of PUBLISHED reports):")
    print(r["headlines"].to_string(index=False))
