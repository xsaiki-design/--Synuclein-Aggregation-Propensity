"""Where does the literature agree about alpha-synuclein mutations, and where
does it contradict itself?

The project's premise was that the literature openly disagrees with itself —
A30P has been reported to slow aggregation, to accelerate it, and to have no
effect — and that the disagreement is the finding. This module measures that
directly, from data already extracted.

Pancoe's Tables S3-S6 report, for each (region x mutation type), the mean
aggregation score, its **standard deviation**, and the number of underlying
experiments. That SD has been treated as noise everywhere else in this repo.
It is not noise: it is how much the published experiments disagreed with each
other about that kind of mutation.

The metric. Scores take values in {-1, 0, +1}, so the spread attainable at a
given mean is bounded. A distribution on {-1, 0, +1} with mean m has maximum
variance 1 - m^2, reached by putting all mass on +/-1 — i.e. by the literature
splitting into "accelerates" and "reduces" camps with nothing in between. So

    discordance = sd / sqrt(1 - mean^2)

runs from 0 (every experiment agrees) to 1 (maximally polarised at that mean).
Dividing by the attainable maximum matters: a type with mean 0 can reach sd 1
while a type with mean 0.9 cannot exceed 0.44, so raw SD would make
strong-effect types look spuriously consistent.
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


def score(df: pd.DataFrame | None = None, min_support: int = 5) -> pd.DataFrame:
    """Per (region x type) discordance, for rows with enough experiments."""
    d = (load_types() if df is None else df).copy()
    d = d[d["average"].notna() & (d["support"] >= min_support)
          & d["type"].isin(TRANSITIONS)].copy()
    mx = _max_sd(d["average"])
    # At |mean| = 1 every experiment agreed, so discordance is 0 by definition.
    safe = np.where(mx <= _TOL, 1.0, mx)
    d["discordance"] = np.where(mx <= _TOL, 0.0, d["sd"] / safe)
    d["discordance"] = d["discordance"].clip(0.0, 1.0).round(3)
    d["verdict"] = np.where(d["discordance"] >= 0.75, "contested",
                     np.where(d["discordance"] <= 0.35, "consensus", "mixed"))
    return (d[["scope", "type", "average", "sd", "support", "discordance", "verdict"]]
            .sort_values(["discordance", "support"], ascending=[False, False])
            .reset_index(drop=True))


def by_region(df: pd.DataFrame | None = None, min_support: int = 5) -> pd.DataFrame:
    """Is the literature more self-consistent about some regions than others?"""
    s = score(df, min_support)
    rows = []
    for scope in DOMAINS + ["whole_protein"]:
        g = s[s.scope == scope]
        if g.empty:
            continue
        rows.append({
            "scope": scope,
            "types": len(g),
            "experiments": int(g.support.sum()),
            "mean_discordance": round(float(np.average(g.discordance,
                                                       weights=g.support)), 3),
            "contested": int((g.verdict == "contested").sum()),
            "consensus": int((g.verdict == "consensus").sum()),
        })
    return pd.DataFrame(rows)


def headline_robustness(df: pd.DataFrame | None = None) -> pd.DataFrame:
    """Do this project's own headline findings rest on agreed-upon data?

    A large mean effect built on a contested literature is a weaker claim than
    the number alone suggests, so the claims are checked against their own
    discordance rather than reported as settled.
    """
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
            "reduce": int(r.n_reduce),
            "split": r.split,
            "verdict": ("unanimous in direction" if r.split == 0
                        else "deadlocked" if r.split >= 0.4 else "majority"),
        })
    return pd.DataFrame(rows)


def run() -> dict:
    t = load_types()
    return {"violations": validate(t), "scored": score(t),
            "by_region": by_region(t), "headlines": headline_robustness(t)}




def vote_split(df: pd.DataFrame | None = None, min_support: int = 5) -> pd.DataFrame:
    """Recover how many experiments said accelerate / no effect / reduce.

    For scores in {-1, 0, +1} the group mean, SD and count determine the vote
    exactly. With n experiments, mean m and SD s:

        E[x]   = (n+ - n-) / n = m
        E[x^2] = (n+ + n-) / n = s^2 + m^2     (since x^2 is 1 unless x is 0)

    so  n+ = n(s^2 + m^2 + m)/2,  n- = n(s^2 + m^2 - m)/2,  n0 = n - n+ - n-.

    This is more interpretable than any dispersion index — it says "5 papers
    found acceleration and 4 found reduction" — and the reconstructed counts
    must come out as non-negative near-integers, which is a further independent
    check on the extraction.
    """
    d = (load_types() if df is None else df).copy()
    d = d[d["average"].notna() & (d["support"] >= min_support)
          & d["type"].isin(TRANSITIONS)].copy()
    n, m, s = d["support"].astype(float), d["average"], d["sd"]
    nz = s ** 2 + m ** 2                      # fraction with a non-zero score
    d["n_accel"] = (n * (nz + m) / 2)
    d["n_reduce"] = (n * (nz - m) / 2)
    d["n_noeffect"] = n - d["n_accel"] - d["n_reduce"]
    # How far the reconstruction sits from whole experiments.
    d["integrality_error"] = np.max(np.abs(
        np.vstack([d[c] - d[c].round() for c in ("n_accel", "n_reduce", "n_noeffect")])
    ), axis=0).round(3)
    for c in ("n_accel", "n_reduce", "n_noeffect"):
        d[c] = d[c].round().astype(int).clip(lower=0)
    # Minority share: 0 = unanimous on direction, 0.5 = evenly split.
    directional = d["n_accel"] + d["n_reduce"]
    d["split"] = np.where(directional > 0,
                          np.minimum(d["n_accel"], d["n_reduce"]) / directional, 0.0)
    d["split"] = d["split"].round(3)
    return (d[["scope", "type", "support", "average", "n_accel", "n_noeffect",
               "n_reduce", "split", "integrality_error"]]
            .sort_values(["split", "support"], ascending=[False, False])
            .reset_index(drop=True))

if __name__ == "__main__":
    r = run()
    print(f"extraction sanity — rows breaching the sd bound: {len(r['violations'])}\n")
    s = r["scored"]
    print("MOST CONTESTED (literature splits):")
    print(s.head(8).to_string(index=False))
    print("\nMOST CONSISTENT:")
    print(s.tail(6).to_string(index=False))
    print("\nBY REGION:")
    print(r["by_region"].to_string(index=False))
    print("\nMOST SPLIT (papers pointing opposite ways):")
    print(vote_split().head(5).to_string(index=False))
    print("\nDO OUR OWN HEADLINES REST ON AGREEMENT?")
    print(r["headlines"].to_string(index=False))
