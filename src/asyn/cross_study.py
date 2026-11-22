"""Compare the familial-mutant effects reported by different studies.

Studies measure different things under different conditions, so this module
never pools them into one number. It lines them up per mutant and asks whether
they agree on *direction*, which is the only comparison their differing
endpoints support.

  Pancoe2022      categorical score in {-1,0,+1}, averaged, normalised to a
                  matched WT control within each study; predominantly bulk
                  ThT assays.
  Flagmeier2016   log10 fold-change vs WT, separated by mechanism, measured
                  with lipid vesicles (nucleation) or seeds (amplification).

The headline comparison is Pancoe's aggregate score against Flagmeier's
lipid-induced nucleation, because those are the two "does this mutation make
aggregation happen sooner" measurements.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from . import si_analysis

CROSS_CSV = Path("data/cross_study_effects.csv")


def load_cross(path: Path = CROSS_CSV) -> pd.DataFrame:
    return pd.read_csv(path, keep_default_na=False, na_values=[""])


def pancoe_per_mutant() -> pd.Series:
    """Support-weighted mean categorical score per familial mutant."""
    d = si_analysis.load()
    return d.groupby("mutant").apply(
        lambda x: np.average(x["score"], weights=x["support"]),
        include_groups=False)


def compare(endpoint: str = "lipid_nucleation") -> pd.DataFrame:
    """Line up Pancoe's score against one Flagmeier mechanism, per mutant."""
    cross = load_cross()
    f = cross[(cross.study == "Flagmeier2016") & (cross.endpoint == endpoint)]
    p = pancoe_per_mutant()
    rows = []
    for _, r in f.iterrows():
        if r.mutant not in p.index:
            continue
        ps = float(p[r.mutant])
        rows.append({
            "mutant": r.mutant,
            "pancoe_score": round(ps, 2),
            "flagmeier_log10": r.value,
            "flagmeier_censored": bool(r.censored),
            "same_direction": _same_sign(ps, r.value),
        })
    return pd.DataFrame(rows).sort_values("mutant").reset_index(drop=True)


def _same_sign(a: float, b: float, tol: float = 1e-9) -> bool | None:
    """None when either side reports no effect — agreement is undefined then."""
    if abs(a) <= tol or abs(b) <= tol:
        return None
    return bool(np.sign(a) == np.sign(b))


def _binom_two_sided(k: int, n: int) -> float:
    """Exact two-sided binomial p against p=0.5 (no scipy dependency)."""
    from math import comb
    if n == 0:
        return float("nan")
    probs = [comb(n, i) / 2 ** n for i in range(n + 1)]
    return float(min(1.0, sum(q for q in probs if q <= probs[k] + 1e-12)))


def _spearman(a, b) -> float:
    """Spearman rho via ranked Pearson; NaN if either side is constant."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 3:
        return float("nan")
    ra, rb = pd.Series(a).rank().values, pd.Series(b).rank().values
    if ra.std() == 0 or rb.std() == 0:
        return float("nan")
    return float(np.corrcoef(ra, rb)[0, 1])


def agreement_summary() -> pd.DataFrame:
    """Directional agreement between Pancoe and each Flagmeier mechanism."""
    out = []
    for ep in ["lipid_nucleation", "secondary_nucleation", "elongation"]:
        c = compare(ep)
        decided = c[c.same_direction.notna()]
        k, n = int(decided.same_direction.sum()), len(decided)
        p = _binom_two_sided(k, n)
        out.append({
            "flagmeier_endpoint": ep,
            "mutants_compared": len(c),
            "decidable": n,
            "agree": k,
            "disagree": n - k,
            "spearman_rho": round(_spearman(c.pancoe_score, c.flagmeier_log10), 3),
            "binom_p": round(p, 3),
            # With a handful of mutants, sign agreement alone proves nothing.
            "distinguishable_from_chance": bool(p < 0.05),
        })
    return pd.DataFrame(out)


def claim_checks() -> pd.DataFrame:
    """Which of the project's three claims this new source can actually test.

    Two of them it cannot touch, and saying so is part of the result.
    """
    cross = load_cross()
    f_mutants = sorted(cross[cross.study == "Flagmeier2016"].mutant.unique())
    # Every familial mutant studied sits in the N-terminus (residues 1-60).
    positions = {"A30P": 30, "E46K": 46, "H50Q": 50, "G51D": 51,
                 "A53E": 53, "A53T": 53}
    all_n_term = all(positions[m] <= 60 for m in f_mutants if m in positions)
    has_proline = any(m.endswith("P") for m in f_mutants)
    agree = agreement_summary()
    lipid = agree[agree.flagmeier_endpoint == "lipid_nucleation"].iloc[0]
    return pd.DataFrame([
        {"claim": "NAC is the most impactful region",
         "testable_here": False,
         "why": f"every mutant studied is N-terminal ({', '.join(f_mutants)}); "
                "no NAC or C-terminal variant is covered"},
        {"claim": "X -> Pro is the strongest single substitution",
         "testable_here": False,
         "why": "only one proline variant is covered (A30P), so the ranking "
                "cannot be reproduced. That single point disagrees: Pancoe "
                "scores N-terminal X -> Pro at -0.81, while this source reports "
                "no notable effect of A30P on lipid-induced nucleation (n=1)"},
        {"claim": "sequence matters more than condition",
         "testable_here": True,
         "why": f"directly contradicted across assay families: only "
                f"{int(lipid.agree)}/{int(lipid.decidable)} "
                "mutants keep their direction between bulk-assay scores and "
                f"lipid-induced nucleation (binomial p={lipid.binom_p}, so the "
                "split itself is not yet distinguishable from chance)"},
    ])


def run() -> dict:
    return {"comparison": compare(),
            "agreement": agreement_summary(),
            "claims": claim_checks()}


if __name__ == "__main__":
    r = run()
    print("Pancoe score vs Flagmeier lipid-induced nucleation:")
    print(r["comparison"].to_string(index=False))
    print("\nDirectional agreement by mechanism:")
    print(r["agreement"].to_string(index=False))
    print("\nWhat this source can test:")
    print(r["claims"].to_string(index=False))


def mechanism_spread(study: str = "Flagmeier2016") -> pd.DataFrame:
    """How much the same mutations differ, mechanism by mechanism, in ONE study.

    This needs no cross-study comparison, so it is immune to the endpoint- and
    condition-mismatch problems that make the Pancoe/Flagmeier agreement
    uninformative. One lab, one protein prep, one set of variants: any spread
    is a sequence effect on that mechanism.
    """
    d = load_cross()
    d = d[d.study == study]
    rows = []
    for ep, g in d.groupby("endpoint"):
        rows.append({
            "endpoint": ep,
            "n_mutants": len(g),
            "min_log10": g.value.min(),
            "max_log10": g.value.max(),
            "spread_log10": round(float(g.value.max() - g.value.min()), 2),
            "fold_range": f"{10 ** float(g.value.max() - g.value.min()):,.0f}x",
            "any_censored": bool(g.censored.any()),
        })
    return (pd.DataFrame(rows)
            .sort_values("spread_log10", ascending=False).reset_index(drop=True))


PANCOE_REFS = Path("data/pancoe_references.csv")


def cited_by_pancoe(doi: str) -> bool:
    """Is this DOI in Pancoe's deposited reference list?

    The independence test for any candidate source. Pancoe is a meta-analysis
    of 375 published examples, so a paper it cites may already be inside the
    compilation; treating it as corroboration then risks counting the same
    underlying experiments twice.

    A citation does not prove the paper's measurements were tabulated — it may
    be cited only in discussion — but it removes the presumption of
    independence, which is the direction that matters here.

    Reference list retrieved by DOI content negotiation (doi.org redirects to
    api.crossref.org); cached in data/pancoe_references.csv.
    """
    if not PANCOE_REFS.exists():
        raise FileNotFoundError(f"{PANCOE_REFS} missing; see cross_study docstring")
    refs = pd.read_csv(PANCOE_REFS, keep_default_na=False)
    return doi.strip().lower() in set(refs["doi"].str.strip().str.lower())
