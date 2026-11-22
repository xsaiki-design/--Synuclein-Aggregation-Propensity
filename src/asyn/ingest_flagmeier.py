"""Flagmeier et al. 2016 — continuous, mechanism-separated effects.

  Flagmeier P, Meisl G, Vendruscolo M, Knowles TPJ, Dobson CM, Buell AK,
  Galvagnion C. "Mutations associated with familial Parkinson's disease alter
  the initiation and amplification steps of alpha-synuclein aggregation."
  PNAS 2016; 113(37):10328-10333. doi:10.1073/pnas.1604645113
  Open access: https://pmc.ncbi.nlm.nih.gov/articles/PMC5027465

Why this source matters: it covers the same familial mutants as the Pancoe
compilation, but reports effects as orders of magnitude relative to WT and
separates them by mechanism (lipid-induced primary nucleation, secondary
nucleation, elongation) instead of collapsing everything into one label.

IMPORTANT — what these values are. The paper states its effects in prose as
approximate orders of magnitude; it does not tabulate rate constants. Each
value below is an integer log10 fold-change transcribed from a quoted sentence,
NOT a digitised measurement. The supporting quote is stored on every row so the
number can be audited against the source. Treat them as directional and
coarse: good for comparing signs and rough magnitudes, not for regression on
exact values.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

STUDY = "Flagmeier2016"
SOURCE = "https://pmc.ncbi.nlm.nih.gov/articles/PMC5027465"
OUT = Path("data/cross_study_effects.csv")

_LIPID = ("Our global analysis shows that the A30P mutation does not notably "
          "affect the lipid-induced aggregation rate, whereas the A53T mutation "
          "increases this rate by approximately two orders of magnitude relative "
          "to that of WT alpha-synuclein.")
_LIPID2 = ("The mutations E46K, H50Q, and G51D all decrease the lipid-induced "
           "aggregation rate, by one (H50Q) or two (E46K, G51D) orders of "
           "magnitude, relative to that of WT alpha-synuclein.")
_SEC = ("the rates of fibril amplification for A53T, A30P, E46K, and WT "
        "alpha-synuclein were observed to be within an order of magnitude of "
        "each other, whereas the rates of secondary nucleation were found to be "
        "undetectably slow for H50Q and G51D")
_ELONG = ("such normalized elongation rates varied by less than one order of "
          "magnitude, except for the value determined for monomeric G51D, which "
          "elongates the fibrils of different variants notably more slowly")

# (mutant, endpoint, condition_class, log10 fold vs WT, supporting quote)
ROWS = [
    ("A30P", "lipid_nucleation", "lipid_vesicles", 0.0, _LIPID),
    ("A53T", "lipid_nucleation", "lipid_vesicles", 2.0, _LIPID),
    ("E46K", "lipid_nucleation", "lipid_vesicles", -2.0, _LIPID2),
    ("H50Q", "lipid_nucleation", "lipid_vesicles", -1.0, _LIPID2),
    ("G51D", "lipid_nucleation", "lipid_vesicles", -2.0, _LIPID2),
    ("A30P", "secondary_nucleation", "seeded_acidic", 0.0, _SEC),
    ("E46K", "secondary_nucleation", "seeded_acidic", 0.0, _SEC),
    ("A53T", "secondary_nucleation", "seeded_acidic", 0.0, _SEC),
    # "undetectably slow" is a bound, not a value; see CENSORED below.
    ("H50Q", "secondary_nucleation", "seeded_acidic", -2.0, _SEC),
    ("G51D", "secondary_nucleation", "seeded_acidic", -2.0, _SEC),
    ("A30P", "elongation", "seeded", 0.0, _ELONG),
    ("E46K", "elongation", "seeded", 0.0, _ELONG),
    ("H50Q", "elongation", "seeded", 0.0, _ELONG),
    ("A53T", "elongation", "seeded", 0.0, _ELONG),
    ("G51D", "elongation", "seeded", -1.0, _ELONG),
]
# Rows the paper reports as "undetectably slow" rather than quantifying: the
# -2.0 is a floor standing in for "at least this slow", not a measurement.
CENSORED = {("H50Q", "secondary_nucleation"), ("G51D", "secondary_nucleation")}


def build() -> pd.DataFrame:
    """Flagmeier's stated effects, one row per (mutant, mechanism)."""
    return pd.DataFrame([{
        "study": STUDY, "year": 2016, "source": SOURCE,
        "mutant": m, "endpoint": e, "condition_class": c,
        "normalized_to_wt": True,        # stated relative to WT
        "value": v, "value_type": "log10_fold_vs_wt",
        "censored": (m, e) in CENSORED,
        "extraction": "text_statement",  # prose, not a table
        "quote": q,
    } for m, e, c, v, q in ROWS])


if __name__ == "__main__":
    df = build()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT, index=False)
    print(f"wrote {OUT} ({len(df)} rows) [Flagmeier 2016, prose-stated magnitudes]")
