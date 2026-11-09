"""Canonical sequence, domain boundaries, and dataset loading.

The schema here is the contract every other module relies on. `ingest.py`
produces a CSV in this shape from the real supplementary PDF; `synth.py`
produces a clearly-labelled synthetic stand-in in the same shape so the whole
pipeline runs before the real data is in place.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

# Canonical human alpha-synuclein (SNCA), UniProt P37840, 140 aa.
# Source: https://www.uniprot.org/uniprotkb/P37840  (verify on first ingest)
ASYN_SEQUENCE = (
    "MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQVTNVGGAVVTGVTAVAQKT"
    "VEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA"
)
assert len(ASYN_SEQUENCE) == 140, f"expected 140 aa, got {len(ASYN_SEQUENCE)}"

# Domain boundaries (1-based, inclusive), the standard partition used in the
# literature. NAC is the aggregation-prone core.
DOMAINS = {
    "amphipathic": (1, 60),   # N-terminal, lipid-binding
    "nac": (61, 95),          # non-amyloid-beta component, fibril core
    "acidic": (96, 140),      # C-terminal, highly acidic
}
NAC_START, NAC_END = DOMAINS["nac"]

# The columns ingest/synth must produce. Anything else is optional metadata.
REQUIRED_COLUMNS = [
    "variant",          # e.g. "A30P", or "WT"
    "study",            # grouping key: the paper the row was measured in
    "position",         # residue number (1-based), NaN for WT
    "wt_residue",       # single-letter, NaN for WT
    "mut_residue",      # single-letter, NaN for WT
    "t_half_h",         # outcome: time to half-maximal ThT fluorescence (hours)
    "t_half_rel_wt",    # t_half normalised to the WT control in the same study
    # experimental conditions (the independent variable of the main experiment):
    "pH",
    "temperature_C",
    "protein_conc_uM",
    "agitation",        # categorical: "quiescent" | "shaking" | "beads"
    "buffer",           # categorical
    "seeds",            # 0/1
    "lipids",           # 0/1
    "ThT_conc_uM",
]

CONDITION_COLUMNS = [
    "pH", "temperature_C", "protein_conc_uM",
    "agitation", "buffer", "seeds", "lipids", "ThT_conc_uM",
]

CURATED_CSV = Path("data/aggregation_curated.csv")


def domain_of(position: float) -> str:
    """Return the domain name for a 1-based residue position (or 'wt' if NaN)."""
    if position is None or pd.isna(position):
        return "wt"
    pos = int(position)
    for name, (lo, hi) in DOMAINS.items():
        if lo <= pos <= hi:
            return name
    raise ValueError(f"position {pos} outside 1..140")


def distance_to_nac(position: float) -> float:
    """Residue distance to the NAC domain; 0 if inside it, NaN for WT."""
    if position is None or pd.isna(position):
        return float("nan")
    pos = int(position)
    if NAC_START <= pos <= NAC_END:
        return 0.0
    return float(NAC_START - pos if pos < NAC_START else pos - NAC_END)


def load_curated(path: Path | str = CURATED_CSV) -> pd.DataFrame:
    """Load the curated dataset and validate the schema.

    Raises a clear error if the file is missing, pointing at how to make it.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found. Build it first:\n"
            "  make data        # from the real supplementary PDF (data/raw/SI.pdf)\n"
            "  make data-synth  # synthetic stand-in, to run the pipeline without the PDF"
        )
    df = pd.read_csv(path, comment="#")  # '#' header marks the synthetic stand-in
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"{path} is missing required columns: {missing}")
    return df
