"""Three feature blocks, kept separate so they can be tested independently:

  1. sequence   — physicochemical change on mutation + positional context
  2. structure  — per-residue fibril structure context (SASA, SS, core)
  3. conditions — the experimental columns (the main experiment's predictor)

Encoding and scaling are deliberately NOT done here. They must happen inside
each cross-validation fold (see methodology guardrail: scale inside the fold),
so this module emits raw numeric/categorical columns and the model pipeline
does the fold-local transforms.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .data import CONDITION_COLUMNS, distance_to_nac, domain_of

# --- amino-acid property tables (standard public constants) ------------------
# Kyte-Doolittle hydropathy.
HYDRO = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}
# Net charge at ~pH 7.
CHARGE = {a: 0 for a in "ACFGILMNPQSTVWY"}
CHARGE.update({"D": -1, "E": -1, "K": 1, "R": 1, "H": 0.1})
# Side-chain volume (A^3), Zamyatnin.
VOLUME = {
    "A": 88.6, "R": 173.4, "N": 114.1, "D": 111.1, "C": 108.5, "Q": 143.8,
    "E": 138.4, "G": 60.1, "H": 153.2, "I": 166.7, "L": 166.7, "K": 168.6,
    "M": 162.9, "F": 189.9, "P": 112.7, "S": 89.0, "T": 116.1, "W": 227.8,
    "Y": 193.6, "V": 140.0,
}
# Chou-Fasman helix (P_a) and sheet (P_b) propensities.
HELIX = {
    "A": 1.42, "R": 0.98, "N": 0.67, "D": 1.01, "C": 0.70, "Q": 1.11, "E": 1.51,
    "G": 0.57, "H": 1.00, "I": 1.08, "L": 1.21, "K": 1.16, "M": 1.45, "F": 1.13,
    "P": 0.57, "S": 0.77, "T": 0.83, "W": 1.08, "Y": 0.69, "V": 1.06,
}
SHEET = {
    "A": 0.83, "R": 0.93, "N": 0.89, "D": 0.54, "C": 1.19, "Q": 1.10, "E": 0.37,
    "G": 0.75, "H": 0.87, "I": 1.60, "L": 1.30, "K": 0.74, "M": 1.05, "F": 1.38,
    "P": 0.55, "S": 0.75, "T": 1.19, "W": 1.37, "Y": 1.47, "V": 1.70,
}

SEQUENCE_COLS = [
    "d_hydro", "d_charge", "d_volume", "d_helix", "d_sheet",
    "position", "dist_to_nac", "dom_amphipathic", "dom_nac", "dom_acidic",
]
STRUCTURE_COLS = ["sasa", "ss_sheet", "in_fibril_core"]
CONDITION_COLS = list(CONDITION_COLUMNS)

STRUCTURE_CSV = Path("data/structure_features.csv")


def _delta(wt, mut, table):
    if pd.isna(wt) or pd.isna(mut):
        return 0.0  # WT row: no substitution
    return float(table[mut] - table[wt])


def sequence_block(df: pd.DataFrame) -> pd.DataFrame:
    out = pd.DataFrame(index=df.index)
    out["d_hydro"] = [_delta(w, m, HYDRO) for w, m in zip(df.wt_residue, df.mut_residue)]
    out["d_charge"] = [_delta(w, m, CHARGE) for w, m in zip(df.wt_residue, df.mut_residue)]
    out["d_volume"] = [_delta(w, m, VOLUME) for w, m in zip(df.wt_residue, df.mut_residue)]
    out["d_helix"] = [_delta(w, m, HELIX) for w, m in zip(df.wt_residue, df.mut_residue)]
    out["d_sheet"] = [_delta(w, m, SHEET) for w, m in zip(df.wt_residue, df.mut_residue)]
    out["position"] = df.position.fillna(0).astype(float)
    out["dist_to_nac"] = [distance_to_nac(p) for p in df.position]
    doms = [domain_of(p) for p in df.position]
    out["dom_amphipathic"] = [int(d == "amphipathic") for d in doms]
    out["dom_nac"] = [int(d == "nac") for d in doms]
    out["dom_acidic"] = [int(d == "acidic") for d in doms]
    out["dist_to_nac"] = out["dist_to_nac"].fillna(0.0)
    return out


def structure_block(df: pd.DataFrame, path: Path | str = STRUCTURE_CSV) -> pd.DataFrame:
    """Per-residue fibril structure context.

    If data/structure_features.csv exists (columns: position, sasa, ss_sheet),
    produced by DSSP on a fibril PDB (6CU7/6A6B), it is joined by position.
    Otherwise a coarse, explicitly-flagged proxy is used: the NAC domain is the
    fibril core, which is true at the domain level but not residue-resolved.
    """
    path = Path(path)
    out = pd.DataFrame(index=df.index)
    if path.exists():
        ref = pd.read_csv(path).set_index("position")
        out["sasa"] = df.position.map(ref["sasa"]).astype(float)
        out["ss_sheet"] = df.position.map(ref["ss_sheet"]).astype(float)
    else:
        out["sasa"] = np.nan  # unknown without DSSP; imputed inside the fold
        out["ss_sheet"] = np.nan
    out["in_fibril_core"] = [int(domain_of(p) == "nac") for p in df.position]
    return out


def conditions_block(df: pd.DataFrame) -> pd.DataFrame:
    return df[CONDITION_COLS].copy()


def build_features(df: pd.DataFrame):
    """Return (X, groups, y, blocks).

    X         : all feature columns (raw; encode/scale inside the fold)
    groups    : the study key for leave-one-study-out splitting
    y         : target = log10(t_half_h)
    blocks    : {"sequence": [...], "structure": [...], "conditions": [...]}
    """
    seq = sequence_block(df)
    struct = structure_block(df)
    cond = conditions_block(df)
    X = pd.concat([seq, struct, cond], axis=1)
    y = np.log10(df["t_half_h"].astype(float))
    groups = df["study"].astype(str)
    blocks = {
        "sequence": list(seq.columns),
        "structure": list(struct.columns),
        "conditions": list(cond.columns),
    }
    return X, groups, y, blocks
