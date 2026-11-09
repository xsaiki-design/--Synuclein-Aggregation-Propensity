"""Synthetic stand-in dataset — NOT REAL DATA.

Every value here is generated from a random seed. It exists only to exercise
the pipeline end-to-end before the real curated CSV (from the Pancoe et al.
2022 supplementary PDF) is in place. No scientific conclusion may be drawn
from it. The real data replaces it via `make data`.

It is written with a `synthetic=1` column and a header comment so it can never
be mistaken for the real thing downstream.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from .data import ASYN_SEQUENCE, REQUIRED_COLUMNS

_AA = "ACDEFGHIKLMNPQRSTVWY"
_HYDRO = {  # Kyte-Doolittle, only used to give synthetic rows some structure
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}


def generate(n_studies: int = 16, seed: int = 0) -> pd.DataFrame:
    """Generate a synthetic dataset of ~130 single-mutant rows across studies."""
    rng = np.random.default_rng(seed)
    rows = []
    for s in range(n_studies):
        study = f"synthlab_{s:02d}_2022"
        # Each study has its own fixed conditions (the whole point: conditions
        # are shared within a paper and drive a study-level offset).
        cond = dict(
            pH=round(float(rng.uniform(6.5, 7.6)), 2),
            temperature_C=float(rng.choice([25, 30, 37])),
            protein_conc_uM=float(rng.choice([20, 50, 70, 100])),
            agitation=str(rng.choice(["quiescent", "shaking", "beads"])),
            buffer=str(rng.choice(["PBS", "Tris", "phosphate"])),
            seeds=int(rng.choice([0, 1], p=[0.8, 0.2])),
            lipids=int(rng.choice([0, 1], p=[0.9, 0.1])),
            ThT_conc_uM=float(rng.choice([10, 20, 25, 50])),
        )
        study_offset = rng.normal(0, 0.6)  # lab-level shift in log t_half
        # WT control for this study.
        wt_log = 2.0 + study_offset + rng.normal(0, 0.1)
        wt_t = float(np.exp(wt_log))
        rows.append(_row(study, "WT", np.nan, np.nan, np.nan, wt_t, 1.0, cond, rng))
        # A handful of single mutants per study.
        for _ in range(int(rng.integers(5, 11))):
            pos = int(rng.integers(1, 141))
            wt_res = ASYN_SEQUENCE[pos - 1]
            mut_res = str(rng.choice([a for a in _AA if a != wt_res]))
            seq_effect = -0.15 * (_HYDRO[mut_res] - _HYDRO[wt_res]) / 9.0
            log_t = wt_log + seq_effect + rng.normal(0, 0.2)
            t = float(np.exp(log_t))
            rows.append(
                _row(study, f"{wt_res}{pos}{mut_res}", pos, wt_res, mut_res,
                     t, t / wt_t, cond, rng)
            )
    df = pd.DataFrame(rows)[REQUIRED_COLUMNS + ["synthetic"]]
    return df


def _row(study, variant, position, wt_res, mut_res, t_half, t_rel, cond, rng):
    return dict(
        variant=variant, study=study, position=position,
        wt_residue=wt_res, mut_residue=mut_res,
        t_half_h=round(t_half, 3), t_half_rel_wt=round(t_rel, 3),
        synthetic=1, **cond,
    )


def write(path: Path | str = "data/aggregation_curated.csv", seed: int = 0) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df = generate(seed=seed)
    with open(path, "w") as fh:
        fh.write("# SYNTHETIC STAND-IN DATA - NOT REAL. Replace via `make data`.\n")
        df.to_csv(fh, index=False)
    return path


if __name__ == "__main__":
    p = write()
    print(f"wrote {p} ({len(pd.read_csv(p, comment='#'))} rows) [SYNTHETIC]")
