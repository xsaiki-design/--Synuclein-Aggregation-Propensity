"""One-command rerun of the whole analysis: `python -m asyn.pipeline`.

Loads the curated dataset, runs the variance decomposition and the validated
predictor, writes the figures and a results summary. Refuses silently nothing:
if the dataset is the synthetic stand-in, every output says so.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from . import decompose, predict, report
from .data import load_curated

RESULTS_MD = Path("reports/results.md")


def _is_synthetic(df: pd.DataFrame) -> bool:
    return "synthetic" in df.columns and bool(df["synthetic"].fillna(0).eq(1).all())


def main() -> None:
    df = load_curated()
    synthetic = _is_synthetic(df)

    dec = decompose.run(df)
    pred = predict.run(df)

    report.fig_between_lab(df)
    report.fig_decomposition(dec["nested_cv_r2"], dec["study_icc"])
    report.fig_model_vs_baselines(pred["metrics"])

    RESULTS_MD.parent.mkdir(parents=True, exist_ok=True)
    with open(RESULTS_MD, "w") as fh:
        fh.write("# Results\n\n")
        if synthetic:
            fh.write("> **SYNTHETIC STAND-IN DATA — not a scientific result.** "
                     "Replace `data/aggregation_curated.csv` with the real curated "
                     "set (`make data`) and rerun.\n\n")
        icc = dec["study_icc"]
        fh.write("## Variance decomposition\n\n")
        fh.write(f"Study-level ICC (share of variance at the lab/condition level): "
                 f"**{icc['icc_study']:.2f}** "
                 f"({icc['n_studies']} studies, n={icc['n_obs']}).\n\n")
        fh.write(dec["nested_cv_r2"].to_markdown(index=False) + "\n\n")
        fh.write("## Prediction vs. baselines\n\n")
        fh.write(f"Target: {pred['target']}. "
                 "Leave-one-study-out, nested tuning, bootstrap 95% CI over studies.\n\n")
        fh.write(pred["metrics"].to_markdown(index=False) + "\n\n")
        if pred["note"]:
            fh.write(f"_{pred['note']}_\n\n")
        fh.write("## Figures\n\n")
        for f in ["fig1_between_lab.png", "fig2_decomposition.png",
                  "fig3_model_vs_baselines.png"]:
            fh.write(f"- `reports/figures/{f}`\n")

    tag = "[SYNTHETIC]" if synthetic else "[REAL]"
    print(f"{tag} wrote {RESULTS_MD} and 3 figures under reports/figures/")
    print(dec["nested_cv_r2"].to_string(index=False))
    print(pred["metrics"].to_string(index=False))


if __name__ == "__main__":
    main()
