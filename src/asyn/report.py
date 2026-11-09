"""Figures and the results summary (Phases 2, 4, 5, 6 outputs).

The four figures that matter, per the plan:
  1. between-lab variability (same variant, different labs);
  2. the variance decomposition (sequence vs structure vs conditions);
  3. model vs. the three baselines, with bootstrap intervals;
  4. external held-out test (when that dataset is present).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

FIG_DIR = Path("reports/figures")


def _ax(title, xlabel=None, ylabel=None):
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.set_title(title)
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    return fig, ax


def fig_between_lab(df: pd.DataFrame) -> Path:
    """Spread in reported rate for variants measured in more than one study."""
    rep = df.groupby("variant").filter(lambda g: g["study"].nunique() > 1)
    fig, ax = _ax("Between-lab variability (same variant, different labs)",
                  "variant", "log10 t_half (h)")
    if rep.empty:
        ax.text(0.5, 0.5, "no variant measured in >1 study", ha="center", va="center")
    else:
        rep = rep.assign(logt=np.log10(rep["t_half_h"]))
        order = rep.groupby("variant")["logt"].std().sort_values(ascending=False).index
        data = [rep.loc[rep.variant == v, "logt"].values for v in order]
        ax.boxplot(data, showfliers=False)
        ax.set_xticks(range(1, len(order) + 1))
        ax.set_xticklabels(list(order), rotation=45, ha="right")
        ax.scatter(np.concatenate([[i + 1] * len(d) for i, d in enumerate(data)]),
                   np.concatenate(data), s=14, alpha=0.6, color="#205493")
    out = _save(fig, "fig1_between_lab.png")
    return out


def fig_decomposition(table: pd.DataFrame, icc: dict) -> Path:
    fig, ax = _ax("Variance explained out of sample (leave-one-study-out)",
                  None, "CV R^2")
    ax.bar(table["model"], table["cv_r2_mean"],
           yerr=table["cv_r2_std"], capsize=4, color="#205493")
    ax.axhline(0, color="k", lw=0.8)
    ax.annotate(f"study-level ICC = {icc['icc_study']:.2f}\n"
                f"({icc['n_studies']} studies, n={icc['n_obs']})",
                xy=(0.98, 0.95), xycoords="axes fraction", ha="right", va="top",
                fontsize=9, bbox=dict(boxstyle="round", fc="#eef3f8"))
    plt.xticks(rotation=20, ha="right")
    return _save(fig, "fig2_decomposition.png")


def fig_model_vs_baselines(metrics: pd.DataFrame) -> Path:
    fig, ax = _ax("Prediction vs. baselines (bootstrap 95% CI over studies)",
                  None, "R^2")
    m = metrics.copy()
    yerr = np.vstack([m["r2"] - m["r2_lo"], m["r2_hi"] - m["r2"]])
    ax.bar(m["model"], m["r2"], yerr=yerr, capsize=4, color="#205493")
    ax.axhline(0, color="k", lw=0.8)
    plt.xticks(rotation=20, ha="right")
    return _save(fig, "fig3_model_vs_baselines.png")


def fig_external(ext: pd.DataFrame | None) -> Path | None:
    if ext is None or ext.empty:
        return None
    fig, ax = _ax("External held-out test (Ulamec P1 scan)",
                  "observed log10 t_half", "predicted")
    ax.scatter(ext["observed"], ext["predicted"], color="#205493")
    lim = [min(ext.observed.min(), ext.predicted.min()),
           max(ext.observed.max(), ext.predicted.max())]
    ax.plot(lim, lim, "k--", lw=0.8)
    return _save(fig, "fig4_external.png")


def _save(fig, name) -> Path:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    out = FIG_DIR / name
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out
