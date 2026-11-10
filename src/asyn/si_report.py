"""Figures for the REAL result (Tables S7–S10).

Two figures carry the finding:
  1. variance explained per condition dimension — mutant identity (sequence)
     against condition category;
  2. every extracted cell plotted by mutant — the visual form of the same
     claim: points for one mutant sit close together across very different
     conditions, while different mutants sit far apart.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from . import si_analysis  # noqa: E402

FIG_DIR = Path("reports/figures")
SEQ_C, COND_C = "#205493", "#d96d2b"


def fig_variance(ve) -> Path:
    d = ve[ve.dimension != "ALL (pooled)"].copy()
    x = np.arange(len(d))
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.bar(x - 0.19, d["r2_sequence_mutant"], 0.38,
           label="sequence (which mutant)", color=SEQ_C)
    ax.bar(x + 0.19, d["r2_conditions_category"], 0.38,
           label="conditions (which category)", color=COND_C)
    ax.set_xticks(x)
    ax.set_xticklabels(d["dimension"])
    ax.set_ylabel("variance explained (support-weighted $R^2$)")
    ax.set_ylim(0, 1)
    ax.set_title("Aggregation score: sequence explains it, conditions barely move it")
    ax.legend(frameon=False)
    for xi, (a, b) in enumerate(zip(d["r2_sequence_mutant"], d["r2_conditions_category"])):
        ax.text(xi - 0.19, a + 0.02, f"{a:.2f}", ha="center", fontsize=9)
        ax.text(xi + 0.19, b + 0.02, f"{b:.2f}", ha="center", fontsize=9)
    return _save(fig, "si_fig1_variance.png")


def fig_scores_by_mutant(df) -> Path:
    mutants = ["A30P", "E46K", "H50Q", "G51D", "A53E", "A53T"]
    marks = {"buffer": "o", "salt": "s", "shaking": "^", "concentration": "D"}
    fig, ax = plt.subplots(figsize=(7.6, 4.4))
    rng = np.random.default_rng(0)
    for i, m in enumerate(mutants):
        sub = df[df.mutant == m]
        for dim, mk in marks.items():
            s = sub[sub.dimension == dim]
            if s.empty:
                continue
            ax.scatter(i + rng.uniform(-0.17, 0.17, len(s)), s["score"],
                       marker=mk, s=28, alpha=0.75, color=SEQ_C,
                       label=dim if i == 0 else None)
        mu = np.average(sub["score"], weights=sub["support"])
        ax.hlines(mu, i - 0.3, i + 0.3, color=COND_C, lw=2.5, zorder=5)
    ax.set_xticks(range(len(mutants)))
    ax.set_xticklabels(mutants)
    ax.set_ylabel("aggregation score  (−1 slows · +1 accelerates)")
    ax.axhline(0, color="k", lw=0.8, ls=":")
    ax.set_title("Each mutant's score holds across buffer, salt, shaking and concentration")
    ax.legend(frameon=False, fontsize=8, ncol=4, loc="lower center")
    ax.set_ylim(-1.35, 1.35)
    return _save(fig, "si_fig2_scores_by_mutant.png")


def _save(fig, name) -> Path:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    out = FIG_DIR / name
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def run() -> list[Path]:
    res = si_analysis.run()
    df = si_analysis.load()
    return [fig_variance(res["variance_explained"]), fig_scores_by_mutant(df)]


if __name__ == "__main__":
    for p in run():
        print(f"wrote {p}")
