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


def fig_variance(perm) -> Path:
    """Observed variance explained, each factor against its own chance level.

    The chance marks matter: with 5 degrees of freedom on ~16 cells a random
    mutant factor already reaches ~0.4, so the bars alone would overstate the
    contrast.
    """
    d = perm[perm.dimension != "ALL (pooled)"].copy()
    x = np.arange(len(d))
    fig, ax = plt.subplots(figsize=(7.8, 4.6))
    ax.bar(x - 0.19, d["r2_mutant"], 0.38, label="sequence (which mutant)", color=SEQ_C)
    ax.bar(x + 0.19, d["r2_category"], 0.38, label="conditions (which category)",
           color=COND_C)
    for xi, (cm, cc) in enumerate(zip(d["chance_mutant"], d["chance_category"])):
        ax.hlines(cm, xi - 0.38, xi, color="k", lw=1.6, ls="--")
        ax.hlines(cc, xi, xi + 0.38, color="k", lw=1.6, ls="--")
    ax.plot([], [], "k--", lw=1.6, label="chance level (shuffled labels)")
    for xi, (a, b, pm) in enumerate(zip(d["r2_mutant"], d["r2_category"], d["p_mutant"])):
        ax.text(xi - 0.19, a + 0.02, f"{a:.2f}", ha="center", fontsize=9)
        ax.text(xi + 0.19, b + 0.02, f"{b:.2f}", ha="center", fontsize=9)
        ax.text(xi - 0.19, a - 0.07, f"p={pm:.3f}", ha="center", fontsize=7.5, color="white")
    ax.set_xticks(x)
    ax.set_xticklabels(d["dimension"])
    ax.set_ylabel("variance explained (support-weighted $R^2$)")
    ax.set_ylim(0, 1.08)
    ax.set_title("Mutant identity beats chance; condition category does not")
    ax.legend(frameon=False, fontsize=8.5, loc="upper right")
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


def run(n_perm: int = 1000) -> list[Path]:
    df = si_analysis.load()
    perm = si_analysis.permutation_test(df, n_perm=n_perm)
    Path("reports").mkdir(exist_ok=True)
    perm.to_csv("reports/permutation_test.csv", index=False)
    figs = [fig_variance(perm), fig_scores_by_mutant(df)]
    types = si_analysis.TYPES_CSV
    if types.exists():
        figs.append(fig_type_by_domain(si_analysis.type_by_domain(si_analysis.load_types())))
    return figs




def fig_type_by_domain(mat) -> Path:
    """16 substitution transitions x 3 regions. Red accelerates, blue slows —
    the colour convention the SI's own heat maps use."""
    fig, ax = plt.subplots(figsize=(6.4, 7.2))
    im = ax.imshow(mat.values, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(mat.shape[1]))
    ax.set_xticklabels(mat.columns, rotation=20, ha="right")
    ax.set_yticks(range(mat.shape[0]))
    ax.set_yticklabels(mat.index)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat.values[i, j]
            ax.text(j, i, "—" if np.isnan(v) else f"{v:.2f}", ha="center",
                    va="center", fontsize=8,
                    color="white" if (not np.isnan(v) and abs(v) > 0.55) else "#222")
    ax.set_title("Which substitutions change aggregation, and where\n"
                 "(blue = slows, red = accelerates; — = no data)", fontsize=11)
    fig.colorbar(im, ax=ax, shrink=0.6, label="aggregation score")
    return _save(fig, "si_fig3_type_by_domain.png")


if __name__ == "__main__":
    for p in run():
        print(f"wrote {p}")
