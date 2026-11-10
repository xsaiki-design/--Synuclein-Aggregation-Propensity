# What the real SI actually contains (and the first real result)

The project plan assumed a 131-row table of single mutants with continuous **T½**
(half-time to max ThT) and full per-row conditions. **That table is not in this
supplementary PDF.** Reading all 12 pages of Pancoe et al. 2022 SI:

- The outcome is a **categorical aggregation score** in {−1, 0, +1} (reduces /
  no effect / accelerates), **normalized to the WT control**. No continuous T½.
- Tables **S3–S6**: scores averaged by mutation *type* and domain.
- Tables **S7–S10**: scores for the 7 familial mutants (A30P, E46K, H50Q, G51D,
  A53E, A53T) × condition category (buffer / salt / shaking / αS concentration),
  each with a support count. **This is the only per-(variant × condition) data
  in the file**, so it is the sequence-vs-conditions dataset.

`data/conditions_familial.csv` (77 rows, 301 underlying measurements) is
extracted from S7–S10 by `asyn.ingest_si` — numbers decoded from the PDF's own
text layer (the font's fixed CID offset), labels read from the rendered tables,
row counts asserted. Real, verified, no fabrication.

## First real result (`asyn.si_analysis`)
Support-weighted variance explained in the aggregation score:

| dimension | mutant (sequence) R² | condition R² |
|---|---|---|
| buffer | 0.94 | 0.05 |
| concentration | 0.83 | 0.07 |
| salt | 0.88 | 0.14 |
| shaking | 0.84 | 0.03 |
| **pooled** | **0.86** | **0.04** |

![variance explained](reports/figures/si_fig1_variance.png)

Within-mutant spread across conditions is **0.20–0.38×** the between-mutant
spread. Plotting every extracted cell makes the same point directly — each
mutant's score holds across very different buffers, salts, shaking speeds and
concentrations, while the mutants themselves sit far apart:

![scores by mutant](reports/figures/si_fig2_scores_by_mutant.png)

The exceptions are visible and worth keeping: A30P in low salt reads +0.67
against its −0.57 mean, and G51D at medium shaking reads +0.50 against −0.49.
Conditions can flip an individual reading even where they explain little
variance overall. The authors' own example: A30P scores −0.75 in PBS vs −0.31 in Tris,
but averaged PBS→Tris is −0.12 ± 0.31 — no generalizable buffer effect.

**Conclusion (for this dataset): sequence dominates conditions** — the opposite
of the plan's hypothesis that conditions would swamp sequence. A legitimate,
reportable negative result.

### The essential caveat
The scores are **already WT-normalized within each study**, which removes much
of the raw between-lab condition variance *by construction*. So this result
says conditions add little **once effects are expressed relative to a matched WT
control** — which is exactly the authors' recommendation. Testing whether raw,
un-normalized aggregation rates are condition-dominated needs the raw kinetics
(T½), which this SI does not provide. The continuous-T½ regression scaffold
(`asyn.features/splits/decompose/predict`) is kept ready for that data.
