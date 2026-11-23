# What the real SI actually contains (and the first real result)

The project plan assumed a 131-row table of single mutants with continuous **T½**
(half-time to max ThT) and full per-row conditions. **That table is not in this
supplementary PDF.** Reading all 12 pages of Pancoe et al. 2022 SI:

- The outcome is a **categorical aggregation score** in {−1, 0, +1} (reduces /
  no effect / accelerates), **normalized to the WT control**. No continuous T½.
- Tables **S3–S6**: scores averaged by mutation *type* and domain.
- Tables **S7–S10**: scores for the 6 familial mutants (A30P, E46K, H50Q, G51D,
  A53E, A53T) × condition category (buffer / salt / shaking / αS concentration),
  each with a support count. **This is the only per-(variant × condition) data
  in the file**, so it is the sequence-vs-conditions dataset.

`data/conditions_familial.csv` (77 rows, 301 underlying measurements) is
extracted from S7–S10 by `asyn.ingest_si` — numbers decoded from the PDF's own
text layer (the font's fixed CID offset), labels read from the rendered tables,
row counts asserted, and 5 values checked against the rendered pages. Real,
verified, no fabrication.

## First real result (`asyn.si_analysis`)
Support-weighted variance explained in the aggregation score:

| dimension | mutant R² | *chance* | p | condition R² | *chance* |
|---|---|---|---|---|---|
| buffer | 0.94 | *0.38* | 0.002 | 0.05 | *0.16* |
| concentration | 0.83 | *0.39* | 0.005 | 0.07 | *0.17* |
| salt | 0.88 | *0.45* | 0.013 | 0.08 | *0.21* |
| shaking | 0.83 | *0.36* | 0.004 | 0.03 | *0.23* |

**What the denominator is.** These R² values are fitted on condition-category
*means*, not on individual experiments: each of the 14–19 points already
averages up to 16 measurements, so within-cell disagreement is gone before the
variance is computed. The same model at experiment level explains **at least
51–58%** (`asyn.si_analysis.variance_explained_bounds`), taking within-cell
variance at its maximum for labels in {−1,0,+1}. The true figure is somewhere
between that and the numbers below, and S7–S10 report no SDs, so it cannot be
narrowed from this source. Quote the range, not the ceiling.

Two further limits on how far these numbers travel: the scores are normalised
to a matched WT control within each study, so this is variance in *normalised
effect* rather than in reported aggregation; and all six familial mutants are
N-terminal (residues 30–53), so this is an N-terminus result.

**Raw R² alone would overstate this, and the chance columns are why.** Mutant
is a 5-degree-of-freedom factor and a condition category has 1–2, on only
14–19 cells, so the larger factor explains more variance mechanically: a
*randomly shuffled* mutant label already reaches R² ≈ 0.4 per dimension. Each
factor is therefore compared against its own permutation null
(`asyn.si_analysis.permutation_test`, 1000 shuffles; `reports/permutation_test.csv`).

One bug had to be fixed before this comparison was fair: `"None"` — a real
category meaning no added salt or no shaking — is in pandas' default NA list,
so those 12 cells were read as missing and silently dropped from the *category*
model while the *mutant* model kept them (58 rows against 66). Since the
quiescent/no-salt levels are where a condition effect should be largest, the
bug favoured the conclusion. Fixed; both models now use all 66 cells, and the
finding holds (salt's condition R² moves 0.14 → 0.08).

Against that baseline the result is cleaner than the raw numbers suggested:
mutant identity beats chance in every table (p = 0.002–0.013),
while **condition category lands *below* its own chance level everywhere**
(0.03–0.08 observed against 0.16–0.23 expected, all p ≥ 0.71). Condition
category carries no detectable signal at all — it does not merely lose to
sequence.

**Why there is no pooled row.** An earlier version reported a pooled result
(R² 0.86 against a 0.10 chance level, p = 0.001). It was invalid. Every
experiment has a buffer *and* a salt *and* a shaking speed *and* a
concentration, so the four tables are four partitions of the **same**
experiments — A30P's buffer cells total 23 measurements and its salt cells total
the same 23. Stacking the tables counted each experiment up to four times (301
summed support against at most 93 in any one table), which is what drove that
pooled chance level down and its p-value to the floor. The four tables are
therefore four consistent views of one dataset, not four independent
replications, and their p-values cannot be combined either.

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

**Conclusion (for this dataset): which mutation it is explains the aggregation
score; which condition category it was measured under does not** — the opposite
of the plan's hypothesis that conditions would swamp sequence. A legitimate,
reportable negative result.

## The sequence side (Tables S3–S6)

`asyn.ingest_types` extracts all four mutation-type tables (104 rows, 866
measurements) — same CID decoding, but **S5 and S6 use a different glyph
offset**, so the offset is calibrated per page and every value is bounds-checked
and spot-checked against the rendered pages (9 checks, in the ingester).

Support-weighted mean score by region — mutations in the fibril core slow
aggregation most, mutations in the N-terminus barely matter:

| region | weighted mean score | measurements |
|---|---|---|
| **NAC (61–95)** | **−0.57** | 145 |
| C-terminus | −0.35 | 110 |
| N-terminus | −0.10 | 178 |

Across the whole protein the largest effects (n ≥ 10) are **X → Pro (−0.88,
n=33)** and **X → PG (−0.85, n=34)** — proline blocks aggregation, the
β-sheet-breaking behaviour the fibril model predicts, recovered from the data.

**But those whole-protein rows aggregate across regions and should not be read
as regional results.** X → Pro's 33 measurements split 16 in the N-terminus and
17 in the NAC, with **none in the C-terminus** — so proline substitution is
simply untested there. Ranked within each region (support floor n ≥ 5,
`asyn.si_analysis.strongest_by_region`):

| region | top-ranked substitution | score | sd | n |
|---|---|---|---|---|
| N-terminus | nβ → β | +1.000 | 0.000 | 5 |
| NAC | **X → Pro** | **−0.941** | 0.235 | **17** |
| C-terminus | B → nB | −1.000 | 0.000 | 7 |

**Only the NAC row should be quoted as a result.** The other two are the
maximum of 5–7 ranked candidates, and the maximum of a small sample is biased
away from zero: both land on exactly ±1.000 with sd = 0 on n = 5 and n = 7.
They are order statistics, not effect estimates. `region_ranking()` prints
every candidate so the selection is visible.

The regional view carries information the aggregate loses: **the N-terminus's
strongest effect accelerates aggregation, while the NAC's and C-terminus's
block it.** Averaging over regions partly cancels those opposing signs.

Two cautions on reading that table. The N-terminus and C-terminus winners rest
on n=5 and n=7 with sd = 0 — unanimous, but small; X → Pro in the NAC is the
only one resting on substantial support. And in the N-terminus X → PG is
numerically identical to X → Pro (−0.813, n=16) because every PG substitution
there is a proline one — X → Gly has no N-terminal measurements. Per-region
rankings: `asyn.si_analysis.region_ranking`.

![type by domain](reports/figures/si_fig3_type_by_domain.png)

### The essential caveat
The scores are **already WT-normalized within each study**, which removes much
of the raw between-lab condition variance *by construction*. So this result
says conditions add little **once effects are expressed relative to a matched WT
control** — which is exactly the authors' recommendation. Testing whether raw,
un-normalized aggregation rates are condition-dominated needs the raw kinetics
(T½), which this SI does not provide. The continuous-T½ regression scaffold
(`asyn.features/splits/decompose/predict`) is kept ready for that data.


## Does the literature agree? (`asyn.discordance`)

Tables S3–S6 report a **standard deviation** beside every mean. Used everywhere
else as noise, it is really a record of how much the published experiments
disagreed. Because the scores are means of {−1, 0, +1} labels, a group's mean,
SD and count determine the vote exactly:

> n₊ = n(s² + m² + m)/2  ·  n₋ = n(s² + m² − m)/2  ·  n₀ = n − n₊ − n₋

So each row can be read back as *"k experiments found acceleration, j found
reduction"* — no modelling, just algebra.

**This is also a third, independent check on the extraction.** The
reconstructed counts come out as whole experiments (max integrality error
**0.026**), and no row breaches the bound sd ≤ √(1−mean²). A mis-decoded mean
or SD would fail both.

### Proline is the one thing nobody disputes

| group | accelerate | no effect | reduce |
|---|---|---|---|
| X → Pro, whole protein | **0** | 4 | **29** |
| X → Pro, NAC | **0** | 1 | **16** |
| X → Pro, N-terminus | **0** | 3 | **13** |
| X → PG, NAC | **0** | 2 | **16** |

**No published experiment in this compilation reports a proline substitution
accelerating aggregation.** That is a stronger statement than the −0.88 mean.

**It is a count of published reports, not of experiments performed.** A result
matching a strong prior expectation — proline breaks β-sheet geometry — is
easier to publish than a null or contradictory one, so 0 of 33 is consistent
with a real effect, with selective publication, or with both. This data cannot
separate them, and no correction is possible without the unpublished record.

### What the field is deadlocked about

| group | accelerate | no effect | reduce | split |
|---|---|---|---|---|
| nB → B, whole protein | 5 | 2 | 5 | **0.50** |
| nP → P, N-terminus | 5 | 0 | 4 | 0.44 |
| B → nB, N-terminus | 3 | 0 | 4 | 0.43 |
| Sm → Lg, whole protein | 5 | 4 | 7 | 0.42 |

**Changes in bulk, hydrophobicity and size split the literature almost evenly.**
These are exactly the substitutions whose effect is most likely to depend on
something unrecorded — protein prep, surface, assay geometry — because the
sequence change alone evidently does not determine the outcome.

![vote split](reports/figures/si_fig4_vote_split.png)

### A regional asymmetry

Measured on the recovered votes (minority share of directional reports,
support-weighted), **regions only** — the whole-protein rows of Table S3 are
the exact sum of the three regional tables for all 26 types, so including them
would double count:

| region | mean split | unanimous types | deadlocked types |
|---|---|---|---|
| N-terminus | 0.170 | 3 of 7 | 2 |
| NAC | 0.046 | 3 of 6 | 0 |
| C-terminus | 0.000 | 5 of 5 | 0 |

The published record is divided mainly about the **N-terminus**; the NAC and
C-terminus are close to undivided.

**This corrects an earlier version of this section**, which used a dispersion
ratio (sd ÷ the maximum attainable at that mean) and reported N-terminus 0.76,
NAC 0.74 and C-terminus 0.34, describing the NAC as "3 contested, 0
consensus". That was wrong: the NAC vote is 3 unanimous and 0 deadlocked. The
ratio systematically misreads groups with extreme means — the same failure that
labelled X → Pro in the NAC "mixed" when its vote is 16-0-1. The ratio has been
removed from the codebase.



## A second source, and why it is not corroboration (`asyn.cross_study`)

[Flagmeier et al. 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC5027465) covers
the same familial mutants but reports orders of magnitude relative to WT,
separated by mechanism. Extracted by `asyn.ingest_flagmeier`; every row carries
the verbatim sentence it came from, because the paper states its effects in
prose and tabulates no rate constants.

**Pancoe cites Flagmeier.** Pancoe's deposited reference list (240 entries,
retrieved by DOI content negotiation, cached in `data/pancoe_references.csv`,
queryable with `cross_study.cited_by_pancoe`) contains it at
`10.1073/pnas.1604645113`, along with Buell 2014, Stephens 2020 and Newberry
2020. A citation does not prove a paper's measurements were tabulated among the
375 examples, but it removes the presumption of independence. **This comparison
is therefore not an external check**, and is reported as case observations.

Of the candidate sources checked, only **Ulamec et al. 2022**
(`10.1038/s41467-022-32687-1`) is *not* cited. It is also the held-out test set
the original plan named, so it is the next source to add.

### Case observations, not a rate
Against Flagmeier's lipid-induced nucleation, A53T (+0.91 vs +2) and G51D
(−0.50 vs −2) agree in direction; **E46K (+0.70 vs −2) and H50Q (+0.85 vs −1)
flip**. At n=4 that split is binomial p = 1.0 — not distinguishable from chance,
and not quotable as an agreement rate.

One point cuts against the circularity concern and is recorded rather than
suppressed: had Flagmeier's labels been folded into Pancoe's scores, the two
would tend to agree, and for E46K and H50Q they disagree sharply.

### The one result needing no comparison
Within Flagmeier alone — same five mutants, one lab, one protein prep — the
spread across those mutants differs by mechanism:

| mechanism | spread | fold range |
|---|---|---|
| lipid-induced nucleation | 4.0 log10 | **10,000×** |
| secondary nucleation | 2.0 log10 (censored) | 100× |
| elongation | 1.0 log10 | 10× |

**How much a mutation matters depends on which mechanism is measured**, by three
orders of magnitude. No cross-study comparison is involved, so none of the
independence or endpoint caveats apply. The values are prose-stated
approximations at n=5, so treat them as directional.

## Status of every claim

| claim | status |
|---|---|
| No *published* experiment reports X → Pro accelerating (0/29) | **holds** — publication bias unmeasurable |
| Bulk/hydrophobicity deadlocked (nB → B 5–5) | **holds** |
| NAC is the most affected region | **holds** — −0.70 transitions-only, −0.57 with controls, at every threshold 1–10 |
| Mechanism dependence (10,000× vs 10×) | **holds** — within-study, n=5, prose-stated |
| Condition *category* carries no signal | **holds** — below its own chance baseline in all four tables |
| Mutant identity explains 83–94% | **qualified** — ceiling; experiment-level floor 51–58%; N-terminus only; normalised effect |
| X → Pro strongest substitution | **qualified** — only the NAC row (n=17); the −0.88/n=33 figure is a cross-region aggregate with no C-terminal data |
| Pooled R² = 0.86, p = 0.001 | **withdrawn** — pseudo-replication across the four condition tables |
| Regional "contested" labels (0.76 / 0.74 / 0.34) | **withdrawn** — dispersion ratio misread extreme means |
| Three-way region ordering (−0.10 / −0.35 / −0.57) | **withdrawn** — N- and C-terminus indistinguishable without control rows |
| Per-region strongest for N- and C-terminus (±1.000) | **withdrawn** — selection on the maximum at n=5 and 7 |
| Flagmeier as external corroboration | **withdrawn** — Pancoe cites it |

Every finding rests on one meta-analysis, and no independent external check has
yet been performed.
