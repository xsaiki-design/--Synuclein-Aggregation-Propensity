# α-Synuclein aggregation: sequence vs. conditions

**Question:** how much of the variation in reported α-synuclein aggregation
behaviour is explained by the protein's **sequence**, and how much by the
**experimental conditions** each measurement was made under?

Data: the open-access supplementary tables of **Pancoe et al. 2022**
(*J Mol Biol* 434(23):167859, [PMC9922159](https://pmc.ncbi.nlm.nih.gov/articles/PMC9922159)),
whose authors explicitly tabulated these data for computational reuse.
Everything below is extracted from that PDF and verified against it — see
[`SI_FINDINGS.md`](SI_FINDINGS.md) for the full write-up and caveats.

## Results

**1. Sequence explains the score; condition category carries no signal.**
Across the 7 familial mutants, *which mutant* explains ~86% of the variance in
aggregation score against a shuffled-label chance level of ~10% (p = 0.001);
*which condition category* (buffer / salt / shaking / concentration) explains
~4% against a chance level of ~18% — i.e. **below chance**. The chance levels
matter: mutant is a 5-df factor on 14–19 cells, so raw R² alone would overstate
the contrast. Every factor is tested against its own permutation null.

![variance explained](reports/figures/si_fig1_variance.png)

Each mutant's score holds across very different conditions, while the mutants
themselves sit far apart — though the real exceptions are visible and kept
(A30P in low salt flips to +0.67; G51D at medium shaking to +0.50):

![scores by mutant](reports/figures/si_fig2_scores_by_mutant.png)

**2. The NAC core is where sequence matters.** Support-weighted mean score by
region: **NAC (61–95) −0.57**, C-terminus −0.35, N-terminus −0.10. The
strongest single substitution is **X → Pro (−0.88, n=33)**, strongest inside
the NAC (−0.94) — the β-sheet-breaking effect the fibril model predicts,
recovered from the data.

![type by domain](reports/figures/si_fig3_type_by_domain.png)

> **The caveat that matters:** the scores are *already WT-normalized within
> each study*, which removes much raw between-lab condition variance by
> construction. So result 1 says conditions add little **once effects are
> expressed against a matched WT control** — which is exactly the authors'
> recommendation, not evidence that raw rates are condition-independent.

## The data, and how it was verified
The outcome in this SI is a **categorical aggregation score** in {−1, 0, +1}
(reduces / no effect / accelerates), **not** continuous T½. T½ is not tabulated
in this file.

| file | source | rows |
|---|---|---|
| `data/conditions_familial.csv` | Tables S7–S10 (mutant × condition) | 77 |
| `data/type_scores.csv` | Tables S3–S6 (mutation type × region) | 104 |

The PDF's label column is a CID-mapped font with no Unicode table, and **the
numeric glyph offset differs between pages** (S5/S6 vs S7–S10) — a fixed map
silently decodes `1.000` as `20111`. So the ingesters calibrate the offset per
page, read columns positionally, bounds-check every value, and assert 9 values
read off the rendered pages. Drift fails loudly rather than producing plausible
garbage. No value is transcribed by hand into the data.

```bash
make setup                      # pip install -e .[dev]
python -m asyn.ingest_si        # Tables S7–S10  -> conditions_familial.csv
python -m asyn.ingest_types     # Tables S3–S6   -> type_scores.csv
python -m asyn.si_analysis      # the decomposition
python -m asyn.si_report        # the three figures
make test                       # 18 tests
```
`data/raw/SI.pdf` is gitignored; download it from the PMC link above if you
need to re-run extraction.

## Also in here: a validated-prediction scaffold
`asyn.features / splits / decompose / predict` implement the full
leave-one-study-out regression pipeline from the original plan — study-grouped
splits, fold-local preprocessing, nested CV, three baselines, bootstrap CIs, no
neural network. **It needs raw kinetics (T½), which this SI does not contain**,
so it currently runs on a labelled synthetic stand-in (`make all`) and produces
no scientific result. It is kept ready for that data.

## Layout
```
src/asyn/   ingest_si · ingest_types · si_analysis · si_report   (real)
            features · splits · decompose · predict · synth      (T½ scaffold)
data/       conditions_familial.csv · type_scores.csv · data_dictionary.md
reports/figures/  si_fig1..3 (real)
SI_FINDINGS.md · preregistration.md
```
