# Data dictionary

Source for everything real in this repo: the open-access supplementary tables of
**Pancoe et al. 2022** (*J Mol Biol* 434(23):167859,
[PMC9922159](https://pmc.ncbi.nlm.nih.gov/articles/PMC9922159)), whose authors
tabulated these data for computational reuse.

**The outcome in this source is a categorical aggregation score**, not a
continuous rate:

| value | meaning |
|---|---|
| `-1` | the change *reduces* aggregation |
| `0` | no effect |
| `+1` | the change *accelerates* aggregation |

Scores are **normalised to a matched WT control in the same study**, and the
tabulated figures are *averages* of that score over a group of experiments. T½
(time to half-maximal ThT fluorescence) is **not tabulated in this SI** — see
"The T½ scaffold" at the bottom.

---

## `conditions_familial.csv` — Tables S7–S10 (77 rows)
One row per (familial mutant × condition category). Built by `asyn.ingest_si`.

| Column | Type | Meaning |
|---|---|---|
| `mutant` | str | One of the 6 familial mutants: A30P, E46K, H50Q, G51D, A53E, A53T. |
| `dimension` | str | Which condition table the row came from: `buffer`, `salt`, `shaking`, `concentration`. |
| `category` | str | The level within that dimension — see below. |
| `score` | float | Mean aggregation score for that cell, in [−1, +1]. |
| `support` | int | Number of underlying experiments averaged into `score`. |

Categories by dimension:

| dimension | levels |
|---|---|
| `buffer` | `PBS`, `Tris`, `Other` |
| `salt` | `High`, `Low`, `None` |
| `shaking` | `Fast`, `Medium`, `Slow`, `None` |
| `concentration` | `High`, `Normal`, `Low` |

> **`"None"` is a real level**, meaning *no added salt* or *no shaking* — not a
> missing value. It is in pandas' default NA list, so it must be read with
> `keep_default_na=False` or those 12 rows silently vanish. `asyn.si_analysis.load`
> does this; a regression test guards it.

> **The four dimensions are not independent samples.** Every experiment has a
> buffer *and* a salt *and* a shaking speed *and* a concentration, so these are
> four partitions of the *same* experiments. A30P's buffer cells total 23
> measurements and its salt cells total the same 23. Never pool the four tables
> — that counts each experiment up to four times. Analyse per dimension.

`support = 0` means the SI printed no measurements for that cell; such rows are
filtered out before analysis rather than treated as a score of zero.

---

## `type_scores.csv` — Tables S3–S6 (104 rows)
One row per (sequence region × mutation type). Built by `asyn.ingest_types`.

| Column | Type | Meaning |
|---|---|---|
| `scope` | str | `whole_protein` (S3), `N-terminus` (S4), `NAC` (S5), `C-terminus` (S6). |
| `type` | str | Mutation type, e.g. `X -> Pro`, `Sm -> Lg`, `nβ -> β`, or a same-class `… Cntrl` row. |
| `average` | float | Mean aggregation score in [−1, +1]; empty when `support` is 0. |
| `sd` | float | Standard deviation of that mean; empty when `support` is 0. |
| `support` | int | Number of underlying experiments. |

Region boundaries follow the SI: N-terminus 1–60, NAC 61–95, C-terminus 96–140.
The 26 `type` values are the same list, in the same order, in all four tables —
16 substitution transitions plus 10 same-class control rows.

### Extraction integrity
Both CSVs are decoded from the PDF's own text layer, not transcribed by hand.
The label column uses a CID-mapped font with no Unicode table, and **the numeric
glyph offset differs between pages** (S5/S6 vs S7–S10) — a fixed map decodes
`1.000` as `20111`. So the ingesters:

1. calibrate the glyph offset **per page**, rejecting any offset that yields a
   value outside known bounds (score ∈ [−1,1], sd ∈ [0,1], support a
   non-negative integer);
2. read the columns **positionally**, because the label column leaks digits;
3. assert values read off the rendered pages — **5** spot checks in
   `ingest_si`, **9** in `ingest_types` — on every build.

A future font or layout change fails loudly instead of producing plausible
garbage.

---

## `cross_study_effects.csv` — effects from other papers (15 rows)
One row per (study × mutant × mechanism). Built by `asyn.ingest_flagmeier`;
one ingester module per paper. Full method and limits:
[`COMPARISON_METHOD.md`](../COMPARISON_METHOD.md).

| Column | Meaning |
|---|---|
| `study`, `year`, `source` | provenance. **`study` is a grouping key — never pool across it.** |
| `mutant` | join key to the Pancoe data |
| `endpoint` | what was actually measured (`lipid_nucleation`, `secondary_nucleation`, `elongation`) — never compared across without a declared mapping |
| `condition_class` | `lipid_vesicles` / `seeded_acidic` / `seeded` / `bulk` / `cellular` |
| `normalized_to_wt` | whether the value is relative to a WT control |
| `value`, `value_type` | e.g. `2.0`, `log10_fold_vs_wt` |
| `censored` | true where the source gives a bound ("undetectably slow"), not a value |
| `extraction` | `text_statement` / `table_decoded` / `figure_digitised` |
| `quote` | verbatim sentence supporting the value, so every number is auditable |

> **Flagmeier values are prose-stated approximations**, not digitised
> measurements — the paper gives orders of magnitude in text and tabulates no
> rate constants. Good for direction and rough magnitude; not for regression on
> exact values. The two `censored` rows report "undetectably slow"; their
> −2.0 is a floor standing in for a bound.

---

## The T½ scaffold — `aggregation_curated.csv` (not real data)
`asyn.features / splits / decompose / predict` implement a leave-one-study-out
regression pipeline whose outcome is continuous **T½ in hours**, one row per
individual experiment with its own conditions. **No such table exists in this
SI**, so that pipeline runs only on a labelled synthetic stand-in
(`asyn.synth`, `synthetic=1`, `make data-synth`) and produces no scientific
result. It is kept ready for raw kinetics if that data is obtained.

Its schema: `variant`, `study`, `position`, `wt_residue`, `mut_residue`,
`t_half_h`, `t_half_rel_wt`, plus per-row conditions (`pH`, `temperature_C`,
`protein_conc_uM`, `agitation`, `buffer`, `seeds`, `lipids`, `ThT_conc_uM`).
`study` is the grouping key for leave-one-study-out validation.

If that data is ever obtained, note that T½ conflates nucleation and elongation
into one number; Ohgita et al. 2022 ([PMC9043213](https://pmc.ncbi.nlm.nih.gov/articles/PMC9043213))
separate those phases. **Missing conditions are data** — a condition a paper
does not report stays NaN rather than being imputed at curation time, and the
frequency of missingness is reported.

## Other provenance
- Canonical sequence for feature construction: UniProt **P37840**, 140 aa.
- Fibril structure features (optional block): PDB **6CU7 / 6A6B** via DSSP →
  `data/structure_features.csv` (`position, sasa, ss_sheet`). Without it, the
  structure block carries only `in_fibril_core` (NAC = residues 61–95).
