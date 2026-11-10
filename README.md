# α-Synuclein aggregation: sequence vs. conditions

**Question:** how much of the variation in reported α-synuclein aggregation
rate (T½, time to half-maximal ThT fluorescence) is explained by the protein's
**sequence**, and how much by the **experimental conditions** each measurement
was made under?

This is deliberately *not* framed as "build a predictor." The literature
disagrees with itself — A30P has been reported to slow, accelerate, and not
affect aggregation — and that disagreement is the finding. The main result is a
**variance decomposition**, with a properly-validated predictor as support.

> **Status:** real data extracted. The SI PDF turned out to hold **categorical
> aggregation scores**, not continuous T½ — see **[`SI_FINDINGS.md`](SI_FINDINGS.md)**
> for what it actually contains and the first real result (sequence explains
> ~86% of score variance, conditions ~4%, pooled; with an important
> normalization caveat). `asyn.ingest_si` builds the real
> `data/conditions_familial.csv`; `asyn.si_analysis` runs the decomposition.
> The continuous-T½ regression modules below are kept ready for raw kinetics,
> which this SI does not provide; they currently run on a synthetic stand-in.

## Get the real data (one manual step)
The curated set comes from the open-access supplementary PDF of Pancoe et al.
2022. Its download is behind a browser proof-of-work gate, so fetch it by hand:

1. Open <https://pmc.ncbi.nlm.nih.gov/articles/PMC9922159/>
2. Download **`Supporting_Information.pdf`** (NIHMS1869546-supplement) and save
   it as `data/raw/SI.pdf`.
3. `make inspect` → dumps the extracted tables to
   `data/raw/si_tables_dump.csv` for you to confirm the layout.
4. Confirm the column mapping in `src/asyn/ingest.py::_map_rows`, then
   `make data` to build `data/aggregation_curated.csv`.

Until then, `make data-synth` writes a labelled synthetic stand-in so the whole
analysis runs. See `data/data_dictionary.md` for the schema and provenance.

## Run it
```bash
make setup        # pip install -e .[dev]
make all          # synthetic data + full analysis (or: make data && make analysis)
make check        # tests + analysis
```
Outputs: `reports/results.md`, and figures under `reports/figures/`.

## What the code does (maps to the phased plan)
| Module | Phase | Role |
|---|---|---|
| `asyn.ingest` | 1 | Build curated CSV from the real SI PDF (no fabrication). |
| `asyn.synth` | — | Synthetic stand-in in the same schema. |
| `asyn.features` | 3 | Three independent blocks: sequence, structure, conditions. |
| `asyn.splits` | guardrails | Leave-one-study-out + **fold-local** impute/scale/encode. |
| `asyn.decompose` | 4 | Mixed-effects study ICC + nested-CV R² per block set. **Main result.** |
| `asyn.predict` | 5 | Nested LOSO predictor, three baselines, bootstrap CIs. |
| `asyn.report` | 2/4/5/6 | The four figures. |
| `asyn.pipeline` | — | One-command rerun → `reports/results.md`. |

## Methodology guardrails (enforced in code)
- **Split by study, not by row** — `LeaveOneGroupOut` on `study`; a leakage
  assert is in the tests.
- **Scale inside the fold** — every transform is in the sklearn `Pipeline`;
  tested that the scaler uses train stats only.
- **Nested CV** — inner `GroupKFold` tunes, outer loop scores.
- **Three baselines, always** — mean / existing-tool score / one feature.
- **Bootstrap CIs over studies** — on every reported metric.
- **No neural network** — n ≈ 130; the default estimator is a small forest.
- **Pre-registered** — see `preregistration.md`, committed before Phase 4.

## Limitations
- Outcome is T½, which conflates nucleation and elongation (Ohgita et al. 2022
  separate them).
- The structure block carries only `in_fibril_core` until a DSSP-derived
  `data/structure_features.csv` (from PDB 6CU7/6A6B) is supplied.
- The existing-tool baseline needs a `tool_score` column (TANGO/WALTZ/AGGRESCAN
  output) added to the data.

## Layout
```
data/   raw/ (SI.pdf, gitignored) · aggregation_curated.csv · data_dictionary.md
src/asyn/   the modules above
tests/   feature/split/decompose unit tests (fast, run on synthetic data)
reports/   results.md + figures/
preregistration.md
```
