# Pre-registration

Written **before** running Phase 4. If any point below is changed after seeing
results, the change and its reason are recorded in the commit history; the
original stays here.

## Question
How much of the variation in reported α-synuclein aggregation rate (T½) is
explained by the protein's sequence versus the experimental conditions each
measurement was made under?

## Primary analysis (Phase 4)
Mixed-effects null model `log10(t_half) ~ 1 + (1 | study)`. The pre-registered
primary quantity is the **study-level ICC** — the fraction of total variance
absorbed by the study random intercept. Reported alongside: leave-one-study-out
CV R² for four nested feature sets (conditions-only, sequence-only,
sequence+structure, everything).

## Prediction (Phase 5)
- Outer validation: **leave-one-study-out**. No random row splits.
- Hyperparameters: tuned by **inner GroupKFold only** (nested CV). Outer folds
  only score.
- All imputation/scaling/encoding fit **inside the training fold**.
- Estimator: a small random forest. **No neural network** (n ≈ 130).
- Baselines, always reported: (1) predict the training-fold mean; (2) an
  existing tool's score (TANGO/WALTZ) when available; (3) a single feature,
  hydrophobicity change.
- Uncertainty: **bootstrap 95% CI resampled over whole studies** on every metric.

## External test (Phase 6)
Ulamec et al. 2022 P1-region scan, held out from the start, never touched
during development.

## Decision rule
If the full model does not clear all available baselines under
leave-one-study-out by a margin outside the bootstrap CI, that is reported as
the result — the complexity is not justified.


---

## Deviations from this pre-registration (recorded after seeing the data)

**1. The planned outcome variable does not exist in the source.** This plan
plus the primary mixed-effects model assumed a continuous T½ per experiment.
The Pancoe et al. 2022 SI tabulates a categorical aggregation score in
{−1, 0, +1}, already normalised to a matched WT control, and no per-experiment
T½. The planned `log10(t_half) ~ 1 + (1|study)` model and the leave-one-study-out
predictor therefore could not be run on this source. They remain implemented
(`asyn.decompose`, `asyn.predict`) for raw kinetics if that data is obtained.

**2. Substituted analysis.** The sequence-vs-conditions question was instead
answered at the level the SI supports: support-weighted variance in the
aggregation score attributed to mutant identity versus condition category,
over Tables S7–S10 (`asyn.si_analysis`).

**3. Added a permutation null (not in the original plan).** Raw R² is not
comparable across those two factors — mutant carries 5 degrees of freedom
against a category's 1–2, on 14–19 cells, so the larger factor explains more
variance mechanically. Each factor is now compared with its own shuffled-label
baseline. This was added because the first pass overstated the contrast, and it
changed the claim: condition category does not merely explain less than
sequence, it falls *below* its own chance level.

**4. Bootstrap CIs not reported for this decomposition.** The plan's bootstrap
requirement was written for the Phase 5 predictor. With 6 mutants, resampling
clusters frequently yields rank-deficient designs, so the permutation test is
reported in its place rather than a fragile interval.
