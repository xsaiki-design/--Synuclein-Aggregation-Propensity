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
