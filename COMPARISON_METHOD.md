# Cross-study comparison: method, and what it can't yet do

How a second source was added, what the comparison showed, and the specific
weaknesses to fix before adding a third. Written after doing it once, so the
limitations are observed rather than anticipated.

## What was added

`asyn.ingest_flagmeier` → `data/cross_study_effects.csv` (15 rows).

[Flagmeier et al. 2016, PNAS 113(37):10328](https://pmc.ncbi.nlm.nih.gov/articles/PMC5027465)
covers the same familial mutants as Pancoe but reports **orders of magnitude
relative to WT**, separated by mechanism. Each row carries the verbatim
sentence it came from, because the paper states these in prose and does not
tabulate rate constants.

| | Pancoe 2022 | Flagmeier 2016 |
|---|---|---|
| value | categorical score in {−1,0,+1}, averaged | log10 fold vs WT |
| mechanism | lumped | nucleation / amplification / elongation, separate |
| conditions | mixed, mostly bulk ThT | lipid vesicles, or seeded at acidic pH |
| normalised to WT | yes, within study | yes |
| extraction | decoded from PDF text layer | transcribed from prose |

## What the comparison showed

Pancoe's per-mutant score against Flagmeier's lipid-induced nucleation:

| mutant | Pancoe | Flagmeier (log10) | same direction |
|---|---|---|---|
| A30P | −0.57 | 0.0 | — (no effect reported) |
| A53T | +0.91 | +2.0 | yes |
| E46K | +0.70 | **−2.0** | **no** |
| G51D | −0.50 | −2.0 | yes |
| H50Q | +0.85 | **−1.0** | **no** |

**Two of four decidable mutants flip sign, and the flips are large** — E46K
scores as an accelerator in the bulk literature but is reported two orders of
magnitude *slower* at lipid-induced nucleation.

**But 2/4 is not a result.** Binomial p = 1.0; Spearman ρ = 0.36 on five
points. The individual discrepancies are real and worth reporting as cases;
the overall pattern is underpowered, and the code now says so
(`distinguishable_from_chance`).

### Effect on the project's three claims

| claim | testable from this source | outcome |
|---|---|---|
| NAC is the most impactful region | **no** | every familial mutant is N-terminal (30, 46, 50, 51, 53) — no NAC or C-terminal variant is covered |
| X → Pro is the strongest substitution | **no** | only A30P is a proline variant. That one point disagrees: Pancoe scores N-terminal X → Pro at −0.81, Flagmeier reports no notable A30P effect (n=1) |
| sequence matters more than condition | **yes** | the sign flips occur *between assay families*, which is the condition effect the WT-normalised analysis could not see — though not yet above chance |

**The first two claims rest on Pancoe alone and remain uncorroborated.** That
is a scope limitation worth stating plainly: the entire familial-mutant
analysis is an N-terminus analysis.

## Decision: per-paper analysis is the primary unit

Cross-study agreement was tried first and is underpowered (p = 1.0 at n = 4),
with endpoint and condition mismatches that cannot be separated from real
disagreement. Analysing each paper **on its own terms** yields more, because
within one paper the lab, protein prep and protocol are held fixed.

The clearest example, from Flagmeier alone — the same five mutants, measured by
the same group:

| mechanism | spread across the 5 mutants | fold range |
|---|---|---|
| lipid-induced nucleation | 4.0 log10 | **10,000×** |
| secondary nucleation | 2.0 log10 (censored) | 100× |
| elongation | 1.0 log10 | 10× |

**How much a mutation matters depends on which mechanism is measured — by three
orders of magnitude.** No cross-study comparison is involved, so none of the
comparability caveats apply.

This reframes the project's third claim. "Does sequence matter more than
condition?" is not well posed against a single lumped outcome: the same
substitutions are decisive for membrane-induced nucleation and nearly
irrelevant for elongation. A lumped score such as Pancoe's averages over that
distinction, which is one concrete mechanism by which a real condition effect
becomes invisible.

**Operating rule going forward.** Each paper gets its own ingester and its own
mechanism-stratified profile. `cross_study_effects.csv` remains the provenance
registry and the basis for `claim_checks` (which claims a source can even
address). Cross-paper discordance is reported as **named case studies** — E46K
and H50Q below — never as an agreement rate, until the pool is large enough for
a rate to mean something.

## Weaknesses to fix before adding a third paper

**1. Independence is unverified — check it first.** If a paper is already
inside Pancoe's 375-example compilation, it is not independent evidence, and
agreement is partly circular. Flagmeier's inclusion status is **unknown**: the
PMC page for Pancoe returned a 167-character stub, so the reference list could
not be read. Resolve this per paper before counting it as corroboration.

**2. Power. Five mutants cannot settle anything.** Sign tests need roughly
15–20 paired points before they discriminate. Two consequences: prefer sources
covering many variants over sources covering the same six, and report small-n
comparisons as case studies, never as rates.

**3. Endpoint equivalence is currently a judgement call.** Pancoe's lumped
score was compared against lipid-induced nucleation because both answer "does
aggregation happen sooner". That choice was not pre-specified, and with three
candidate mechanisms there were three chances to find agreement. Declare an
explicit compatibility mapping before looking at results, and record it.

**4. Condition matching is absent.** Pancoe's score averages over buffer, salt,
shaking and concentration; Flagmeier's is lipid-specific. Like is not being
compared with like, which is itself a plausible explanation for the sign flips
and cannot currently be separated from a real disagreement.

**5. Censoring is faked.** "Undetectably slow" became −2.0. That is a floor
standing in for a bound, flagged by `censored` but not yet honoured in the
statistics. Either handle bounds properly or exclude those rows from numeric
comparisons.

**6. Uncertainty is asymmetric.** Pancoe rows carry support counts and SDs;
Flagmeier rows carry nothing, because prose has no error bars. Weighted
comparison is impossible until a source provides them, so weighting is
currently uniform, which overstates the prose-derived values.

**7. Extraction reliability differs and is not used.** `extraction` records
`text_statement` vs decoded table, but nothing down-weights the less reliable
kind. A figure-digitised value, when one arrives, will be a third tier again.

## Schema for adding a paper

Every row in `data/cross_study_effects.csv`:

| column | why it exists |
|---|---|
| `study`, `year`, `source` | provenance; `study` is the grouping key — never pool across it |
| `mutant` | join key |
| `endpoint` | **the thing actually measured** — never compare across these without a declared mapping |
| `condition_class` | bulk / lipid_vesicles / seeded_acidic / cellular |
| `normalized_to_wt` | the whole project turned on this; mixing raw and normalised silently reproduces the original error |
| `value`, `value_type` | the number and its units |
| `censored` | true when the source gives a bound, not a value |
| `extraction` | `text_statement` / `table_decoded` / `figure_digitised` |
| `quote` | verbatim sentence supporting the value, so every number is auditable |

**Add one module per paper** (`ingest_<study>.py`), each exposing `build()`,
each carrying its quotes. That pattern scales; a single hand-maintained CSV
would not.

### Papers queued, with the endpoint caveat that governs each

| paper | endpoint | note |
|---|---|---|
| Stephens 2020 (Nat Commun) | ThT kinetics + N-terminus exposure | raw data deposited; exposure is a *structural feature* the predictor lacks |
| Lázaro 2014 (PLoS Genet) | cellular inclusions, 19 variants | **cellular, not in vitro** — separate axis, never pooled |
| Oh 2020 (Front Cell Neurosci) | transmission, aggregation, 8 mutants | independently reports H50Q promotes / G51D delays, matching our extremes |
| Buell 2014 (PNAS) | condition dependence of mechanism | the control on claim 3 — conditions set absolute rate |

The recurring hazard is endpoint drift: toxicity, inclusion counts and
transmission are *not* aggregation kinetics. Pooling them would grow n while
quietly averaging incomparable quantities — a worse error than small n.
