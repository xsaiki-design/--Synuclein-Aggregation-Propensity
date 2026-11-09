# Data dictionary — `aggregation_curated.csv`

One row per *reported aggregation experiment*. Built from the open-access
supplementary tables of **Pancoe et al. 2022** (J Mol Biol 434(23):167859,
[PMC9922159](https://pmc.ncbi.nlm.nih.gov/articles/PMC9922159)), whose authors
explicitly intended the tabulated data for computational reuse.

The tabulated outcome is **T½** — the time to half-maximal ThT fluorescence.
Note for the write-up: T½ conflates nucleation and elongation into one number;
papers that separate those phases (e.g. Ohgita et al. 2022) are a known
limitation of using a single half-time as the outcome.

| Column | Type | Meaning |
|---|---|---|
| `variant` | str | Mutation in WT-position-MUT form (`A30P`), or `WT`. |
| `study` | str | The paper/lab the row was measured in. **The grouping key** for leave-one-study-out validation. |
| `position` | int | 1-based residue number; empty for WT. |
| `wt_residue` | str | Wild-type residue (1-letter); empty for WT. |
| `mut_residue` | str | Mutant residue (1-letter); empty for WT. |
| `t_half_h` | float | **Outcome.** Time to half-maximal ThT signal, hours. |
| `t_half_rel_wt` | float | `t_half_h` normalised to the WT control measured in the same study. |
| `pH` | float | Assay pH. |
| `temperature_C` | float | Assay temperature (°C). |
| `protein_conc_uM` | float | Monomer concentration (µM). |
| `agitation` | str | `quiescent` / `shaking` / `beads`. |
| `buffer` | str | Buffer system (e.g. PBS, Tris, phosphate). |
| `seeds` | 0/1 | Pre-formed fibril seeds added. |
| `lipids` | 0/1 | Lipid/membrane present. |
| `ThT_conc_uM` | float | Thioflavin T concentration (µM). |
| `synthetic` | 0/1 | **Present only in the stand-in.** `1` ⇒ generated, not real. |

**Missing conditions are data.** If a paper does not report a condition for a
row, that cell stays empty (NaN) — never imputed at curation time, and the
frequency of missingness is itself reported.

## Provenance / target variable
- Outcome metric and its caveat: T½ (half-time), conflates nucleation +
  elongation. Separate-phase kinetics: Ohgita et al. 2022 (PMC9043213).
- Canonical sequence for feature construction: UniProt **P37840**, 140 aa.
- Fibril structure features (optional block): PDB **6CU7 / 6A6B** via DSSP →
  `data/structure_features.csv` (`position, sasa, ss_sheet`). Without it, the
  structure block carries only `in_fibril_core` (NAC = residues 61–95).
