"""alpha-synuclein aggregation propensity: sequence vs. conditions.

See README.md for the question and the phased plan. Public entry points:
  asyn.data       schema, canonical sequence, domains, load_curated()
  asyn.ingest     build the curated CSV from the real supplementary PDF
  asyn.synth      synthetic stand-in dataset (clearly labelled, not real)
  asyn.features   the three feature blocks
  asyn.splits     study-grouped CV + fold-local preprocessing
  asyn.decompose  Phase 4 variance decomposition
  asyn.predict    Phase 5 validated predictor + baselines
  asyn.report     figures
  asyn.pipeline   one-command rerun
"""
__all__ = [
    "data", "ingest", "synth", "features", "splits",
    "decompose", "predict", "report", "pipeline",
]
