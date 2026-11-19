"""Guards on the real SI data path. Skips cleanly if the PDF-derived CSV
hasn't been built (data/conditions_familial.csv is committed, so it normally
runs)."""
from pathlib import Path

import pytest

from asyn import si_analysis

pytestmark = pytest.mark.skipif(
    not Path("data/conditions_familial.csv").exists(),
    reason="run `python -m asyn.ingest_si` first",
)


def test_a30p_buffer_matches_paper():
    out = si_analysis.run()
    # The paper's own headline numbers (SI Table S7).
    assert out["a30p_buffer"]["PBS"] == -0.75
    assert out["a30p_buffer"]["Tris"] == -0.31


def test_sequence_explains_more_than_conditions():
    ve = si_analysis.run()["variance_explained"]
    assert set(ve["dimension"]) == {"buffer", "salt", "shaking", "concentration"}
    assert (ve["r2_sequence_mutant"] > ve["r2_conditions_category"]).all()


def test_si_report_writes_figures():
    from asyn import si_report
    paths = si_report.run(n_perm=50)  # figures only; rigour lives in the CLI run
    assert len(paths) >= 2 and all(p.exists() and p.stat().st_size > 1000 for p in paths)


def test_type_scores_match_rendered_pages():
    """The ingester's own spot checks must hold (9 values read off the SI)."""
    from asyn import ingest_types
    import pandas as pd
    df = pd.read_csv(ingest_types.OUT_CSV)
    assert len(df) == 104 and df["scope"].nunique() == 4
    ingest_types._verify(df)  # raises on CID drift


def test_nac_is_the_most_aggregation_critical_region():
    d = si_analysis.load_types()
    m = si_analysis.domain_means(d).set_index("scope")["weighted_mean_score"]
    assert m["NAC"] < m["C-terminus"] < m["N-terminus"]


def test_condition_category_is_below_its_chance_level():
    """The corrected claim: conditions carry no signal, not merely less signal."""
    t = si_analysis.permutation_test(si_analysis.load(), n_perm=200, seed=1)
    # checked per table: the tables share experiments, so there is no valid pool
    assert (t["r2_mutant"] > t["chance_mutant"]).all()       # sequence beats chance
    assert (t["r2_category"] < t["chance_category"]).all()   # conditions do not
    assert (t["p_mutant"] < 0.05).all()


def test_none_category_is_not_parsed_as_missing():
    """'None' (no added salt / no shaking) is a real level and is in pandas'
    default NA list; losing it silently shrank the category model."""
    d = si_analysis.load()
    assert "None" in set(d["category"])
    assert d["category"].isna().sum() == 0


def test_conditions_data_matches_rendered_pages():
    from asyn import ingest_si
    ingest_si._verify(si_analysis.load())  # raises on CID drift
