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
    ve = si_analysis.run()["variance_explained"].set_index("dimension")
    row = ve.loc["ALL (pooled)"]
    assert row["r2_sequence_mutant"] > row["r2_conditions_category"]


def test_si_report_writes_figures():
    from asyn import si_report
    paths = si_report.run()
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
