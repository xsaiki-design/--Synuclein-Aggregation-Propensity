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
