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


def test_si_report_writes_figures(tmp_path, monkeypatch):
    from asyn import si_report
    # Write into a temp dir: tests must not clobber the repo's real figures,
    # which are generated at full permutation count by the CLI.
    monkeypatch.setattr(si_report, "FIG_DIR", tmp_path)
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


def test_whole_protein_rows_partition_across_regions():
    """X -> Pro's n=33 is an aggregate: 16 N-term + 17 NAC + 0 C-term."""
    t = si_analysis.load_types()
    x = t[t["type"] == "X -> Pro"].set_index("scope")["support"]
    assert x["whole_protein"] == 33
    assert x["N-terminus"] + x["NAC"] + x["C-terminus"] == x["whole_protein"]
    assert x["C-terminus"] == 0  # proline substitution untested in the C-terminus


def test_strongest_by_region_differs_in_sign():
    top = si_analysis.strongest_by_region(si_analysis.load_types()).set_index("scope")
    assert top.loc["NAC", "type"] == "X -> Pro"
    assert top.loc["N-terminus", "average"] > 0   # N-term strongest accelerates
    assert top.loc["C-terminus", "average"] < 0   # C-term strongest blocks
    assert (top["support"] >= 5).all()


def test_flagmeier_rows_carry_auditable_quotes():
    from asyn import ingest_flagmeier
    df = ingest_flagmeier.build()
    assert len(df) == 15
    assert (df["quote"].str.len() > 40).all()        # every value traceable
    assert (df["extraction"] == "text_statement").all()
    assert df["censored"].sum() == 2                 # the two "undetectably slow"


def test_cross_study_small_n_is_not_called_significant():
    """2/4 sign agreement must not be reported as a finding."""
    from asyn import cross_study
    a = cross_study.agreement_summary()
    assert not a["distinguishable_from_chance"].any()
    assert (a["binom_p"] >= 0.05).all()


def test_claims_untestable_from_flagmeier_are_marked_so():
    from asyn import cross_study
    c = cross_study.claim_checks().set_index("claim")
    assert not c.loc["NAC is the most impactful region", "testable_here"]
    assert not c.loc["X -> Pro is the strongest single substitution", "testable_here"]


def test_mechanism_spread_is_largest_for_nucleation():
    """The within-study finding: sequence sensitivity is mechanism-dependent."""
    from asyn import cross_study
    m = cross_study.mechanism_spread().set_index("endpoint")["spread_log10"]
    assert m["lipid_nucleation"] > m["secondary_nucleation"] > m["elongation"]
    assert m["lipid_nucleation"] >= 4.0   # 10,000-fold range


def test_vote_reconstruction_yields_whole_experiments():
    """mean+sd+n must decompose into integer counts — a third extraction check."""
    from asyn import discordance
    v = discordance.vote_split()
    assert v["integrality_error"].max() < 0.05
    assert (v[["n_accel", "n_noeffect", "n_reduce"]] >= 0).all().all()
    assert (v["n_accel"] + v["n_noeffect"] + v["n_reduce"]
            - v["support"]).abs().max() <= 1


def test_sd_never_breaches_the_categorical_bound():
    from asyn import discordance
    assert len(discordance.validate(discordance.load_types())) == 0


def test_proline_substitution_is_unanimous_in_direction():
    """No experiment anywhere reports X -> Pro accelerating aggregation."""
    from asyn import discordance
    v = discordance.vote_split()
    pro = v[v["type"] == "X -> Pro"]
    assert len(pro) >= 3
    assert (pro["n_accel"] == 0).all()
    assert (pro["split"] == 0).all()


@pytest.mark.parametrize("module", [
    "ingest_si", "ingest_types", "ingest_flagmeier",
    "si_analysis", "cross_study", "discordance",
])
def test_module_runs_as_a_script(module):
    """Each module documents `python -m asyn.<name>`; check it actually runs.

    si_report is excluded: its CLI writes the repo's real artifacts.
    """
    import subprocess, sys, pathlib
    r = subprocess.run([sys.executable, "-m", f"asyn.{module}"],
                       cwd=pathlib.Path(__file__).resolve().parent.parent,
                       env={"PYTHONPATH": "src", "PATH": "/usr/bin:/bin:/usr/local/bin"},
                       capture_output=True, text=True, timeout=300)
    assert r.returncode == 0, f"{module} CLI failed:\n{r.stderr[-800:]}"


def test_reported_r2_is_an_upper_bound():
    """Cell-mean R² must exceed the experiment-level lower bound, and the gap
    must be large enough that the headline cannot be quoted without it."""
    b = si_analysis.variance_explained_bounds()
    assert (b["r2_cell_means"] > b["r2_experiment_lower"]).all()
    assert (b["r2_experiment_lower"] < 0.6).all()
    assert (b["r2_cell_means"] >= 0.82).all()


def test_by_region_excludes_the_nested_whole_protein_row():
    """whole_protein is the sum of the three regions; listing it invites
    double counting (verified: exact for all 26 types)."""
    from asyn import discordance
    scopes = set(discordance.by_region()["scope"])
    assert scopes == {"N-terminus", "NAC", "C-terminus"}


def test_whole_protein_is_the_sum_of_the_regions():
    t = si_analysis.load_types()
    w = t[t.scope == "whole_protein"].set_index("type")["support"]
    r = t[t.scope != "whole_protein"].groupby("type")["support"].sum()
    assert (w.sort_index() == r.sort_index()).all()


def test_dispersion_ratio_is_gone():
    """score() misread extreme means as contested; vote_split replaces it."""
    from asyn import discordance
    assert not hasattr(discordance, "score")


def test_nac_is_robust_but_nterm_cterm_ordering_is_not():
    """NAC lowest holds with or without control rows; the other two do not."""
    t = si_analysis.load_types()
    withc = si_analysis.domain_means(t).set_index("scope")["weighted_mean_score"]
    without = si_analysis.domain_means(t, transitions_only=True).set_index(
        "scope")["weighted_mean_score"]
    assert withc["NAC"] < withc["C-terminus"] and withc["NAC"] < withc["N-terminus"]
    assert without["NAC"] < without["C-terminus"] and without["NAC"] < without["N-terminus"]
    # the N-term/C-term gap collapses once controls are excluded
    assert abs(without["N-terminus"] - without["C-terminus"]) < 0.05
    assert abs(withc["N-terminus"] - withc["C-terminus"]) > 0.2
