from asyn import decompose, synth


def test_study_icc_is_a_fraction():
    df = synth.generate(n_studies=8, seed=5)
    icc = decompose.study_icc(df)
    assert 0.0 <= icc["icc_study"] <= 1.0
    assert icc["n_studies"] == 8
    assert icc["n_obs"] == len(df)


def test_nested_cv_table_has_all_block_sets():
    df = synth.generate(n_studies=8, seed=6)
    table = decompose.nested_cv_r2(df)
    assert list(table["model"]) == list(decompose.BLOCK_SETS)
    assert (table["n_features"] > 0).all()
