import numpy as np
import pandas as pd

from asyn import data, features, synth


def test_canonical_sequence_length():
    assert len(data.ASYN_SEQUENCE) == 140


def test_domain_and_distance():
    assert data.domain_of(30) == "amphipathic"
    assert data.domain_of(70) == "nac"
    assert data.domain_of(130) == "acidic"
    assert data.domain_of(float("nan")) == "wt"
    assert data.distance_to_nac(70) == 0.0          # inside NAC
    assert data.distance_to_nac(60) == 1.0          # one before NAC start (61)
    assert data.distance_to_nac(96) == 1.0          # one after NAC end (95)


def test_sequence_block_wt_has_zero_deltas():
    df = pd.DataFrame([{"wt_residue": np.nan, "mut_residue": np.nan, "position": np.nan}])
    out = features.sequence_block(df)
    for c in ["d_hydro", "d_charge", "d_volume", "d_helix", "d_sheet"]:
        assert out[c].iloc[0] == 0.0


def test_sequence_block_known_mutation():
    # A30P: Ala -> Pro. Hydropathy drops (1.8 -> -1.6), charge unchanged.
    df = pd.DataFrame([{"wt_residue": "A", "mut_residue": "P", "position": 30}])
    out = features.sequence_block(df)
    assert out["d_hydro"].iloc[0] == features.HYDRO["P"] - features.HYDRO["A"]
    assert out["d_charge"].iloc[0] == 0.0
    assert out["dom_amphipathic"].iloc[0] == 1


def test_build_features_shapes_and_blocks():
    df = synth.generate(n_studies=4, seed=1)
    X, groups, y, blocks = features.build_features(df)
    assert len(X) == len(df) == len(y) == len(groups)
    assert set(blocks) == {"sequence", "structure", "conditions"}
    assert list(X.columns) == blocks["sequence"] + blocks["structure"] + blocks["conditions"]
    assert np.isfinite(y).all()
