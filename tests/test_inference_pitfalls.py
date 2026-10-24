import numpy as np
import pytest
import pandas as pd
from scipy import stats
from inference_pitfalls import target_count, paired_summary, make_pairs, leakage_case


def test_counts_hand_calculated():
    assert target_count(1_000_000, 20_000, 100_000) == 200_000
    assert target_count(400_000, 35_000, 100_000) == 140_000


@pytest.mark.parametrize("args", [(1, 0, 1), (0, 1, 1), (1, 1, 1)])
def test_valid_boundaries(args):
    assert 0 <= target_count(*args) <= args[0]


@pytest.mark.parametrize(
    "args",
    [
        (1, 1, 0),
        (-1, 1, 2),
        (1, -1, 2),
        (1, 3, 2),
        (1, 0.5, 2),
        (1, 1, 2.5),
        (True, 1, 2),
        (1, None, 2),
        (np.inf, 1, 2),
        (1, np.nan, 2),
    ],
)
def test_bad_counts(args):
    with pytest.raises(ValueError):
        target_count(*args)


def test_pair_units_and_independent_arithmetic():
    result, means = paired_summary(make_pairs())
    d = means.treated - means.control
    assert result["n_donors"] == 8
    assert result["df"] == 7
    assert np.isclose(result["t"], stats.ttest_1samp(d, 0).statistic)


def test_row_order_invariant():
    a, _ = paired_summary(make_pairs())
    b, _ = paired_summary(make_pairs().sample(frac=1, random_state=2))
    assert a == b


@pytest.mark.parametrize(
    "problem", ["duplicate", "missing_pair", "nan", "condition", "empty"]
)
def test_bad_pairs(problem):
    d = make_pairs()
    if problem == "duplicate":
        d = pd.concat([d, d.iloc[:1]], ignore_index=True)
    if problem == "missing_pair":
        d = d[~((d.donor == "D01") & (d.condition == "treated"))]
    if problem == "nan":
        d.loc[0, "value"] = np.nan
    if problem == "condition":
        d.loc[0, "condition"] = "unknown"
    if problem == "empty":
        d = d.iloc[:0]
    with pytest.raises(ValueError):
        paired_summary(d)


def test_wrong_unit_changes_uncertainty():
    d = make_pairs()
    right, _ = paired_summary(d)
    w = d.pivot(index=["donor", "replicate"], columns="condition", values="value")
    wrong = stats.ttest_rel(w.treated, w.control)
    assert wrong.pvalue < right["p_t"]
    assert len(w) > right["n_donors"]


def test_donor_split_rejects_identity_shortcut():
    x = leakage_case()
    assert max(x["grouped_overlap_counts"]) == 0
    assert (
        x["sample_split_balanced_accuracy"] - x["donor_split_balanced_accuracy"] > 0.2
    )


def test_conversion_overflow_rejected():
    with pytest.raises(ValueError):
        target_count(10**400, 1, 1)
