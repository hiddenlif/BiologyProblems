import numpy as np
import pytest
from scipy import optimize, stats
from analysis import (
    complex_concentration,
    occupancy,
    ratio_interval,
    bead_count,
    minimum_cells_no_background,
    rare_event_design,
)


@pytest.mark.parametrize(
    "r,l,k",
    [
        (100, 10, 1),
        (1, 1, 1),
        (0, 2, 3),
        (2, 0, 3),
        (2, 3, 0),
        (1e-12, 1e-9, 1e-10),
        (1e150, 1e149, 1e148),
    ],
)
def test_binding_constraints_and_mass_action(r, l, k):
    c = complex_concentration(r, l, k)
    assert 0 <= c <= min(r, l)
    if max(r, l, k) > 0:
        assert np.isclose(
            ((r - c) / max(r, l, k)) * ((l - c) / max(r, l, k)),
            (k / max(r, l, k)) * (c / max(r, l, k)),
            atol=1e-12,
        )


def test_binding_independent_root():
    r, l, k = 100, 10, 1
    root = optimize.brentq(lambda c: (r - c) * (l - c) - k * c, 0, 10)
    assert np.isclose(complex_concentration(r, l, k), root)


def test_ligand_depletion_exposes_wrong_shortcut():
    assert occupancy(100, 10, 1) < 0.1 and 10 / 11 > 0.9


@pytest.mark.parametrize("bad", [True, None, -1, np.nan, np.inf, "1"])
def test_bad_binding_input(bad):
    with pytest.raises(ValueError):
        complex_concentration(bad, 1, 1)


def test_zero_receptor_occupancy_undefined():
    with pytest.raises(ValueError):
        occupancy(0, 1, 1)


def test_count_hand_solution():
    assert bead_count(120, 800, 10000, 100)["cells_per_ul"] == 15


@pytest.mark.parametrize("args", [(0, 800, 10000, 100), (120, 800, 10000, 100)])
def test_count_interval_contains_estimate(args):
    x = bead_count(*args)
    assert x["ci95"][0] <= x["cells_per_ul"] <= x["ci95"][1]


@pytest.mark.parametrize(
    "args", [(1, 0, 1, 1), (1, 1, 0, 1), (1, 1, 1, 0), (0.5, 1, 1, 1), (True, 1, 1, 1)]
)
def test_bad_count(args):
    with pytest.raises(ValueError):
        bead_count(*args)


@pytest.mark.parametrize("n,p", [(5, 0.1), (10, 0.5), (30, 0.9)])
def test_exact_conditional_coverage(n, p):
    covered = 0
    for t in range(n + 1):
        lo, hi = ratio_interval(t, n - t)
        truth = p / (1 - p)
        if lo <= truth <= hi:
            covered += stats.binom.pmf(t, n, p)
    assert covered >= 0.95 - 1e-12


def test_sample_size_minimal():
    n = minimum_cells_no_background(0.0001, 0.95, 0.95)
    q = 0.0001 * 0.95
    assert stats.binom.sf(0, n, q) >= 0.95
    assert stats.binom.sf(0, n - 1, q) < 0.95


@pytest.mark.parametrize("p,s", [(0, 1), (1, 0)])
def test_impossible_detection(p, s):
    with pytest.raises(ValueError):
        minimum_cells_no_background(p, s)


def test_certain_detection():
    assert minimum_cells_no_background(1, 1) == 1


def test_background_control_and_power():
    d = rare_event_design(30000, 0.0001, 0.95, 0.0001)
    assert d["blank_false_positive_probability"] <= 0.05
    assert 0 <= d["power"] <= 1
    assert d["naive_one_positive_blank_rate"] > 0.9


@pytest.mark.parametrize("p", [np.nan, -0.1, 1.1, True])
def test_invalid_probability(p):
    with pytest.raises(ValueError):
        rare_event_design(100, p, 0.9, 0.001)


def test_conversion_overflow_rejected():
    with pytest.raises(ValueError):
        complex_concentration(10**400, 1, 1)
