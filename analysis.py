"""Reference problems in cell-therapy assay mathematics; synthetic values only."""

from pathlib import Path
import json, math
import numpy as np
import pandas as pd
from scipy import stats, optimize
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).parent


def nonnegative(x):
    """Accept finite real numbers; reject booleans, strings and conversion overflow."""
    if isinstance(x, (bool, np.bool_)) or not isinstance(
        x, (int, float, np.integer, np.floating)
    ):
        raise ValueError("Expected finite nonnegative numeric value")
    try:
        value = float(x)
    except (ValueError, TypeError, OverflowError) as exc:
        raise ValueError("Numeric conversion failed") from exc
    if not math.isfinite(value) or value < 0:
        raise ValueError("Expected finite nonnegative numeric value")
    return value


def probability(x, closed=False):
    x = nonnegative(x)
    if not (0 <= x <= 1 if closed else 0 < x < 1):
        raise ValueError("Invalid probability")
    return x


def count(x):
    x = nonnegative(x)
    if not x.is_integer():
        raise ValueError("Counts must be integers")
    return int(x)


def complex_concentration(receptor, ligand, kd):
    r, l, k = map(nonnegative, [receptor, ligand, kd])
    if min(r, l) == 0:
        return 0.0
    if k == 0:
        return min(r, l)
    scale = max(r, l, k)
    rs, ls, ks = r / scale, l / scale, k / scale
    s = rs + ls + ks
    discriminant = (rs - ls) ** 2 + ks * ks + 2 * ks * (rs + ls)
    return scale * (2 * rs * ls) / (s + math.sqrt(discriminant))


def occupancy(receptor, ligand, kd):
    r = nonnegative(receptor)
    if r == 0:
        raise ValueError("Occupancy undefined without receptor")
    return complex_concentration(r, ligand, kd) / r


def ratio_interval(target, beads, confidence=0.95):
    t, b = count(target), count(beads)
    confidence = probability(confidence)
    if t + b == 0:
        raise ValueError("No observed events")
    alpha = 1 - confidence
    lo = 0.0 if t == 0 else stats.beta.ppf(alpha / 2, t, b + 1)
    hi = 1.0 if b == 0 else stats.beta.ppf(1 - alpha / 2, t + 1, b)
    return lo / (1 - lo), np.inf if hi == 1 else hi / (1 - hi)


def bead_count(target, bead_events, beads_added, volume_ul):
    t, b = count(target), count(bead_events)
    added, vol = nonnegative(beads_added), nonnegative(volume_ul)
    if b == 0 or added == 0 or vol == 0:
        raise ValueError("Positive bead count, beads added and volume required")
    factor = added / vol
    lo, hi = ratio_interval(t, b)
    return {"cells_per_ul": t / b * factor, "ci95": [lo * factor, hi * factor]}


def minimum_cells_no_background(prevalence, sensitivity=0.95, confidence=0.95):
    p = probability(prevalence, True)
    s = probability(sensitivity, True)
    c = probability(confidence)
    effective = p * s
    if effective == 0:
        raise ValueError("No finite sample size can detect a zero-probability event")
    if effective == 1:
        return 1
    n = math.ceil(math.log1p(-c) / math.log1p(-effective))
    # Correct possible final-bit rounding without changing the statistical model.
    while -math.expm1(n * math.log1p(-effective)) < c:
        n += 1
    while n > 1 and -math.expm1((n - 1) * math.log1p(-effective)) >= c:
        n -= 1
    return n


def rare_event_design(n, prevalence, sensitivity, false_positive, alpha=0.05):
    n = count(n)
    p = probability(prevalence, True)
    s = probability(sensitivity, True)
    f = probability(false_positive, True)
    alpha = probability(alpha)
    if n == 0:
        raise ValueError("Need observed cells")
    threshold = int(stats.binom.ppf(1 - alpha, n, f)) + 1
    observed_positive = p * s + (1 - p) * f
    return {
        "n_cells": n,
        "positive_threshold": threshold,
        "blank_false_positive_probability": float(stats.binom.sf(threshold - 1, n, f)),
        "power": float(stats.binom.sf(threshold - 1, n, observed_positive)),
        "naive_one_positive_blank_rate": float(stats.binom.sf(0, n, f)),
    }


def run(output=ROOT / "artifacts/run"):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    grid = np.logspace(-2, 3, 150)
    r, k = 100.0, 1.0
    bound = np.array([occupancy(r, l, k) for l in grid])
    approx = grid / (k + grid)
    # Independent numerical root checks use the original mass-action equation.
    for ligand in [0.01, 1, 10, 100, 1000]:
        exact = complex_concentration(r, ligand, k)
        root = optimize.brentq(
            lambda c: (r - c) * (ligand - c) - k * c, 0, min(r, ligand)
        )
        assert np.isclose(exact, root, rtol=1e-9, atol=1e-10)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.semilogx(grid, bound, label="Finite total ligand: mass balance")
    ax.semilogx(
        grid, approx, "--", label="Incorrectly substitute total for free ligand"
    )
    ax.set(
        xlabel="Total ligand (nM)",
        ylabel="Receptor occupancy",
        ylim=(0, 1.05),
        title="Ligand depletion can invalidate the usual occupancy shortcut",
    )
    ax.legend()
    fig.tight_layout()
    fig.savefig(output / "binding.png", dpi=160)
    plt.close(fig)
    point = bead_count(120, 800, 10_000, 100)
    rare = rare_event_design(30000, 0.0001, 0.95, 0.0001)
    no_bg = minimum_cells_no_background(0.0001, 0.95, 0.95)
    rows = []
    for n in [5000, 10000, 20000, 30000, 50000, 100000]:
        rows.append(rare_event_design(n, 0.0001, 0.95, 0.0001))
    pd.DataFrame(rows).to_csv(output / "rare_event_power.csv", index=False)
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(
        [r["n_cells"] for r in rows],
        [r["power"] for r in rows],
        marker="o",
        label="Blank-calibrated detection power",
    )
    ax.plot(
        [r["n_cells"] for r in rows],
        [r["naive_one_positive_blank_rate"] for r in rows],
        marker="x",
        label="False-positive rate: any positive",
    )
    ax.set(
        xlabel="Cells analyzed",
        ylabel="Probability",
        ylim=(0, 1),
        title="Background changes the rare-cell detection problem",
    )
    ax.legend()
    fig.tight_layout()
    fig.savefig(output / "rare_event_power.png", dpi=160)
    plt.close(fig)
    result = {
        "source": "Synthetic reference problems, not product-specific assay specifications",
        "binding_example": {
            "receptor_nM": 100,
            "total_ligand_nM": 10,
            "kd_nM": 1,
            "complex_nM": complex_concentration(100, 10, 1),
            "exact_occupancy": occupancy(100, 10, 1),
            "incorrect_free_ligand_shortcut": 10 / 11,
        },
        "bead_count_example": point,
        "no_background_cells_for_95_percent_detection": no_bg,
        "background_aware_design": rare,
        "assumptions": [
            "Binding assumes one-to-one equilibrium in a well-mixed volume and consistent concentration units; no avidity or kinetics.",
            "Bead interval conditions on total observed target+bead events; known bead addition/volume, independent Poisson sampling and representative mixing assumed. It excludes pipetting and viability uncertainty.",
            "Rare-event calculations assume independent Bernoulli cell classifications, fixed sensitivity and false-positive rates. They do not establish an experimental LOD, product release threshold, or prevalence estimate.",
        ],
    }
    (output / "summary.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n"
    )
    return result


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
