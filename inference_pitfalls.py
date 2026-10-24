"""Inference pitfalls in cell-therapy analysis: absolute counts, technical replication and donor leakage. Controlled synthetic cases; not patient data or evidence of treatment efficacy."""

from pathlib import Path
import itertools, json
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import GroupKFold, StratifiedKFold, cross_val_predict
from sklearn.metrics import balanced_accuracy_score

ROOT = Path(__file__).parent


def number(x, positive=False):
    """Validate physical nonnegative values and strictly positive denominators."""
    if isinstance(x, (bool, np.bool_)) or not isinstance(
        x, (int, float, np.integer, np.floating)
    ):
        raise ValueError("Expected a real numeric value")
    try:
        value = float(x)
    except (ValueError, TypeError, OverflowError) as exc:
        raise ValueError("Numeric conversion failed") from exc
    if not np.isfinite(value) or (value <= 0 if positive else value < 0):
        raise ValueError(
            "Expected finite nonnegative value; denominator must be positive"
        )
    return value


def target_count(total_viable_cells, target_events, parent_events):
    total = number(total_viable_cells)
    target = number(target_events)
    parent = number(parent_events, True)
    if not target.is_integer() or not parent.is_integer() or target > parent:
        raise ValueError("Invalid event counts")
    return total * target / parent


def paired_summary(df):
    required = {"donor", "replicate", "condition", "value"}
    if not required <= set(df) or df.empty or df[list(required)].isna().any().any():
        raise ValueError("Missing paired data")
    if df.duplicated(["donor", "replicate", "condition"]).any() or set(
        df.condition
    ) != {"control", "treated"}:
        raise ValueError("Bad pairing keys")
    if not pd.api.types.is_numeric_dtype(df.value) or not np.isfinite(df.value).all():
        raise ValueError("Invalid measurements")
    # Biological donor, not individual technical well, is the independent unit.
    donor = df.groupby(["donor", "condition"]).value.mean().unstack()
    if donor.isna().any().any() or len(donor) < 2:
        raise ValueError("Need complete biological donor pairs")
    d = donor.treated - donor.control
    n = len(d)
    se = d.std(ddof=1) / np.sqrt(n)
    if se == 0:
        raise ValueError("Zero between-donor variation: ordinary t inference undefined")
    t = d.mean() / se
    ci = stats.t.interval(0.95, n - 1, loc=d.mean(), scale=se)
    scipy = stats.ttest_rel(donor.treated, donor.control)
    if not np.isclose(t, scipy.statistic):
        raise AssertionError("Independent t arithmetic disagrees")
    signs = np.array(list(itertools.product([-1, 1], repeat=n)))
    p = (np.abs(signs @ d.to_numpy() / n) >= abs(d.mean()) - 1e-12).mean()
    return {
        "n_donors": n,
        "mean_difference": float(d.mean()),
        "t": float(t),
        "df": n - 1,
        "p_t": float(scipy.pvalue),
        "ci95": list(ci),
        "exact_sign_flip_p": float(p),
    }, donor


def make_pairs(seed=20261022):
    rng = np.random.default_rng(seed)
    rows = []
    effects = [-2, -1, 0, 0, 0, 1, 2, 8]
    for i, effect in enumerate(effects):
        baseline = 50 + rng.normal(0, 4)
        for rep in range(12):
            for condition, offset in [("control", 0), ("treated", effect)]:
                rows.append(
                    {
                        "donor": f"D{i+1:02}",
                        "replicate": rep,
                        "condition": condition,
                        "value": baseline + offset + rng.normal(0, 0.08),
                    }
                )
    return pd.DataFrame(rows)


def leakage_case(seed=20261022):
    rng = np.random.default_rng(seed)
    n = 80
    reps = 8
    # Donor labels are deliberately unrelated to any across-donor biological signal.
    labels = np.tile([0, 1], n // 2)
    rng.shuffle(labels)
    centers = rng.normal(size=(n, 20))
    groups = np.repeat(np.arange(n), reps)
    y = np.repeat(labels, reps)
    x = np.repeat(centers, reps, axis=0) + rng.normal(0, 0.02, (n * reps, 20))
    model = make_pipeline(StandardScaler(), KNeighborsClassifier(3, n_jobs=1))
    random = StratifiedKFold(4, shuffle=True, random_state=seed)
    grouped = GroupKFold(4)
    rp = cross_val_predict(model, x, y, cv=random)
    gp = cross_val_predict(model, x, y, groups=groups, cv=grouped)
    overlaps = []
    for tr, te in grouped.split(x, y, groups):
        overlap = len(set(groups[tr]) & set(groups[te]))
        assert overlap == 0
        overlaps.append(overlap)
    return {
        "n_donors": n,
        "technical_rows": len(y),
        "sample_split_balanced_accuracy": balanced_accuracy_score(y, rp),
        "donor_split_balanced_accuracy": balanced_accuracy_score(y, gp),
        "grouped_overlap_counts": overlaps,
        "data_generating_truth": "No population-level feature/label relation; repeated donor fingerprints can leak donor-specific random labels.",
    }


def run(output=ROOT / "artifacts/run"):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid")
    cells = pd.DataFrame(
        [
            {
                "condition": "Before",
                "total_viable_cells": 1_000_000,
                "target_events": 20_000,
                "parent_events": 100_000,
            },
            {
                "condition": "After",
                "total_viable_cells": 400_000,
                "target_events": 35_000,
                "parent_events": 100_000,
            },
        ]
    )
    cells["target_fraction"] = cells.target_events / cells.parent_events
    cells["absolute_target_cells"] = [
        target_count(r.total_viable_cells, r.target_events, r.parent_events)
        for r in cells.itertuples()
    ]
    cells.to_csv(output / "cell_counts.csv", index=False)
    fig, axs = plt.subplots(1, 2, figsize=(9, 4))
    sns.barplot(
        data=cells, x="condition", y="target_fraction", ax=axs[0], color="#2b8cbe"
    )
    sns.barplot(
        data=cells, x="condition", y="absolute_target_cells", ax=axs[1], color="#238b45"
    )
    axs[0].set(title="Frequency rises", ylabel="Fraction of parent events")
    axs[1].set(title="Absolute count falls", ylabel="Estimated target cells")
    fig.suptitle("Synthetic example: the denominator changes the conclusion")
    fig.tight_layout()
    fig.savefig(output / "counts_vs_percentages.png", dpi=160)
    plt.close(fig)
    pairs = make_pairs()
    pairs.to_csv(output / "technical_measurements.csv", index=False)
    correct, donor = paired_summary(pairs)
    donor.to_csv(output / "donor_means.csv")
    wells = pairs.pivot(
        index=["donor", "replicate"], columns="condition", values="value"
    )
    wrong = stats.ttest_rel(wells.treated, wells.control)
    wrong_ci = wrong.confidence_interval()
    naive = {
        "n_pretended_independent_pairs": len(wells),
        "p": float(wrong.pvalue),
        "ci95": [float(wrong_ci.low), float(wrong_ci.high)],
    }
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot([0, 1], donor[["control", "treated"]].T, alpha=0.6)
    ax.set(
        xticks=[0, 1],
        xticklabels=["Control", "Treated"],
        ylabel="Synthetic assay response (arbitrary units)",
        title="Eight biological donors; twelve technical pairs each",
    )
    fig.tight_layout()
    fig.savefig(output / "donor_pairs.png", dpi=160)
    plt.close(fig)
    leak = leakage_case()
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(
        ["Split technical rows", "Hold out donors"],
        [leak["sample_split_balanced_accuracy"], leak["donor_split_balanced_accuracy"]],
        color=["#d95f02", "#1b9e77"],
    )
    ax.axhline(0.5, color="black", ls="--", label="Null balanced-accuracy reference")
    ax.set(
        ylim=(0, 1.05),
        ylabel="Balanced accuracy",
        title="Synthetic null-signal leakage control",
    )
    ax.legend()
    fig.tight_layout()
    fig.savefig(output / "leakage.png", dpi=160)
    plt.close(fig)
    result = {
        "source": "Entirely synthetic, seeded teaching examples. No actual patient or product data.",
        "cell_count_case": cells.to_dict("records"),
        "correct_donor_analysis": correct,
        "deliberately_wrong_technical_pair_analysis": naive,
        "leakage_negative_control": leak,
        "assumptions": [
            "Absolute target estimate assumes representative sampling and an independently valid total viable count from the same sample.",
            "Paired t interval assumes independent donors and an adequate model for donor differences; exact sign flips assume symmetric/exchangeable signs under the null.",
            "One fixed null simulation is a negative control, not a population generalization guarantee.",
        ],
    }
    (output / "summary.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
