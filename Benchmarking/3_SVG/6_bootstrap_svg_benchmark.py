
"""
Bootstrap SVG benchmark results across synthetic datasets.

For each metric, datasets are sampled with replacement and the mean score
is calculated for each method. The best and second-best methods are
recorded for every bootstrap iteration.
"""

from pathlib import Path

import numpy as np
import pandas as pd


# --------------------------------------------------
# Settings
# --------------------------------------------------

results_dir = Path("results")
output_dir = Path("bootstrap_results")

n_bootstrap = 100
random_seed = 42

metrics = {
    "kendall": {
        "file": results_dir / "kendall_results.csv",
        "value_column": "correlation",
    },
    "auprc": {
        "file": results_dir / "aupr.csv",
        "value_column": "aupr",
    },
}

output_dir.mkdir(parents=True, exist_ok=True)


def bootstrap_metric(
    input_file,
    value_column,
    metric_name,
):
    """Bootstrap method rankings for one benchmark metric."""

    data = pd.read_csv(input_file)

    required_columns = {
        "dataset",
        "method",
        value_column,
    }

    missing = required_columns.difference(data.columns)

    if missing:
        raise ValueError(
            f"{input_file} is missing columns: {sorted(missing)}"
        )

    scores = data.pivot_table(
        index="dataset",
        columns="method",
        values=value_column,
        aggfunc="mean",
    )

    if scores.shape[1] < 2:
        raise ValueError(
            f"At least two methods are required for {metric_name}."
        )

    rng = np.random.default_rng(random_seed)
    datasets = scores.index.to_numpy()
    records = []

    for iteration in range(1, n_bootstrap + 1):
        sampled_datasets = rng.choice(
            datasets,
            size=len(datasets),
            replace=True,
        )

        mean_scores = (
            scores.loc[sampled_datasets]
            .mean()
            .sort_values(ascending=False)
        )

        record = {
            "iteration": iteration,
            "winner": mean_scores.index[0],
            "winner_score": mean_scores.iloc[0],
            "second": mean_scores.index[1],
            "second_score": mean_scores.iloc[1],
            "margin": mean_scores.iloc[0] - mean_scores.iloc[1],
        }

        record.update(mean_scores.to_dict())
        records.append(record)

    bootstrap_results = pd.DataFrame(records)

    winner_summary = (
        bootstrap_results["winner"]
        .value_counts()
        .rename_axis("method")
        .reset_index(name="wins")
    )

    winner_summary["proportion"] = (
        winner_summary["wins"] / n_bootstrap
    )

    bootstrap_results.to_csv(
        output_dir / f"{metric_name}_bootstrap_results.csv",
        index=False,
    )

    winner_summary.to_csv(
        output_dir / f"{metric_name}_winner_summary.csv",
        index=False,
    )

    print(f"\n{metric_name.upper()}")
    print(winner_summary)


def main():
    for metric_name, settings in metrics.items():
        bootstrap_metric(
            input_file=settings["file"],
            value_column=settings["value_column"],
            metric_name=metric_name,
        )


if __name__ == "__main__":
    main()