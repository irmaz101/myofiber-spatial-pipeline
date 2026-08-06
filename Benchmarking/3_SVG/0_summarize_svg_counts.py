"""
Summarize the number of SVGs detected by each method.

The script counts result rows per sample, excludes designated technical
replicates, saves per-sample and per-method summaries, and creates a boxplot.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


# --------------------------------------------------
# Settings
# --------------------------------------------------

results_dir = Path(".")
output_dir = Path("summary")
output_dir.mkdir(parents=True, exist_ok=True)

methods = [
    "hotspot",
    "moran_i",
    "scgco",
    "smash",
    "somde",
    "spagft",
    "spatialde",
    "svgbit",
]

# --------------------------------------------------
# Collect SVG counts
# --------------------------------------------------

records = []

for method in methods:
    method_dir = results_dir / method

    if not method_dir.is_dir():
        print(f"Skipping missing directory: {method_dir}")
        continue

    for csv_path in sorted(method_dir.glob("*.csv")):
        sample = csv_path.stem

        n_svgs = len(pd.read_csv(csv_path))

        records.append({
            "method": method,
            "sample": sample,
            "n_svgs": n_svgs,
        })

counts = pd.DataFrame(records)

if counts.empty:
    raise ValueError("No SVG result files were found.")

counts["method"] = pd.Categorical(
    counts["method"],
    categories=methods,
    ordered=True,
)

counts = counts.sort_values(["method", "sample"])
counts.to_csv(
    output_dir / "svg_counts_per_sample.csv",
    index=False,
)

# --------------------------------------------------
# Summary statistics
# --------------------------------------------------
summary = (
    counts
    .groupby("method", observed=True)["n_svgs"]
    .agg(
        median="median",
        q1=lambda x: x.quantile(0.25),
        q3=lambda x: x.quantile(0.75),
        mean="mean",
        std="std",
        minimum="min",
        maximum="max",
        n_samples="count",
    )
    .reset_index()
)

summary.to_csv(
    output_dir / "svg_counts_summary.csv",
    index=False,
)

print(summary)

# --------------------------------------------------
# Plot
# --------------------------------------------------
available_methods = [
    method
    for method in methods
    if method in counts["method"].astype(str).unique()
]

plot_data = [
    counts.loc[counts["method"] == method, "n_svgs"]
    for method in available_methods
]

plt.figure(figsize=(10, 5))

plt.boxplot(
    plot_data,
    tick_labels=available_methods,
    showfliers=False,
    widths=0.6,
)

plt.ylabel("Number of SVGs")
plt.xlabel("Method")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()

plt.savefig(
    output_dir / "svg_counts_boxplot.png",
    dpi=300,
    bbox_inches="tight",
)

plt.close()