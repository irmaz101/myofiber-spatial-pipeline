"""
Calculate pairwise overlap between SVG methods.

For each pair of methods, the script compares SVG gene sets from matching
samples, calculates the overlap coefficient, saves per-sample and summary
tables, and generates a heatmap of median overlap.
"""

from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# --------------------------------------------------
# Settings
# --------------------------------------------------

results_dir = Path(".")
output_dir = Path("summary")
replicate_file = Path("technical_replicates.csv")

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

gene_column_candidates = [
    "gene",
    "gene_component",
    "feature",
    "feature_name",
    "symbol",
    "Gene",
    "genes",
]

output_dir.mkdir(parents=True, exist_ok=True)

excluded_samples = set()

if replicate_file.exists():
    replicate_pairs = pd.read_csv(replicate_file)
    required_columns = {"sample_1", "sample_2"}
    missing_columns = required_columns.difference(replicate_pairs.columns)

    if missing_columns:
        raise ValueError(
            f"{replicate_file} is missing columns: "
            f"{sorted(missing_columns)}"
        )

    excluded_samples = set(
        replicate_pairs["sample_2"]
        .dropna()
        .astype(str)
    )


# --------------------------------------------------
# Helper functions
# --------------------------------------------------

def find_gene_column(dataframe):
    """Identify the column containing gene identifiers."""

    if dataframe.shape[1] == 1:
        return dataframe.columns[0]

    for column in gene_column_candidates:
        if column in dataframe.columns:
            return column

    raise ValueError(
        f"Gene column not found. Available columns: "
        f"{list(dataframe.columns)}"
    )


def load_gene_set(csv_path):
    """Load unique gene identifiers from one result file."""

    dataframe = pd.read_csv(csv_path)
    gene_column = find_gene_column(dataframe)

    return set(
        dataframe[gene_column]
        .dropna()
        .astype(str)
        .str.strip()
    )


def overlap_coefficient(set_1, set_2):
    """Calculate intersection size divided by the smaller set size."""

    if not set_1 or not set_2:
        return np.nan

    return len(set_1 & set_2) / min(
        len(set_1),
        len(set_2),
    )


# --------------------------------------------------
# Load SVG gene sets
# --------------------------------------------------

svg_sets = {}

for method in methods:
    method_dir = results_dir / method

    if not method_dir.is_dir():
        print(f"Skipping missing directory: {method_dir}")
        continue

    svg_sets[method] = {
        csv_path.stem: load_gene_set(csv_path)
        for csv_path in sorted(method_dir.glob("*.csv"))
        if csv_path.stem not in excluded_samples
    }

    if not svg_sets[method]:
        print(f"No CSV files found for {method}.")

available_methods = [
    method
    for method in methods
    if method in svg_sets and svg_sets[method]
]

if len(available_methods) < 2:
    raise ValueError(
        "Results from at least two SVG methods are required."
    )


# --------------------------------------------------
# Pairwise overlap
# --------------------------------------------------

records = []

for method_1, method_2 in combinations(
    available_methods,
    2,
):
    common_samples = sorted(
        set(svg_sets[method_1])
        & set(svg_sets[method_2])
    )

    if not common_samples:
        print(
            f"No matching samples between "
            f"{method_1} and {method_2}."
        )
        continue

    for sample in common_samples:
        genes_1 = svg_sets[method_1][sample]
        genes_2 = svg_sets[method_2][sample]

        records.append({
            "method_1": method_1,
            "method_2": method_2,
            "sample": sample,
            "n_svg_1": len(genes_1),
            "n_svg_2": len(genes_2),
            "intersection": len(genes_1 & genes_2),
            "overlap_coefficient": overlap_coefficient(
                genes_1,
                genes_2,
            ),
        })

overlap_results = pd.DataFrame(records)

if overlap_results.empty:
    raise ValueError(
        "No overlap comparisons were created. "
        "Check that sample filenames match across method directories."
    )

for column in ["method_1", "method_2"]:
    overlap_results[column] = pd.Categorical(
        overlap_results[column],
        categories=methods,
        ordered=True,
    )

overlap_results = overlap_results.sort_values(
    ["method_1", "method_2", "sample"]
)

overlap_results.to_csv(
    output_dir / "svg_overlap_per_sample.csv",
    index=False,
)


# --------------------------------------------------
# Summary statistics
# --------------------------------------------------

summary = (
    overlap_results
    .groupby(
        ["method_1", "method_2"],
        observed=True,
    )["overlap_coefficient"]
    .agg(
        median="median",
        mean="mean",
        std="std",
        n_samples="count",
    )
    .reset_index()
    .sort_values(["method_1", "method_2"])
)

summary.to_csv(
    output_dir / "svg_overlap_summary.csv",
    index=False,
)

print(summary)


# --------------------------------------------------
# Heatmap
# --------------------------------------------------

heatmap = pd.DataFrame(
    np.nan,
    index=available_methods,
    columns=available_methods,
)

for row in summary.itertuples(index=False):
    heatmap.loc[
        row.method_1,
        row.method_2,
    ] = row.median

    heatmap.loc[
        row.method_2,
        row.method_1,
    ] = row.median

for method in available_methods:
    heatmap.loc[method, method] = 1.0

plt.figure(figsize=(9, 7))

image = plt.imshow(
    heatmap,
    aspect="auto",
    vmin=0,
    vmax=1,
)

plt.colorbar(
    image,
    label="Median overlap coefficient",
)

plt.xticks(
    np.arange(len(heatmap.columns)),
    heatmap.columns,
    rotation=45,
    ha="right",
)

plt.yticks(
    np.arange(len(heatmap.index)),
    heatmap.index,
)

plt.tight_layout()

plt.savefig(
    output_dir / "svg_overlap_heatmap.png",
    dpi=300,
    bbox_inches="tight",
)

plt.close()
