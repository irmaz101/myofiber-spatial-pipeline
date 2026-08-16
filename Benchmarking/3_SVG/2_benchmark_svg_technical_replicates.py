"""
Compare SVG sets across technical-replicate pairs.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


results_dir = Path(".")
output_dir = Path("technical_replicate_summary")
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

gene_columns = [
    "gene",
    "gene_name",
    "gene_component",
    "feature",
    "feature_name",
    "symbol",
    "Gene",
    "genes",
]

output_dir.mkdir(exist_ok=True)


def load_genes(path):
    df = pd.read_csv(path)

    if df.shape[1] == 1:
        column = df.columns[0]
    else:
        column = next(
            (c for c in gene_columns if c in df.columns),
            None,
        )

    if column is None:
        raise ValueError(f"No gene column found in {path}")

    return set(
        df[column]
        .dropna()
        .astype(str)
        .str.strip()
    )


def overlap(a, b):
    return np.nan if not a or not b else len(a & b) / min(len(a), len(b))


def jaccard(a, b):
    return np.nan if not (a or b) else len(a & b) / len(a | b)


pairs = pd.read_csv(replicate_file)
records = []

for method in methods:
    for row in pairs.itertuples(index=False):
        file_1 = results_dir / method / f"{row.sample_1}.csv"
        file_2 = results_dir / method / f"{row.sample_2}.csv"

        if not file_1.exists() or not file_2.exists():
            continue

        genes_1 = load_genes(file_1)
        genes_2 = load_genes(file_2)

        records.append({
            "method": method,
            "sample_1": row.sample_1,
            "sample_2": row.sample_2,
            "n_svg_1": len(genes_1),
            "n_svg_2": len(genes_2),
            "intersection": len(genes_1 & genes_2),
            "overlap_coefficient": overlap(genes_1, genes_2),
            "jaccard_index": jaccard(genes_1, genes_2),
        })

results = pd.DataFrame(records)

if results.empty:
    raise ValueError("No valid replicate comparisons were found.")

results.to_csv(
    output_dir / "svg_replicate_metrics.csv",
    index=False,
)

summary = (
    results.groupby("method", as_index=False)
    .agg(
        n_pairs=("method", "size"),
        overlap_median=("overlap_coefficient", "median"),
        jaccard_median=("jaccard_index", "median"),
    )
)

summary.to_csv(
    output_dir / "svg_replicate_summary.csv",
    index=False,
)

results["pair"] = results["sample_1"] + " vs " + results["sample_2"]

for metric in ["overlap_coefficient", "jaccard_index"]:
    heatmap = results.pivot(
        index="method",
        columns="pair",
        values=metric,
    ).reindex(methods)

    plt.figure(figsize=(8, 5))
    image = plt.imshow(
        heatmap,
        aspect="auto",
        vmin=0,
        vmax=1,
    )

    plt.colorbar(image, label=metric.replace("_", " ").title())
    plt.xticks(
        range(len(heatmap.columns)),
        heatmap.columns,
        rotation=45,
        ha="right",
    )
    plt.yticks(
        range(len(heatmap.index)),
        heatmap.index,
    )
    plt.tight_layout()
    plt.savefig(
        output_dir / f"{metric}_heatmap.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()