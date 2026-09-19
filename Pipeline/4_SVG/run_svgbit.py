"""Run SVGbit on one AnnData file."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import svgbit as sb
from scipy import sparse


def parse_args():
    parser = argparse.ArgumentParser(description="Run SVGbit")
    parser.add_argument("-i", "--input", type=Path, required=True)
    parser.add_argument("-o", "--output", type=Path, required=True)
    parser.add_argument("--all-output", type=Path)
    parser.add_argument("--neighbors", type=int, default=6)
    parser.add_argument("--top", type=int, default=1000)
    return parser.parse_args()


def main():
    args = parse_args()
    adata = sc.read_h5ad(args.input)

    if adata.raw is None:
        raise ValueError("adata.raw with raw counts is missing")

    expression = adata.raw.X

    if sparse.issparse(expression):
        expression = expression.toarray()

    expression = np.asarray(expression)

    # Remove genes with zero counts or zero variance.
    keep = (
        (expression.sum(axis=0) > 0)
        & (expression.var(axis=0) > 0)
    )

    expression = pd.DataFrame(
        expression[:, keep],
        index=adata.obs_names,
        columns=adata.raw.var_names[keep],
    )

    coordinates = pd.DataFrame(
        adata.obs[
            ["x_centroid", "y_centroid"]
        ].to_numpy(),
        index=adata.obs_names,
        columns=["x", "y"],
    )

    dataset = sb.STDataset(
        count_df=expression,
        coordinate_df=coordinates,
    )

    dataset = sb.normalizers.logcpm_normalizer(dataset)
    dataset = sb.filters.low_variance_filter(dataset, var=0)

    dataset.acquire_weight(
        k=min(args.neighbors, adata.n_obs - 1)
    )
    dataset.acquire_hotspot()
    dataset.acquire_density()
    dataset.find_clusters()

    results = pd.DataFrame({
        "gene": dataset.genes,
        "AI": np.asarray(dataset.AI).ravel(),
        "Di_mean": np.asarray(
            dataset.Di.mean(axis=0)
        ).ravel(),
    })

    results = (
        results
        .sort_values("AI", ascending=False)
        .reset_index(drop=True)
    )

    results["rank_ai"] = np.arange(
        1,
        len(results) + 1,
    )
    results["is_svg"] = (
        results["rank_ai"] <= args.top
    )

    if args.all_output:
        args.all_output.parent.mkdir(
            parents=True,
            exist_ok=True,
        )
        results.to_csv(
            args.all_output,
            index=False,
        )

    selected = results[
        results["is_svg"]
    ]

    args.output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )
    selected.to_csv(
        args.output,
        index=False,
    )

    print(
        f"Saved {len(selected)} SVGs "
        f"to {args.output}"
    )


if __name__ == "__main__":
    main()
