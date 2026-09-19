"""Run SMASH on one AnnData file using normalized raw counts."""

import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from statsmodels.stats.multitest import multipletests


def parse_args():
    parser = argparse.ArgumentParser(description="Run SMASH")
    parser.add_argument("-i", "--input", type=Path, required=True)
    parser.add_argument("-o", "--output", type=Path, required=True)
    parser.add_argument("--all-output", type=Path)
    parser.add_argument("--source", type=Path, default=Path("SMASH/SMASH"))
    parser.add_argument("--target-sum", type=float, default=10_000)
    parser.add_argument("--fdr", type=float, default=0.05)
    return parser.parse_args()


def main():
    args = parse_args()
    sys.path.insert(0, str(args.source.resolve()))
    import SMASH

    sample = args.input.stem
    original = sc.read_h5ad(args.input)

    if original.raw is None:
        raise ValueError(
            f"{sample}: adata.raw with raw counts is missing."
        )

    if not {"x_centroid", "y_centroid"}.issubset(
        original.obs.columns
    ):
        raise ValueError(
            f"{sample}: x_centroid and y_centroid are missing."
        )

    # Reconstruct expression from raw counts.
    adata = original.raw.to_adata()
    adata.obs = original.obs.copy()

    if sparse.issparse(adata.X):
        stored_values = adata.X.data
        minimum = (
            stored_values.min()
            if stored_values.size
            else 0
        )
    else:
        minimum = np.asarray(adata.X).min()

    if minimum < 0:
        raise ValueError(
            f"{sample}: adata.raw.X contains negative values."
        )

    sc.pp.normalize_total(
        adata,
        target_sum=args.target_sum,
    )
    sc.pp.log1p(adata)

    expression = adata.X

    if sparse.issparse(expression):
        expression = expression.toarray()
    else:
        expression = np.asarray(expression)

    expression = pd.DataFrame(
        expression,
        index=adata.obs_names,
        columns=adata.var_names,
    )

    coordinates = original.obs[
        ["x_centroid", "y_centroid"]
    ].copy()
    coordinates.columns = ["x", "y"]

    run_result = SMASH.SMASH(
        expression,
        coordinates,
        mean_only=False,
        kernel_covariance="All",
        forcePD=False,
    )

    results = run_result["SMASH"].copy()

    if "gene" not in results.columns:
        results = (
            results
            .rename_axis("gene")
            .reset_index()
        )
    else:
        results = results.reset_index(drop=True)

    results["Adjusted p-val"] = multipletests(
        results["p-val"].to_numpy(),
        method="fdr_by",
    )[1]

    results["is_svg"] = (
        results["Adjusted p-val"] <= args.fdr
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

    selected = (
        results[
            results["is_svg"]
        ]
        .sort_values(
            "Adjusted p-val",
            ascending=True,
        )
    )

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
