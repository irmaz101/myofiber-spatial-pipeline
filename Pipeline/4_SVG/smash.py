"""Run SMASH on one AnnData file."""

import argparse
import sys
from pathlib import Path

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
    parser.add_argument("--fdr", type=float, default=0.01)
    return parser.parse_args()


def main():
    args = parse_args()
    sys.path.insert(0, str(args.source.resolve()))
    import SMASH

    adata = sc.read_h5ad(args.input)
    expression = adata.X
    if sparse.issparse(expression):
        expression = expression.toarray()

    expression = pd.DataFrame(
        expression,
        index=adata.obs_names,
        columns=adata.var_names,
    )
    coordinates = adata.obs[["x_centroid", "y_centroid"]]

    results = SMASH.SMASH(
        expression,
        coordinates,
        mean_only=False,
        kernel_covariance="All",
        forcePD=False,
    )["SMASH"].copy()

    results["adjusted_p_value"] = multipletests(
        results["p-val"],
        method="fdr_by",
    )[1]
    if "gene" not in results.columns:
        results.insert(0, "gene", results.index)
    results["is_svg"] = results["adjusted_p_value"] < args.fdr

    if args.all_output:
        args.all_output.parent.mkdir(parents=True, exist_ok=True)
        results.to_csv(args.all_output, index=False)

    selected = results[results["is_svg"]]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(args.output, index=False)
    print(f"Saved {len(selected)} SVGs to {args.output}")


if __name__ == "__main__":
    main()