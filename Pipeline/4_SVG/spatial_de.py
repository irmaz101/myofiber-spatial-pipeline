"""Run SpatialDE on one AnnData file."""

import argparse
from pathlib import Path

import NaiveDE
import numpy as np
import pandas as pd
import scanpy as sc
import SpatialDE
from scipy import sparse


def parse_args():
    parser = argparse.ArgumentParser(description="Run SpatialDE")
    parser.add_argument("-i", "--input", type=Path, required=True)
    parser.add_argument("-o", "--output", type=Path, required=True)
    parser.add_argument("--all-output", type=Path)
    parser.add_argument("--qvalue", type=float, default=0.05)
    return parser.parse_args()


def main():
    args = parse_args()
    adata = sc.read_h5ad(args.input)
    adata.var_names_make_unique()

    counts = adata.layers["counts"]
    if sparse.issparse(counts):
        counts = counts.toarray()

    counts = pd.DataFrame(
        counts,
        index=adata.obs_names,
        columns=adata.var_names,
    )
    counts = counts.loc[:, counts.sum() > 0]

    sample_info = pd.DataFrame({
        "x": adata.obs["x_centroid"].to_numpy(),
        "y": adata.obs["y_centroid"].to_numpy(),
        "total_counts": np.asarray(counts.sum(axis=1)),
    })

    normalized = NaiveDE.stabilize(counts.T).T
    residuals = NaiveDE.regress_out(
        sample_info,
        normalized.T,
        "np.log(total_counts)",
    ).T

    results = SpatialDE.run(
        sample_info[["x", "y"]],
        residuals,
    )
    results = (
        results.sort_values("LLR", ascending=False)
        .drop_duplicates("g")
    )
    results["is_svg"] = results["qval"] < args.qvalue

    if args.all_output:
        args.all_output.parent.mkdir(parents=True, exist_ok=True)
        results.to_csv(args.all_output, index=False)

    selected = results[results["is_svg"]]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(args.output, index=False)
    print(f"Saved {len(selected)} SVGs to {args.output}")


if __name__ == "__main__":
    main()
