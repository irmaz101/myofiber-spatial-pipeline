"""Run SpaGFT on one AnnData file."""

import argparse
from pathlib import Path

import numpy as np
import scanpy as sc
import SpaGFT as sg
from scipy.sparse import issparse


def parse_args():
    parser = argparse.ArgumentParser(description="Run SpaGFT")
    parser.add_argument("-i", "--input", type=Path, required=True)
    parser.add_argument("-o", "--output", type=Path, required=True)
    parser.add_argument("--all-output", type=Path)
    parser.add_argument("--fdr", type=float, default=0.05)
    return parser.parse_args()


def main():
    args = parse_args()
    original = sc.read_h5ad(args.input)

    if original.raw is None:
        raise ValueError("adata.raw with raw counts is missing")

    if not {"x_centroid", "y_centroid"}.issubset(original.obs.columns):
        raise ValueError("x_centroid and y_centroid are missing")

    # Reconstruct normalized and log-transformed expression from raw counts.
    adata = original.raw.to_adata()
    adata.obs = original.obs.copy()

    if issparse(adata.X):
        minimum = adata.X.data.min() if adata.X.data.size else 0
    else:
        minimum = np.asarray(adata.X).min()

    if minimum < 0:
        raise ValueError("adata.raw.X contains negative values")

    adata.uns.pop("log1p", None)

    sc.pp.normalize_total(
        adata,
        target_sum=10_000,
    )
    sc.pp.log1p(adata)

    adata.var_names_make_unique()

    results = sg.detect_svg(
        adata,
        spatial_info=["x_centroid", "y_centroid"],
        ratio_neighbors=1,
        filter_peaks=True,
        S=6,
    )

    results = (
        results
        .rename_axis("gene")
        .reset_index()
    )

    results["is_svg"] = (
        results["cutoff_gft_score"]
        & (results["fdr"] <= args.fdr)
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
            "gft_score",
            ascending=False,
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
