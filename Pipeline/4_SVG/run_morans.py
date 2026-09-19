"""Run Moran's I on one AnnData file."""

import argparse
from pathlib import Path

import scanpy as sc
import squidpy as sq


def parse_args():
    parser = argparse.ArgumentParser(description="Run Moran's I")
    parser.add_argument("-i", "--input", type=Path, required=True)
    parser.add_argument("-o", "--output", type=Path, required=True)
    parser.add_argument("--all-output", type=Path)
    parser.add_argument("--neighbors", type=int, default=6)
    parser.add_argument("--fdr", type=float, default=0.05)
    return parser.parse_args()


def main():
    args = parse_args()

    original = sc.read_h5ad(args.input)

    if original.raw is None:
        raise ValueError(
            f"{args.input.name}: raw counts are missing from adata.raw."
        )

    required_columns = {"x_centroid", "y_centroid"}

    if not required_columns.issubset(original.obs.columns):
        raise ValueError(
            f"{args.input.name}: x_centroid and y_centroid are missing."
        )

    # Reconstruct AnnData from raw counts.
    adata = original.raw.to_adata()
    adata.obs = original.obs.copy()

    adata.obsm["spatial"] = adata.obs[
        ["x_centroid", "y_centroid"]
    ].to_numpy()

    # Normalize and log-transform raw counts.
    sc.pp.normalize_total(
        adata,
        target_sum=10_000,
    )
    sc.pp.log1p(adata)

    sq.gr.spatial_neighbors(
        adata,
        coord_type="generic",
        n_neighs=args.neighbors,
    )

    sq.gr.spatial_autocorr(
        adata,
        genes=adata.var_names.tolist(),
        mode="moran",
        attr="X",
    )

    results = adata.uns["moranI"].copy()
    results.index.name = "gene"

    required_result_columns = {
        "I",
        "pval_norm_fdr_bh",
    }

    missing_columns = required_result_columns.difference(
        results.columns
    )

    if missing_columns:
        raise ValueError(
            f"{args.input.name}: missing result columns: "
            f"{sorted(missing_columns)}"
        )

    results["is_svg"] = (
        results["pval_norm_fdr_bh"] < args.fdr
    )

    if args.all_output:
        args.all_output.parent.mkdir(
            parents=True,
            exist_ok=True,
        )
        results.to_csv(args.all_output)

    selected = (
        results[
            results["is_svg"]
        ]
        .sort_values(
            "I",
            ascending=False,
        )
    )

    args.output.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    selected.to_csv(args.output)

    print(
        f"Saved {len(selected)} SVGs "
        f"to {args.output}"
    )


if __name__ == "__main__":
    main()
