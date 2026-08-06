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
    adata = sc.read_h5ad(args.input)

    adata.obsm["spatial"] = adata.obs[
        ["x_centroid", "y_centroid"]
    ].to_numpy()

    sq.gr.spatial_neighbors(
        adata,
        coord_type="generic",
        n_neighs=args.neighbors,
    )
    sq.gr.spatial_autocorr(adata, mode="moran")

    results = adata.uns["moranI"].copy()
    results.index.name = "gene"
    results["is_svg"] = results["pval_norm_fdr_bh"] < args.fdr

    if args.all_output:
        args.all_output.parent.mkdir(parents=True, exist_ok=True)
        results.to_csv(args.all_output)

    selected = results[results["is_svg"]]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(args.output)
    print(f"Saved {len(selected)} SVGs to {args.output}")


if __name__ == "__main__":
    main()