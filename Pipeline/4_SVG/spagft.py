"""Run SpaGFT on one AnnData file."""

import argparse
from pathlib import Path

import scanpy as sc
import SpaGFT as sg


def parse_args():
    parser = argparse.ArgumentParser(description="Run SpaGFT")
    parser.add_argument("-i", "--input", type=Path, required=True)
    parser.add_argument("-o", "--output", type=Path, required=True)
    parser.add_argument("--all-output", type=Path)
    parser.add_argument("--fdr", type=float, default=0.05)
    return parser.parse_args()


def main():
    args = parse_args()
    adata = sc.read_h5ad(args.input)
    adata.var_names_make_unique()

    results = sg.detect_svg(
        adata,
        spatial_info=["x_centroid", "y_centroid"],
        ratio_neighbors=1,
        filter_peaks=True,
        S=6,
    )
    results["is_svg"] = (
        results["cutoff_gft_score"]
        & (results["fdr"] < args.fdr)
    )
    results.index.name = "gene"

    if args.all_output:
        args.all_output.parent.mkdir(parents=True, exist_ok=True)
        results.to_csv(args.all_output)

    selected = results[results["is_svg"]]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(args.output)
    print(f"Saved {len(selected)} SVGs to {args.output}")


if __name__ == "__main__":
    main()