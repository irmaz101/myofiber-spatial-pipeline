"""Run Hotspot on one AnnData file."""

import argparse
import os
from pathlib import Path

import hotspot
import pandas as pd
import scanpy as sc


def parse_args():
    parser = argparse.ArgumentParser(description="Run Hotspot")
    parser.add_argument("-i", "--input", type=Path, required=True)
    parser.add_argument("-o", "--output", type=Path, required=True)
    parser.add_argument("--all-output", type=Path)
    parser.add_argument("--neighbors", type=int, default=20)
    parser.add_argument("--fdr", type=float, default=0.05)
    parser.add_argument("--jobs", type=int, default=os.cpu_count() or 1)
    return parser.parse_args()


def main():
    args = parse_args()
    source = sc.read_h5ad(args.input)
    counts = source.layers["counts"]

    adata = sc.AnnData(
        X=counts,
        obs=source.obs.copy(),
        var=pd.DataFrame(index=source.var_names),
    )
    adata.layers["counts"] = counts
    sc.pp.filter_genes(adata, min_cells=1)
    adata.obsm["spatial"] = source.obs[
        ["x_centroid", "y_centroid"]
    ].to_numpy()

    model = hotspot.Hotspot(
        adata,
        model="bernoulli",
        latent_obsm_key="spatial",
        layer_key="counts",
    )
    model.create_knn_graph(
        weighted_graph=False,
        n_neighbors=min(args.neighbors, adata.n_obs - 1),
    )

    results = model.compute_autocorrelations(jobs=args.jobs)
    results.index.name = "gene"
    results = results.sort_values("Z", ascending=False)
    results["is_svg"] = results["FDR"] < args.fdr

    if args.all_output:
        args.all_output.parent.mkdir(parents=True, exist_ok=True)
        results.to_csv(args.all_output)

    selected = results[results["is_svg"]]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(args.output)
    print(f"Saved {len(selected)} SVGs to {args.output}")


if __name__ == "__main__":
    main()