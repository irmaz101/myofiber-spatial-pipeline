"""Run scGCO on one AnnData file."""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


def parse_args():
    parser = argparse.ArgumentParser(description="Run scGCO")
    parser.add_argument("-i", "--input", type=Path, required=True)
    parser.add_argument("-o", "--output", type=Path, required=True)
    parser.add_argument("--all-output", type=Path)
    parser.add_argument("--source", type=Path, default=Path("scgco_utils"))
    parser.add_argument("--min-objects", type=int, default=2)
    parser.add_argument("--fdr", type=float, default=0.05)
    return parser.parse_args()


def main():
    args = parse_args()
    sys.path.insert(0, str(args.source.resolve()))
    from scGCO_simple import (
        create_graph_with_weight,
        gmm_model,
        identify_spatial_genes,
        normalize_count_cellranger,
    )

    adata = sc.read_h5ad(args.input)
    counts = adata.layers["counts"]
    if sparse.issparse(counts):
        counts = counts.toarray()

    expression = pd.DataFrame(
        np.asarray(counts),
        index=adata.obs_names,
        columns=adata.var_names,
    )
    expression = expression.loc[
        :,
        (expression > 0).sum() >= args.min_objects,
    ]
    expression = normalize_count_cellranger(expression)
    expression = expression.loc[:, expression.var() > 0]

    coordinates = adata.obs[
        ["x_centroid", "y_centroid"]
    ].to_numpy()
    graph = create_graph_with_weight(
        coordinates,
        expression.iloc[:, 0],
    )
    models = gmm_model(expression)
    results = identify_spatial_genes(
        coordinates,
        expression,
        graph,
        models,
    )

    results["fdr"] = pd.to_numeric(results["fdr"], errors="coerce")
    results.insert(0, "gene_component", results.index)
    results["is_svg"] = results["fdr"] < args.fdr

    if args.all_output:
        args.all_output.parent.mkdir(parents=True, exist_ok=True)
        results.to_csv(args.all_output, index=False)

    selected = results[results["is_svg"]]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(args.output, index=False)
    print(f"Saved {len(selected)} SVGs to {args.output}")


if __name__ == "__main__":
    main()
