"""Run SOMDE on one AnnData file."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from somde import SomNode


def parse_args():
    parser = argparse.ArgumentParser(description="Run SOMDE")
    parser.add_argument("-i", "--input", type=Path, required=True)
    parser.add_argument("-o", "--output", type=Path, required=True)
    parser.add_argument("--all-output", type=Path)
    parser.add_argument("--grid-size", type=int, default=10)
    parser.add_argument("--fdr", type=float, default=0.05)
    return parser.parse_args()


def main():
    args = parse_args()
    adata = sc.read_h5ad(args.input)

    counts = adata.layers["counts"]
    if sparse.issparse(counts):
        counts = counts.toarray()

    expressed = np.asarray(counts).sum(axis=0) > 0
    expression = pd.DataFrame(
        np.asarray(counts)[:, expressed].T,
        index=adata.var_names[expressed],
        columns=adata.obs_names,
    )
    coordinates = adata.obs[
        ["x_centroid", "y_centroid"]
    ].to_numpy(dtype=np.float32)

    som = SomNode(coordinates, k=args.grid_size)
    som.mtx(expression)
    som.norm()
    results, _ = som.run()

    results = results.sort_values("LLR", ascending=False)
    results = results.rename(columns={"g": "gene"})
    results["is_svg"] = results["qval"] < args.fdr

    if args.all_output:
        args.all_output.parent.mkdir(parents=True, exist_ok=True)
        results.to_csv(args.all_output, index=False)

    selected = results[results["is_svg"]]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    selected.to_csv(args.output, index=False)
    print(f"Saved {len(selected)} SVGs to {args.output}")


if __name__ == "__main__":
    main()
