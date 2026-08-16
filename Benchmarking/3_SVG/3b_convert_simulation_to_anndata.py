"""Convert one scDesign3 simulation from CSV files to AnnData."""

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert one scDesign3 simulation to AnnData"
    )
    parser.add_argument("--counts", type=Path, required=True)
    parser.add_argument("--locations", type=Path, required=True)
    parser.add_argument("-o", "--output", type=Path, required=True)
    parser.add_argument("--target-sum", type=float, default=10_000)
    return parser.parse_args()


def main():
    args = parse_args()

    locations = pd.read_csv(args.locations, index_col=0)
    counts = pd.read_csv(args.counts, index_col=0).transpose()

    required_coordinates = {"spatial1", "spatial2"}
    missing_coordinates = required_coordinates.difference(locations.columns)

    if missing_coordinates:
        raise ValueError(
            f"{args.locations} is missing columns: "
            f"{sorted(missing_coordinates)}"
        )

    if counts.index.equals(locations.index):
        pass
    elif set(counts.index) == set(locations.index):
        locations = locations.loc[counts.index]
    else:
        raise ValueError(
            "Cell identifiers differ between the count and location files."
        )

    feature_metadata = pd.DataFrame(index=counts.columns)
    feature_metadata["feature_name"] = feature_metadata.index.astype(str)

    parsed_features = feature_metadata["feature_name"].str.rsplit(
        "_",
        n=1,
        expand=True,
    )

    if parsed_features.shape[1] != 2:
        raise ValueError(
            "Feature identifiers must end with '_<spatial_var>'."
        )

    feature_metadata["gene"] = parsed_features.iloc[:, 0]
    feature_metadata["spatial_var"] = pd.to_numeric(
        parsed_features.iloc[:, 1],
        errors="raise",
    )

    observation_metadata = locations.copy()
    observation_metadata["x_centroid"] = observation_metadata["spatial1"]
    observation_metadata["y_centroid"] = observation_metadata["spatial2"]

    adata = ad.AnnData(
        X=sparse.csr_matrix(counts.to_numpy(dtype=np.float32)),
        obs=observation_metadata,
        var=feature_metadata,
    )

    adata.obsm["spatial"] = observation_metadata[
        ["x_centroid", "y_centroid"]
    ].to_numpy()

    # Preserve raw counts for all method-specific SVG wrappers.
    adata.layers["counts"] = adata.X.copy()
    adata.raw = adata.copy()

    sc.pp.calculate_qc_metrics(
        adata,
        percent_top=[10],
        inplace=True,
    )
    sc.pp.normalize_total(
        adata,
        target_sum=args.target_sum,
    )
    sc.pp.log1p(adata)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.output)

    print(f"Saved {adata.n_obs} objects and {adata.n_vars} features")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
