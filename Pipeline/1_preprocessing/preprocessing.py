"""
Assign transcripts to segmented objects, compare count- and density-based
QC filters, apply transcript-density filtering, and preprocess each sample.
"""

import os

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import tifffile
from scipy import ndimage


def get_file_pairs(transcript_dir="transcripts", mask_dir="masks_mf"):
    """Match transcript CSV files and mask TIFF files by sample name."""

    pairs = []

    for csv_file in sorted(os.listdir(transcript_dir)):
        if not csv_file.endswith(".csv"):
            continue

        sample = os.path.splitext(csv_file)[0]
        mask_path = os.path.join(mask_dir, f"{sample}.tif")

        if os.path.exists(mask_path):
            pairs.append(
                (
                    os.path.join(transcript_dir, csv_file),
                    mask_path,
                )
            )
        else:
            print(f"WARNING: No matching mask found for {csv_file}")

    return pairs


def preprocess_sample(
    csv_path,
    mask_path,
    output_dir="anndata_mf_density",
):
    """Create and preprocess one AnnData object."""

    transcripts = pd.read_csv(csv_path)
    mask = tifffile.imread(mask_path)

    required_columns = {
        "feature_name",
        "x_location",
        "y_location",
    }

    missing_columns = required_columns.difference(
        transcripts.columns
    )

    if missing_columns:
        raise ValueError(
            f"{csv_path} is missing columns: "
            f"{sorted(missing_columns)}"
        )

    if mask.ndim != 2:
        raise ValueError(
            f"Expected a 2D mask, found shape {mask.shape}"
        )

    # Convert transcript coordinates to integer indices.
    transcripts["x_location_int"] = (
        transcripts["x_location"].astype(int)
    )
    transcripts["y_location_int"] = (
        transcripts["y_location"].astype(int)
    )

    height, width = mask.shape

    in_bounds = (
        transcripts["x_location_int"].between(0, width - 1)
        & transcripts["y_location_int"].between(0, height - 1)
    )

    transcripts = transcripts.loc[in_bounds].copy()

    # Assign transcripts to mask objects.
    transcripts["cell_id"] = mask[
        transcripts["y_location_int"].to_numpy(),
        transcripts["x_location_int"].to_numpy(),
    ]

    transcripts = transcripts[
        transcripts["cell_id"] != 0
    ].copy()

    if transcripts.empty:
        raise ValueError(
            f"No transcripts were assigned to objects for {csv_path}"
        )

    # Create object-by-gene count matrix.
    counts = pd.crosstab(
        transcripts["cell_id"],
        transcripts["feature_name"],
    )

    counts.index = counts.index.astype(str)

    # Calculate object centroids and areas.
    cell_ids = np.unique(mask)
    cell_ids = cell_ids[cell_ids != 0]

    centroids = ndimage.center_of_mass(
        mask,
        labels=mask,
        index=cell_ids,
    )

    areas = np.bincount(mask.ravel())[cell_ids]

    obs = pd.DataFrame(
        {
            "cell_id": cell_ids.astype(str),
            "y_centroid": [
                centroid[0] for centroid in centroids
            ],
            "x_centroid": [
                centroid[1] for centroid in centroids
            ],
            "cell_area": areas,
        }
    ).set_index("cell_id")

    # Retain only objects containing assigned transcripts.
    obs = obs.loc[counts.index]

    adata = ad.AnnData(
        X=sp.csr_matrix(counts.to_numpy()),
        obs=obs,
    )

    adata.var_names = counts.columns.astype(str)

    # Preserve the coordinate convention used in the analysis.
    adata.obsm["spatial"] = obs[
        ["y_centroid", "x_centroid"]
    ].to_numpy()

    # Calculate QC metrics.
    sc.pp.calculate_qc_metrics(
        adata,
        inplace=True,
        percent_top=None,
    )

    adata.obs["transcript_density"] = (
        adata.obs["total_counts"]
        / adata.obs["cell_area"]
    )

    # Compare count- and density-based QC filters.
    count_low, count_high = np.percentile(
        adata.obs["total_counts"],
        [5, 95],
    )

    density_low, density_high = np.percentile(
        adata.obs["transcript_density"],
        [5, 95],
    )

    keep_count = adata.obs["total_counts"].between(
        count_low,
        count_high,
    )

    keep_density = adata.obs["transcript_density"].between(
        density_low,
        density_high,
    )

    sample = os.path.splitext(
        os.path.basename(csv_path)
    )[0]

    os.makedirs(output_dir, exist_ok=True)

    comparison = adata.obs.copy()
    comparison["keep_count_filter"] = keep_count
    comparison["keep_density_filter"] = keep_density

    comparison_path = os.path.join(
        output_dir,
        f"{sample}_filtering_comparison.csv",
    )

    comparison.to_csv(comparison_path)

    print("\n--- QC FILTER COMPARISON ---")
    print(
        f"Count cutoffs: "
        f"{count_low:.2f}–{count_high:.2f}"
    )
    print(
        f"Density cutoffs: "
        f"{density_low:.6f}–{density_high:.6f}"
    )
    print(
        f"Objects retained by count filter: "
        f"{int(keep_count.sum())}"
    )
    print(
        f"Objects retained by density filter: "
        f"{int(keep_density.sum())}"
    )
    print(f"Saved comparison: {comparison_path}")

    # Apply transcript-density filtering.
    adata = adata[keep_density].copy()

    sc.pp.filter_genes(
        adata,
        min_cells=5,
    )

    print(
        f"After QC: {adata.n_obs} objects, "
        f"{adata.n_vars} genes"
    )

    # Preserve raw counts.
    adata.layers["counts"] = adata.X.copy()

    # Identify HVGs from raw counts.
    sc.pp.highly_variable_genes(
        adata,
        layer="counts",
        n_top_genes=1000,
        flavor="seurat_v3",
        subset=False,
    )

    # Normalize and log-transform.
    sc.pp.normalize_total(
        adata,
        target_sum=10_000,
    )
    sc.pp.log1p(adata)

    # Preserve normalized, unscaled expression.
    adata.raw = adata.copy()

    # Scale and calculate dimensionality reduction.
    sc.pp.scale(
        adata,
        max_value=10,
    )

    sc.tl.pca(
        adata,
        n_comps=40,
        use_highly_variable=True,
    )

    sc.pp.neighbors(
        adata,
        n_neighbors=15,
        n_pcs=30,
    )

    sc.tl.umap(adata)

    output_path = os.path.join(
        output_dir,
        f"{sample}.h5ad",
    )

    adata.write_h5ad(output_path)

    print(f"Saved AnnData: {output_path}")


def process_all(
    transcript_dir="transcripts",
    mask_dir="masks_mf",
    output_dir="anndata_mf_density",
):
    """Process all matched transcript and mask files."""

    pairs = get_file_pairs(
        transcript_dir=transcript_dir,
        mask_dir=mask_dir,
    )

    if not pairs:
        raise FileNotFoundError(
            "No matching transcript CSV and mask TIFF files were found."
        )

    print(
        f"Found {len(pairs)} matched transcript/mask pairs."
    )

    for csv_path, mask_path in pairs:
        print(
            f"\nProcessing "
            f"{os.path.basename(csv_path)}"
        )

        preprocess_sample(
            csv_path=csv_path,
            mask_path=mask_path,
            output_dir=output_dir,
        )


if __name__ == "__main__":
    process_all()