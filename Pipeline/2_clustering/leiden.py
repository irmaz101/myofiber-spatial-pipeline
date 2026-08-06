"""
Run Leiden clustering on myofibre AnnData objects.

For each sample, the script performs Leiden clustering, calculates
silhouette scores, calculates the Adjusted Rand Index when ground-truth
labels are available, generates spatial and UMAP plots, and saves the
clustered AnnData object.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
from sklearn.metrics import adjusted_rand_score, silhouette_score


# --------------------------------------------------
# Settings
# --------------------------------------------------

anndata_dir = Path("anndata_myofibres")
ground_truth_file = Path("ground_truth_fibers.csv")

resolution = 0.2
label_column = f"leiden_r{resolution}"

output_csv = Path(f"leiden_r{resolution}_results.csv")
plot_dir = Path(f"leiden_r{resolution}_plots")
anndata_output_dir = Path(f"leiden_r{resolution}_anndata")

plot_dir.mkdir(parents=True, exist_ok=True)
anndata_output_dir.mkdir(parents=True, exist_ok=True)


# --------------------------------------------------
# Load ground-truth annotations
# --------------------------------------------------

if ground_truth_file.exists():
    ground_truth = pd.read_csv(ground_truth_file)

    required_gt_columns = {"sample", "cell_id", "fiber_type"}
    missing_columns = required_gt_columns.difference(ground_truth.columns)

    if missing_columns:
        raise ValueError(
            "The ground-truth file is missing the following columns: "
            f"{sorted(missing_columns)}"
        )

    ground_truth["cell_id"] = ground_truth["cell_id"].astype(str)
else:
    print(
        f"Ground-truth file not found ({ground_truth_file}); "
        "ARI will not be calculated."
    )
    ground_truth = pd.DataFrame(
        columns=["sample", "cell_id", "fiber_type"]
    )


# --------------------------------------------------
# Process samples
# --------------------------------------------------

anndata_files = sorted(anndata_dir.glob("*.h5ad"))

if not anndata_files:
    raise FileNotFoundError(
        f"No AnnData files were found in {anndata_dir.resolve()}."
    )

all_results = []

for path in anndata_files:

    sample = path.stem

    print(f"\nProcessing {sample}")

    adata = sc.read_h5ad(path)
    adata.obs.index = adata.obs.index.astype(str)
    adata.obs.index.name = "cell_id"

    if "X_pca" not in adata.obsm:
        print(f"Skipping {sample}: X_pca not found.")
        continue

    if "neighbors" not in adata.uns:
        print(f"Computing neighbourhood graph for {sample}.")
        sc.pp.neighbors(
            adata,
            use_rep="X_pca"
        )

    if "X_umap" not in adata.obsm:
        print(f"Computing UMAP for {sample}.")
        sc.tl.umap(adata)

    sc.tl.leiden(
        adata,
        resolution=resolution,
        key_added=label_column
    )

    cluster_labels = adata.obs[label_column].astype(str)
    number_of_clusters = cluster_labels.nunique()

    silhouette_pca = None
    silhouette_umap = None

    if number_of_clusters > 1:
        silhouette_pca = silhouette_score(
            adata.obsm["X_pca"],
            cluster_labels
        )

        silhouette_umap = silhouette_score(
            adata.obsm["X_umap"],
            cluster_labels
        )

    # --------------------------------------------------
    # Calculate ARI for samples with ground truth
    # --------------------------------------------------

    ground_truth_sample = ground_truth[
        ground_truth["sample"] == sample
    ].copy()

    ari = None
    number_of_gt_cells = len(ground_truth_sample)
    number_of_matched_cells = 0

    if not ground_truth_sample.empty:

        predicted_labels = (
            adata.obs[[label_column]]
            .copy()
            .reset_index()
            .rename(columns={label_column: "leiden_cluster"})
        )

        predicted_labels["cell_id"] = (
            predicted_labels["cell_id"].astype(str)
        )

        predicted_labels["sample"] = sample

        matched_cells = ground_truth_sample.merge(
            predicted_labels,
            on=["sample", "cell_id"],
            how="inner"
        )

        number_of_matched_cells = len(matched_cells)

        if not matched_cells.empty:
            ari = adjusted_rand_score(
                matched_cells["fiber_type"],
                matched_cells["leiden_cluster"]
            )
        else:
            print(f"No matching ground-truth cells for {sample}.")

    else:
        print(
            f"No ground-truth labels for {sample}; "
            "silhouette scores calculated only."
        )

    all_results.append({
        "sample": sample,
        "resolution": resolution,
        "n_gt_cells": number_of_gt_cells,
        "n_matched_cells": number_of_matched_cells,
        "n_clusters": number_of_clusters,
        "ari_leiden": ari,
        "silhouette_pca": silhouette_pca,
        "silhouette_umap": silhouette_umap
    })

    print(
        f"Clusters={number_of_clusters}, "
        f"ARI={ari}, "
        f"SC_PCA={silhouette_pca}, "
        f"SC_UMAP={silhouette_umap}"
    )

    # --------------------------------------------------
    # Save spatial plot
    # --------------------------------------------------

    sc.pl.spatial(
        adata,
        color=label_column,
        spot_size=20,
        show=False,
        title=f"{sample}: Leiden resolution {resolution}"
    )

    spatial_plot_path = (
        plot_dir /
        f"{sample}_leiden_r{resolution}_spatial.png"
    )

    plt.savefig(
        spatial_plot_path,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

    # --------------------------------------------------
    # Save UMAP plot
    # --------------------------------------------------

    sc.pl.umap(
        adata,
        color=label_column,
        show=False,
        title=f"{sample}: Leiden resolution {resolution}"
    )

    umap_plot_path = (
        plot_dir /
        f"{sample}_leiden_r{resolution}_umap.png"
    )

    plt.savefig(
        umap_plot_path,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

    # --------------------------------------------------
    # Save clustered AnnData
    # --------------------------------------------------

    anndata_output_path = (
        anndata_output_dir /
        f"{sample}_leiden_r{resolution}.h5ad"
    )

    adata.write_h5ad(anndata_output_path)


# --------------------------------------------------
# Save evaluation results
# --------------------------------------------------

results = pd.DataFrame(all_results)
results.to_csv(output_csv, index=False)

print(f"\nSaved results to: {output_csv}")
print(f"Saved plots to: {plot_dir}")
print(f"Saved AnnData files to: {anndata_output_dir}")
print(results)
