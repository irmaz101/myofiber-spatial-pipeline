"""Integrate nuclear count matrices with BBKNN, then run UMAP and Leiden.

Inputs must contain unnormalized counts in .X."""

from pathlib import Path

import anndata as ad
import bbknn
import pandas as pd
import scanpy as sc

# Settings

anndata_dir = Path("nuclei_anndata")
results_dir = Path("results")

target_sum = 10_000
n_highly_variable_genes = 2_000
n_pcs = 30
neighbors_within_batch = 3
leiden_resolution = 0.2
random_state = 0

cluster_key = f"leiden_bbknn_r{leiden_resolution}"

results_dir.mkdir(parents=True, exist_ok=True)
sc.settings.figdir = results_dir
sc.settings.autoshow = False


# Load and combine samples


def load_and_combine():
    """Load all sample-level AnnData files and combine them."""

    h5ad_files = sorted(anndata_dir.glob("*.h5ad"))

    if not h5ad_files:
        raise FileNotFoundError(
            f"No AnnData files were found in {anndata_dir.resolve()}."
        )

    adata_list = []

    for h5ad_path in h5ad_files:
        sample = h5ad_path.stem.replace(".filtered", "")

        print(f"Loading {sample}")

        adata = sc.read_h5ad(h5ad_path)
        adata.obs["sample"] = sample
        adata.obs_names = [f"{sample}_{object_id}" for object_id in adata.obs_names]

        adata_list.append(adata)

    combined = ad.concat(adata_list, join="outer", merge="same", index_unique=None)

    print(f"Combined dataset: {combined.n_obs} nuclei × " f"{combined.n_vars} genes")

    return combined


# BBKNN integration and clustering


def integrate_nuclei(adata):
    """Normalize, integrate and cluster the combined nucleus dataset."""

    # Preserve raw counts.
    adata.layers["counts"] = adata.X.copy()

    # Select HVGs from counts, accounting for sample.
    sc.pp.highly_variable_genes(
        adata,
        layer="counts",
        n_top_genes=n_highly_variable_genes,
        flavor="seurat_v3",
        batch_key="sample",
        subset=False,
    )

    print(f"HVGs detected: " f"{int(adata.var['highly_variable'].sum())}")

    sc.pp.normalize_total(adata, target_sum=target_sum)
    sc.pp.log1p(adata)

    # Retain all genes in .raw for marker and signature scoring.
    adata.raw = adata.copy()

    # Use HVGs for dimensionality reduction and integration.
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()

    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, n_comps=n_pcs, random_state=random_state)

    bbknn.bbknn(
        adata_hvg,
        batch_key="sample",
        use_rep="X_pca",
        neighbors_within_batch=neighbors_within_batch,
        n_pcs=n_pcs,
    )

    sc.tl.umap(adata_hvg, random_state=random_state)

    sc.tl.leiden(
        adata_hvg,
        resolution=leiden_resolution,
        key_added=cluster_key,
        random_state=random_state,
    )

    print(f"Leiden clusters: " f"{adata_hvg.obs[cluster_key].nunique()}")

    # Keep the integrated graph and embeddings with the full gene matrix.
    adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"].copy()
    adata.obsm["X_umap_bbknn"] = adata_hvg.obsm["X_umap"].copy()
    adata.obsp["connectivities"] = adata_hvg.obsp["connectivities"].copy()
    adata.obsp["distances"] = adata_hvg.obsp["distances"].copy()
    adata.uns["neighbors"] = adata_hvg.uns["neighbors"].copy()
    adata.obs[cluster_key] = adata_hvg.obs[cluster_key].copy()
    adata.obs["leiden"] = adata.obs[cluster_key].copy()

    return adata


def save_results(adata):
    """Save the integrated object, cluster counts and UMAP plots."""

    cluster_counts = pd.crosstab(adata.obs[cluster_key], adata.obs["sample"])

    counts_path = results_dir / "cluster_counts_by_sample.csv"
    cluster_counts.to_csv(counts_path)

    sc.pl.embedding(
        adata,
        basis="X_umap_bbknn",
        color=cluster_key,
        frameon=False,
        size=5,
        show=False,
        save="_bbknn_clusters.png",
    )

    sc.pl.embedding(
        adata,
        basis="X_umap_bbknn",
        color="sample",
        frameon=False,
        size=5,
        show=False,
        save="_bbknn_samples.png",
    )

    output_path = results_dir / "combined_nuclei_bbknn.h5ad"
    adata.write_h5ad(output_path)

    print(f"Saved integrated AnnData: {output_path}")
    print(f"Saved cluster counts: {counts_path}")


def main():
    adata = load_and_combine()
    adata = integrate_nuclei(adata)
    save_results(adata)


if __name__ == "__main__":
    main()
