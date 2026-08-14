"""
Annotate skeletal muscle fibre types using canonical Myh markers.

The script expects normalized and log-transformed expression values
in adata.X.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler


ANNDATA_DIR = Path("anndata")
RESULTS_DIR = Path("annotated_anndata")

N_CLUSTERS = 4
RANDOM_STATE = 0
N_INIT = 50

MARKER_GENES = [
    "Myh7",
    "Myh2",
    "Myh1",
    "Myh4",
]

VALID_LABELS = {
    "type_1",
    "type_2a",
    "type_2x",
    "type_2b",
    "other",
}

RESULTS_DIR.mkdir(parents=True, exist_ok=True)


def extract_marker_expression(adata):
    """Extract marker expression from adata.X."""
    gene_lookup = {
        gene.lower(): gene
        for gene in adata.var_names
    }

    missing_genes = [
        gene
        for gene in MARKER_GENES
        if gene.lower() not in gene_lookup
    ]

    if missing_genes:
        raise ValueError(
            f"Missing Myh markers: {missing_genes}"
        )

    actual_genes = [
        gene_lookup[gene.lower()]
        for gene in MARKER_GENES
    ]

    expression = adata[:, actual_genes].X

    if hasattr(expression, "toarray"):
        expression = expression.toarray()

    return np.asarray(expression), actual_genes


def request_cluster_mapping(cluster_ids):
    """Request a fibre-type label for each cluster."""
    print(
        "\nAllowed labels: "
        + ", ".join(sorted(VALID_LABELS))
    )

    cluster_mapping = {}

    for cluster_id in cluster_ids:
        while True:
            fibre_type = input(
                f"Cluster {cluster_id} fibre type: "
            ).strip().lower()

            if fibre_type in VALID_LABELS:
                cluster_mapping[cluster_id] = fibre_type
                break

            print("Invalid label.")

    return cluster_mapping


def annotate_sample(h5ad_path):
    """Cluster and annotate one myofibre dataset."""
    sample = h5ad_path.stem
    adata = sc.read_h5ad(h5ad_path)

    marker_expression, actual_genes = (
        extract_marker_expression(adata)
    )

    scaled_expression = StandardScaler().fit_transform(
        marker_expression
    )

    kmeans = KMeans(
        n_clusters=N_CLUSTERS,
        random_state=RANDOM_STATE,
        n_init=N_INIT,
    )

    cluster_ids = kmeans.fit_predict(
        scaled_expression
    ).astype(str)

    adata.obs["myh_cluster"] = pd.Categorical(
        cluster_ids
    )

    dotplot_path = (
        RESULTS_DIR
        / f"{sample}_myh_dotplot.png"
    )

    dotplot = sc.pl.dotplot(
        adata,
        var_names=actual_genes,
        groupby="myh_cluster",
        standard_scale="var",
        return_fig=True,
        show=False,
    )

    dotplot.savefig(
        dotplot_path,
        dpi=300,
        bbox_inches="tight",
    )

    print(f"\nInspect dotplot: {dotplot_path}")

    unique_clusters = sorted(
        adata.obs["myh_cluster"]
        .astype(str)
        .unique()
    )

    cluster_mapping = request_cluster_mapping(
        unique_clusters
    )

    adata.obs["fiber_type_from_myh"] = (
        adata.obs["myh_cluster"]
        .astype(str)
        .map(cluster_mapping)
        .astype("category")
    )

    adata.obsm["X_myh_scaled"] = scaled_expression

    output_path = (
        RESULTS_DIR
        / f"{sample}_with_myh_fiber_types.h5ad"
    )

    adata.write_h5ad(output_path)

    print(f"\n{sample} fibre-type counts:")
    print(
        adata.obs[
            "fiber_type_from_myh"
        ].value_counts()
    )
    print(f"Saved AnnData: {output_path}")


def main():
    h5ad_files = sorted(
        ANNDATA_DIR.glob("*.h5ad")
    )

    if not h5ad_files:
        raise FileNotFoundError(
            f"No AnnData files found in "
            f"{ANNDATA_DIR.resolve()}."
        )

    for h5ad_path in h5ad_files:
        annotate_sample(h5ad_path)


if __name__ == "__main__":
    main()
