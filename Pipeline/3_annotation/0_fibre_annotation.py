"""
Annotate skeletal muscle fibre types using canonical Myh markers.

For each sample, fibres are clustered using z-scaled expression of
Myh1, Myh2, Myh4, and Myh7. Cluster identities are assigned manually
after inspection of a marker-expression dotplot.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler


# --------------------------------------------------
# Settings
# --------------------------------------------------

anndata_dir = Path("anndata")
results_dir = Path("results")

n_clusters = 4
random_state = 0
n_init = 50

markers = {
    "type_1": ["Myh7"],
    "type_2a": ["Myh2"],
    "type_2x": ["Myh1"],
    "type_2b": ["Myh4"],
}

marker_genes = [
    gene
    for genes in markers.values()
    for gene in genes
]

valid_labels = {
    "type_1",
    "type_2a",
    "type_2x",
    "type_2b",
    "other",
}

results_dir.mkdir(parents=True, exist_ok=True)


# --------------------------------------------------
# Helper functions
# --------------------------------------------------

def extract_marker_expression(adata, genes):
    """Extract marker expression using case-insensitive gene matching."""

    if adata.raw is not None:
        gene_names = pd.Index(adata.raw.var_names)
        source = adata.raw
    else:
        gene_names = pd.Index(adata.var_names)
        source = adata

    lowercase_names = {
        gene.lower(): gene
        for gene in gene_names
    }

    expression = np.zeros(
        (adata.n_obs, len(genes)),
        dtype=float,
    )

    found_genes = {}

    for column, requested_gene in enumerate(genes):
        actual_gene = lowercase_names.get(
            requested_gene.lower()
        )

        if actual_gene is None:
            continue

        values = source[:, actual_gene].X

        if hasattr(values, "toarray"):
            values = values.toarray()

        expression[:, column] = np.asarray(values).ravel()
        found_genes[requested_gene] = actual_gene

    if not found_genes:
        raise ValueError(
            "None of the requested Myh markers were found."
        )

    missing_genes = [
        gene
        for gene in genes
        if gene not in found_genes
    ]

    if missing_genes:
        print(f"Missing markers: {missing_genes}")

    return expression


def request_cluster_labels(cluster_ids):
    """Ask the user to assign a fibre type to each k-means cluster."""

    cluster_mapping = {}

    print(
        "\nAllowed labels: "
        "type_1, type_2a, type_2x, type_2b, other"
    )

    for cluster_id in cluster_ids:
        while True:
            fibre_type = input(
                f"Cluster {cluster_id} fibre type: "
            ).strip().lower()

            if fibre_type in valid_labels:
                cluster_mapping[str(cluster_id)] = fibre_type
                break

            print("Invalid label. Please try again.")

    return cluster_mapping


def annotate_sample(h5ad_path):
    """Cluster and manually annotate one myofibre dataset."""

    sample = h5ad_path.stem

    print(f"\nProcessing {sample}")

    adata = sc.read_h5ad(h5ad_path)

    marker_expression = extract_marker_expression(
        adata,
        marker_genes,
    )

    scaled_expression = StandardScaler().fit_transform(
        marker_expression
    )

    kmeans = KMeans(
        n_clusters=n_clusters,
        random_state=random_state,
        n_init=n_init,
    )

    adata.obs["myh_cluster"] = pd.Categorical(
        kmeans.fit_predict(scaled_expression).astype(str)
    )

    # Display marker profiles before manual annotation.
    sc.pl.dotplot(
        adata,
        var_names=markers,
        groupby="myh_cluster",
        standard_scale="var",
        show=False,
    )

    dotplot_path = results_dir / f"{sample}_myh_dotplot.png"

    plt.savefig(
        dotplot_path,
        dpi=300,
        bbox_inches="tight",
    )
    plt.show()
    plt.close()

    cluster_ids = sorted(
        adata.obs["myh_cluster"]
        .astype(str)
        .unique()
    )

    cluster_mapping = request_cluster_labels(
        cluster_ids
    )

    adata.obs["fiber_type"] = (
        adata.obs["myh_cluster"]
        .astype(str)
        .map(cluster_mapping)
        .fillna("other")
        .astype("category")
    )

    adata.obsm["X_myh_scaled"] = scaled_expression

    mapping_path = (
        results_dir /
        f"{sample}_cluster_mapping.csv"
    )

    pd.DataFrame(
        cluster_mapping.items(),
        columns=["myh_cluster", "fiber_type"],
    ).to_csv(mapping_path, index=False)

    output_path = (
        results_dir /
        f"{sample}_fiber_annotated.h5ad"
    )

    adata.write_h5ad(output_path)

    print("\nFibre-type counts:")
    print(adata.obs["fiber_type"].value_counts())
    print(f"Saved dotplot: {dotplot_path}")
    print(f"Saved mapping: {mapping_path}")
    print(f"Saved AnnData: {output_path}")


# --------------------------------------------------
# Main analysis
# --------------------------------------------------

def main():
    """Annotate all AnnData files in the input directory."""

    h5ad_files = sorted(anndata_dir.glob("*.h5ad"))

    if not h5ad_files:
        raise FileNotFoundError(
            f"No AnnData files were found in "
            f"{anndata_dir.resolve()}."
        )

    for h5ad_path in h5ad_files:
        output_path = (
            results_dir /
            f"{h5ad_path.stem}_fiber_annotated.h5ad"
        )

        if output_path.exists():
            print(
                f"Skipping {h5ad_path.stem}: "
                "annotated output already exists."
            )
            continue

        try:
            annotate_sample(h5ad_path)

        except Exception as error:
            print(
                f"Annotation failed for "
                f"{h5ad_path.stem}: {error}"
            )


if __name__ == "__main__":
    main()