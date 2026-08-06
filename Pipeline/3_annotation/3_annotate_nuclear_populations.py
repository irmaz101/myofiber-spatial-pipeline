"""
Annotate nuclear populations in the integrated BBKNN dataset.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc


# --------------------------------------------------
# Settings
# --------------------------------------------------

input_file = Path("combined_nuclei_bbknn.h5ad")
results_dir = Path("results")

cluster_key = "leiden_bbknn_r0.2"

cluster_labels = {
    "0": "Myonuclei",
    "1": "Endothelial",
    "2": "FAPs",
    "3": "Vascular mural",
    "4": "Macrophages",
    "5": "Schwann",
    "6": "Endothelial",
    "7": "Satellite / myogenic",
    "8": "Lymphatic endothelial",
}

marker_genes = {
    "Myonuclei": ["Ttn", "Mybpc1", "Mybpc2", "Myh1", "Myh2", "Myh4"],
    "Endothelial": ["Pecam1", "Cdh5", "Kdr", "Esam"],
    "FAPs": ["Dcn", "Col1a1", "Col3a1", "Pdgfra"],
    "Vascular mural": ["Myh11", "Notch3", "Pdgfrb"],
    "Macrophages": ["Mrc1", "Csf1r", "Adgre1", "Cd163"],
    "Schwann": ["Cnp", "Itgb4", "Cd9"],
    "Satellite / myogenic": ["Pax7", "Myf5", "Fgfr4"],
    "Lymphatic endothelial": ["Lyve1", "Ccl21a", "Flt4"],
}

results_dir.mkdir(parents=True, exist_ok=True)


def main():
    adata = sc.read_h5ad(input_file)

    if cluster_key not in adata.obs:
        raise KeyError(f"Cluster column not found: {cluster_key}")

    # Differential expression by Leiden cluster
    sc.tl.rank_genes_groups(
        adata,
        groupby=cluster_key,
        method="wilcoxon",
    )

    de_results = sc.get.rank_genes_groups_df(
        adata,
        group=None,
    )

    de_results.to_csv(
        results_dir / "nuclear_cluster_de_genes.csv",
        index=False,
    )

    # Assign manually curated population labels
    adata.obs["nuclear_population"] = (
        adata.obs[cluster_key]
        .astype(str)
        .map(cluster_labels)
        .astype("category")
    )

    # Keep only markers present in the dataset
    markers_present = {
        population: [
            gene for gene in genes
            if gene in adata.var_names
        ]
        for population, genes in marker_genes.items()
    }

    markers_present = {
        population: genes
        for population, genes in markers_present.items()
        if genes
    }

    sc.pl.dotplot(
        adata,
        var_names=markers_present,
        groupby="nuclear_population",
        standard_scale="var",
        show=False,
    )

    plt.savefig(
        results_dir / "nuclear_marker_dotplot.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()

    population_summary = (
        adata.obs["nuclear_population"]
        .value_counts()
        .rename_axis("population")
        .reset_index(name="n_nuclei")
    )

    population_summary["percentage"] = (
        population_summary["n_nuclei"]
        / adata.n_obs
        * 100
    )

    population_summary.to_csv(
        results_dir / "nuclear_population_summary.csv",
        index=False,
    )

    output_file = results_dir / "annotated_nuclei.h5ad"
    adata.write_h5ad(output_file)

    print(f"Saved annotated object: {output_file}")


if __name__ == "__main__":
    main()