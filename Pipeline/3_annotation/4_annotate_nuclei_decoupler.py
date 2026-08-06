"""
Run decoupler ULM annotation using curated PanglaoDB markers.
"""

from pathlib import Path

import decoupler as dc
import pandas as pd
import scanpy as sc


# --------------------------------------------------
# Settings
# --------------------------------------------------

input_file = Path("combined_nuclei_bbknn.h5ad")
results_dir = Path("results")

cluster_key = "leiden_bbknn_r0.2"
minimum_targets = 3

cell_types = [
    "Endothelial cells",
    "Fibroblasts",
    "Macrophages",
    "Pericytes",
    "Smooth muscle cells",
    "Myoblasts",
    "Schwann cells",
]

results_dir.mkdir(parents=True, exist_ok=True)


def human_to_mouse_symbol(gene):
    """Convert uppercase human gene symbols to mouse-style capitalization."""

    gene = str(gene).strip()

    if not gene:
        return gene

    return gene[0].upper() + gene[1:].lower()


def main():
    adata = sc.read_h5ad(input_file)

    if cluster_key not in adata.obs:
        raise KeyError(f"Cluster column not found: {cluster_key}")

    markers = dc.op.resource(
        "PanglaoDB",
        organism="human",
        verbose=True,
    )

    markers = markers[
        markers["mouse"].astype(bool)
        & markers["canonical_marker"].astype(bool)
        & (markers["mouse_sensitivity"].astype(float) > 0.5)
    ][["cell_type", "genesymbol"]].dropna()

    markers = markers.drop_duplicates(
        ["cell_type", "genesymbol"]
    )

    markers["genesymbol"] = (
        markers["genesymbol"]
        .map(human_to_mouse_symbol)
    )

    markers = markers.rename(
        columns={
            "cell_type": "source",
            "genesymbol": "target",
        }
    )

    markers = markers[
        markers["source"].isin(cell_types)
        & markers["target"].isin(adata.var_names)
    ][["source", "target"]].copy()

    if markers.empty:
        raise ValueError(
            "No PanglaoDB markers matched the AnnData genes."
        )

    markers.to_csv(
        results_dir / "panglaodb_markers_used.csv",
        index=False,
    )

    dc.mt.ulm(
        data=adata,
        net=markers,
        tmin=minimum_targets,
    )

    score_adata = dc.pp.get_obsm(
        adata,
        key="score_ulm",
    )

    score_df = score_adata.to_df()
    score_df[cluster_key] = adata.obs[cluster_key].astype(str).values

    cluster_scores = (
        score_df
        .groupby(cluster_key)
        .mean()
    )

    cluster_scores.to_csv(
        results_dir / "decoupler_scores_by_cluster.csv"
    )

    top_predictions = (
        cluster_scores
        .stack()
        .rename("score")
        .reset_index()
        .rename(columns={"level_1": "predicted_cell_type"})
        .sort_values(
            [cluster_key, "score"],
            ascending=[True, False],
        )
    )

    top_predictions["rank"] = (
        top_predictions
        .groupby(cluster_key)
        .cumcount()
        + 1
    )

    top_predictions = top_predictions[
        top_predictions["rank"] <= 5
    ]

    top_predictions.to_csv(
        results_dir / "decoupler_top5_predictions.csv",
        index=False,
    )

    cluster_to_label = (
        top_predictions[
            top_predictions["rank"] == 1
        ]
        .set_index(cluster_key)["predicted_cell_type"]
        .to_dict()
    )

    adata.obs["decoupler_label"] = (
        adata.obs[cluster_key]
        .astype(str)
        .map(cluster_to_label)
        .fillna("Unannotated")
        .astype("category")
    )

    adata.write_h5ad(
        results_dir / "nuclei_decoupler_annotated.h5ad"
    )

    print(
        f"Matched {len(markers)} markers across "
        f"{markers['source'].nunique()} cell types."
    )


if __name__ == "__main__":
    main()