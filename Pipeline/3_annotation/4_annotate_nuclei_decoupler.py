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

input_file = Path("results/annotated_nuclei.h5ad")
results_dir = Path("results/decoupler")

cluster_key = "leiden_bbknn_r0.2"
minimum_targets = 3

cell_types = [
    "Endothelial cells",
    "Fibroblasts",
    "Macrophages",
    "Pericytes",
    "Smooth muscle cells",
    "Myocytes",
    "Myoblasts",
    "Satellite cells",
    "Schwann cells",
]

results_dir.mkdir(parents=True, exist_ok=True)


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

    if adata.raw is None:
        raise ValueError("adata.raw is required for ULM scoring.")

    scoring = adata.raw.to_adata()
    scoring.obs = adata.obs.copy()

    # Retain mouse-supported PanglaoDB entries and match identical
    # symbols case-insensitively to genes measured by Xenium.
    gene_lookup = {
        str(gene).upper(): str(gene)
        for gene in scoring.var_names
    }

    markers["panglaodb_symbol"] = (
        markers["genesymbol"].astype(str).str.strip()
    )
    markers["target"] = (
        markers["panglaodb_symbol"].str.upper().map(gene_lookup)
    )

    markers = markers.rename(
        columns={
            "cell_type": "source",
        }
    )

    markers = markers[
        markers["source"].isin(cell_types)
        & markers["target"].notna()
    ][["source", "target", "panglaodb_symbol"]].drop_duplicates(
        ["source", "target"]
    ).copy()

    coverage = (
        markers.groupby("source")
        .size()
        .rename("n_matched_markers")
        .reset_index()
    )
    coverage.to_csv(
        results_dir / "panglaodb_marker_counts_in_xenium.csv",
        index=False,
    )

    valid_sources = coverage.loc[
        coverage["n_matched_markers"] >= minimum_targets,
        "source",
    ]
    markers = markers[
        markers["source"].isin(valid_sources)
    ].copy()

    if markers.empty:
        raise ValueError(
            "No PanglaoDB markers matched the AnnData genes."
        )

    markers.to_csv(
        results_dir / "panglaodb_mouse_supported_markers_used.csv",
        index=False,
    )

    dc.mt.ulm(
        data=scoring,
        net=markers[["source", "target"]],
        tmin=minimum_targets,
    )

    score_adata = dc.pp.get_obsm(
        scoring,
        key="score_ulm",
    )

    score_df = score_adata.to_df()
    score_df[cluster_key] = adata.obs[cluster_key].astype(str).values

    adata.obsm["panglaodb_ulm_scores"] = score_df.drop(
        columns=[cluster_key]
    ).copy()

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

    adata.obs["panglaodb_top_signature"] = (
        adata.obs[cluster_key]
        .astype(str)
        .map(cluster_to_label)
        .fillna("Unannotated")
        .astype("category")
    )

    adata.write_h5ad(
        results_dir / "nuclei_panglaodb_support.h5ad"
    )

    print(
        f"Matched {len(markers)} markers across "
        f"{markers['source'].nunique()} cell types."
    )


if __name__ == "__main__":
    main()
