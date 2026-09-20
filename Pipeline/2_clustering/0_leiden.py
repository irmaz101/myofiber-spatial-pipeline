"""Benchmark Leiden resolutions against PCA separation and manual fibre labels.

The reference classes are IIa, IIb and double-negative fibres. Resolution
0.2 is retained in the output; annotation and SVG analyses use the original
preprocessed files separately."""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse
from sklearn.metrics import adjusted_rand_score, silhouette_score

# Settings

INPUT_DIR = Path("anndata_mf_density")
GROUND_TRUTH_FILE = Path("ground_truth_fibers.csv")
OUTPUT_DIR = Path("leiden_resolution_results_3class")

RESOLUTIONS = [0.1, 0.2, 0.3, 0.4, 0.5]
SELECTED_RESOLUTION = 0.2

N_TOP_GENES = 1000
TARGET_SUM = 10_000
N_PCS = 40
N_NEIGHBORS = 15
N_PCS_NEIGHBORS = 30
RANDOM_STATE = 0

LABEL_ALIASES = {
    "type_2a": "2a iia type_iia type_2a".split(),
    "type_2b": "2b iib type_iib type_2b".split(),
    "double_negative": (
        "2x iix type_iix type_2x negative double-negative "
        "double_negative double_negative_fibre"
    ).split(),
}

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def normalize_ids(values):
    """Convert object identifiers to comparable strings."""
    values = pd.Series(values, dtype="object")
    numeric = pd.to_numeric(values, errors="coerce")
    normalized = values.astype(str).str.strip()
    numeric_mask = numeric.notna()
    normalized.loc[numeric_mask] = numeric.loc[numeric_mask].astype("Int64").astype(str)
    return normalized.to_numpy()


def normalize_sample_names(values):
    """Standardize sample names."""
    return (
        pd.Series(values, dtype="object")
        .astype(str)
        .str.strip()
        .str.removesuffix(".h5ad")
        .to_numpy()
    )


def load_ground_truth():
    """Load and validate the manual fibre labels."""
    if not GROUND_TRUTH_FILE.exists():
        print(
            f"Ground-truth file not found: {GROUND_TRUTH_FILE}. ARI will not be calculated."
        )
        return None
    truth = pd.read_csv(GROUND_TRUTH_FILE)
    columns = ["sample", "cell_id", "fiber_type"]
    missing_columns = set(columns) - set(truth.columns)
    if missing_columns:
        raise ValueError(
            f"Ground-truth file is missing columns: {sorted(missing_columns)}"
        )
    truth = truth[columns].copy()
    truth["sample"] = normalize_sample_names(truth["sample"])
    truth["cell_id"] = normalize_ids(truth["cell_id"])
    truth["fiber_type"] = (
        truth["fiber_type"]
        .astype(str)
        .str.strip()
        .str.lower()
        .str.replace(" ", "_")
        .replace(
            {
                alias: label
                for label, aliases in LABEL_ALIASES.items()
                for alias in aliases
            }
        )
    )

    allowed_labels = set(LABEL_ALIASES)
    unexpected_labels = set(truth["fiber_type"]) - allowed_labels
    if unexpected_labels:
        raise ValueError(f"Unexpected ground-truth labels: {sorted(unexpected_labels)}")

    # A reference fibre must not have conflicting labels.
    conflicting = truth.groupby(["sample", "cell_id"])["fiber_type"].nunique().gt(1)
    if conflicting.any():
        conflicts = conflicting[conflicting].reset_index()[["sample", "cell_id"]]
        raise ValueError(
            "Conflicting ground-truth labels for these sample–cell pairs:\n"
            f"{conflicts.to_string(index=False)}"
        )

    n_before = len(truth)
    truth = truth.drop_duplicates(
        subset=["sample", "cell_id"], keep="first"
    ).reset_index(drop=True)
    if len(truth) < n_before:
        print(f"Removed {n_before - len(truth)} identical duplicate ground-truth rows.")

    reference_classes = sorted(truth["fiber_type"].unique())
    if set(reference_classes) != allowed_labels:
        raise ValueError(
            f"Expected {sorted(allowed_labels)} in the ground truth, but found: {reference_classes}"
        )
    print(f"Ground-truth classes: {reference_classes}")
    print("Ground-truth counts:")
    print(truth["fiber_type"].value_counts().sort_index().to_string())
    return truth


def get_raw_count_adata(original, sample):
    """Construct a new AnnData object from raw integer counts."""
    if original.raw is None:
        raise ValueError(f"{sample}: adata.raw is missing.")
    counts = original.raw.X.copy()
    if issparse(counts):
        stored_values = counts.data
    else:
        counts = np.asarray(counts)
        stored_values = counts.ravel()
    if stored_values.size:
        if stored_values.min() < 0:
            raise ValueError(f"{sample}: adata.raw.X contains negative values.")
        if not np.allclose(stored_values, np.round(stored_values)):
            raise ValueError(
                f"{sample}: adata.raw.X is not integer-valued and may not contain raw counts."
            )

    adata = ad.AnnData(
        X=counts,
        obs=original.obs.copy(),
        var=pd.DataFrame(index=original.raw.var_names.copy()),
    )
    adata.obs_names = normalize_ids(original.obs_names)
    if adata.obs_names.duplicated().any():
        duplicated = adata.obs_names[adata.obs_names.duplicated(keep=False)].unique()
        raise ValueError(f"{sample}: duplicated object IDs: {duplicated.tolist()}")

    if "spatial" in original.obsm:
        adata.obsm["spatial"] = np.asarray(original.obsm["spatial"]).copy()
    elif {"x_centroid", "y_centroid"}.issubset(original.obs.columns):
        adata.obsm["spatial"] = original.obs[["x_centroid", "y_centroid"]].to_numpy(
            dtype=float
        )
    adata.var_names_make_unique()
    adata.layers["counts"] = adata.X.copy()
    return adata


def preprocess_for_clustering(adata):
    """Select HVGs and calculate the expression graph for Leiden."""
    sc.pp.highly_variable_genes(
        adata, layer="counts", n_top_genes=N_TOP_GENES, flavor="seurat_v3", subset=False
    )

    sc.pp.normalize_total(adata, target_sum=TARGET_SUM)
    sc.pp.log1p(adata)

    # Keep unscaled log expression in this clustering output.
    adata.raw = adata.copy()

    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=N_PCS, use_highly_variable=True, random_state=RANDOM_STATE)
    sc.pp.neighbors(
        adata, n_neighbors=N_NEIGHBORS, n_pcs=N_PCS_NEIGHBORS, random_state=RANDOM_STATE
    )
    sc.tl.umap(adata, random_state=RANDOM_STATE)


def calculate_silhouette(adata, cluster_key):
    """Calculate the silhouette score in PCA space."""
    labels = adata.obs[cluster_key].astype(str).to_numpy()
    n_clusters = len(np.unique(labels))
    if not 1 < n_clusters < adata.n_obs:
        return np.nan
    return silhouette_score(adata.obsm["X_pca"][:, :N_PCS_NEIGHBORS], labels)


def calculate_ground_truth_ari(adata, sample_ground_truth, cluster_key):
    """Calculate ARI on matched manually annotated fibres."""
    cluster_table = pd.DataFrame(
        {
            "cell_id": adata.obs_names.astype(str),
            "cluster": (adata.obs[cluster_key].astype(str).to_numpy()),
        }
    )
    matched = sample_ground_truth.merge(
        cluster_table, on="cell_id", how="inner", validate="one_to_one"
    )
    if (
        len(matched) < 2
        or matched["fiber_type"].nunique() < 2
        or matched["cluster"].nunique() < 2
    ):
        ari = np.nan
    else:
        ari = adjusted_rand_score(matched["fiber_type"], matched["cluster"])
    return ari, matched


def main():
    ground_truth = load_ground_truth()
    input_files = sorted(INPUT_DIR.glob("*.h5ad"))
    if not input_files:
        raise FileNotFoundError(f"No .h5ad files found in {INPUT_DIR.resolve()}.")

    result_rows = []
    cluster_count_rows = []
    matching_rows = []

    for position, input_file in enumerate(input_files, start=1):
        sample = input_file.stem
        print(f"\n[{position}/{len(input_files)}] Processing {sample}")
        original = sc.read_h5ad(input_file)
        adata = get_raw_count_adata(original, sample)
        preprocess_for_clustering(adata)
        if ground_truth is None:
            sample_ground_truth = pd.DataFrame(
                columns=["sample", "cell_id", "fiber_type"]
            )
        else:
            sample_ground_truth = ground_truth[ground_truth["sample"] == sample].copy()
        matched_ids = set(adata.obs_names.astype(str)).intersection(
            sample_ground_truth["cell_id"]
        )
        matching_rows.append(
            {
                "sample": sample,
                "n_objects": adata.n_obs,
                "n_ground_truth": len(sample_ground_truth),
                "n_ground_truth_matched": len(matched_ids),
                "percent_ground_truth_matched": (
                    100 * len(matched_ids) / len(sample_ground_truth)
                    if len(sample_ground_truth)
                    else np.nan
                ),
            }
        )
        for resolution in RESOLUTIONS:
            cluster_key = f"leiden_res_{resolution}"
            sc.tl.leiden(
                adata,
                resolution=resolution,
                key_added=cluster_key,
                random_state=RANDOM_STATE,
                flavor="igraph",
                n_iterations=2,
                directed=False,
            )
            n_clusters = adata.obs[cluster_key].nunique()
            silhouette = calculate_silhouette(adata, cluster_key)
            ari, matched = calculate_ground_truth_ari(
                adata, sample_ground_truth, cluster_key
            )
            result_rows.append(
                {
                    "sample": sample,
                    "resolution": resolution,
                    "is_selected_resolution": (resolution == SELECTED_RESOLUTION),
                    "n_objects": adata.n_obs,
                    "n_clusters": n_clusters,
                    "silhouette_score": silhouette,
                    "n_ground_truth": len(sample_ground_truth),
                    "n_ground_truth_matched": len(matched),
                    "ground_truth_ari": ari,
                }
            )
            cluster_sizes = (
                adata.obs[cluster_key].astype(str).value_counts().sort_index()
            )
            cluster_count_rows.extend(
                {
                    "sample": sample,
                    "resolution": resolution,
                    "cluster": cluster,
                    "n_objects": int(size),
                }
                for cluster, size in cluster_sizes.items()
            )
            ari_text = f"{ari:.4f}" if np.isfinite(ari) else "NA"
            print(
                f"  Resolution {resolution}: {n_clusters} clusters, "
                f"silhouette={silhouette:.4f}, ARI={ari_text}, matched={len(matched)}"
            )
        selected_key = f"leiden_res_{SELECTED_RESOLUTION}"
        if selected_key not in adata.obs:
            raise KeyError(f"Selected Leiden result {selected_key} was not generated.")
        adata.obs["leiden"] = adata.obs[selected_key].copy()
        adata.uns["selected_leiden_resolution"] = SELECTED_RESOLUTION
        output_path = OUTPUT_DIR / f"{sample}_leiden_resolutions.h5ad"
        adata.write_h5ad(output_path)
        print(f"  Saved: {output_path}")

    results = pd.DataFrame(result_rows)
    summary = results.groupby("resolution", as_index=False).agg(
        n_samples=("sample", "nunique"),
        n_clusters_mean=("n_clusters", "mean"),
        n_clusters_sd=("n_clusters", "std"),
        silhouette_mean=("silhouette_score", "mean"),
        silhouette_sd=("silhouette_score", "std"),
        ari_mean=("ground_truth_ari", "mean"),
        ari_sd=("ground_truth_ari", "std"),
        n_samples_with_ground_truth=("ground_truth_ari", "count"),
    )

    tables = {
        "leiden_resolution_results_per_sample.csv": results,
        "leiden_cluster_counts_by_resolution.csv": pd.DataFrame(cluster_count_rows),
        "ground_truth_matching_summary.csv": pd.DataFrame(matching_rows),
        "leiden_resolution_summary.csv": summary,
    }
    for filename, table in tables.items():
        table.to_csv(OUTPUT_DIR / filename, index=False)
    print("\nResolution summary:")
    print(summary.round(4).to_string(index=False))
    print(f"\nResults saved in: {OUTPUT_DIR.resolve()}")


if __name__ == "__main__":
    main()
