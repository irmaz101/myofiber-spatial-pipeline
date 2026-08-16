"""Compare Leiden resolutions on an existing BBKNN graph."""

from pathlib import Path

import scanpy as sc
import pandas as pd
from sklearn.metrics import adjusted_rand_score, silhouette_score


INPUT_FILE = Path("results/combined_nuclei_bbknn.h5ad")
OUTPUT_FILE = Path("results/leiden_resolution_comparison.csv")
RESOLUTIONS = (0.1, 0.2, 0.3, 0.4, 0.5)
ORIGINAL_KEY = "leiden_bbknn_r0.2"
SAMPLE_SIZE = 10_000
RANDOM_STATE = 0


def main():
    adata = sc.read_h5ad(INPUT_FILE)

    if "connectivities" not in adata.obsp:
        raise KeyError("The saved BBKNN connectivity graph is missing.")

    rows = []

    for resolution in RESOLUTIONS:
        key = f"leiden_test_r{resolution}"
        sc.tl.leiden(
            adata,
            resolution=resolution,
            adjacency=adata.obsp["connectivities"],
            key_added=key,
            random_state=RANDOM_STATE,
        )

        labels = adata.obs[key].astype(str)
        score = silhouette_score(
            adata.obsm["X_pca"],
            labels,
            sample_size=min(SAMPLE_SIZE, adata.n_obs),
            random_state=RANDOM_STATE,
        )

        ari = None
        if resolution == 0.2 and ORIGINAL_KEY in adata.obs:
            ari = adjusted_rand_score(
                adata.obs[ORIGINAL_KEY].astype(str),
                labels,
            )

        rows.append({
            "resolution": resolution,
            "n_clusters": labels.nunique(),
            "silhouette_score": score,
            "ari_vs_original_r0.2": ari,
        })

    results = pd.DataFrame(rows)
    results.to_csv(OUTPUT_FILE, index=False)
    print(results.to_string(index=False))


if __name__ == "__main__":
    main()
