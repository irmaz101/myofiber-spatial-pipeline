"""Apply per-sample total-count QC to nucleus-level AnnData objects."""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc


INPUT_DIR = Path("nuclei_anndata_unfiltered")
OUTPUT_DIR = Path("nuclei_anndata")
LOW_PERCENTILE = 5
HIGH_PERCENTILE = 95
MIN_CELLS_PER_GENE = 5


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    input_files = sorted(INPUT_DIR.glob("*.h5ad"))

    if not input_files:
        raise FileNotFoundError(
            f"No .h5ad files found in {INPUT_DIR.resolve()}."
        )

    summaries = []

    for path in input_files:
        adata = ad.read_h5ad(path)
        sample = path.stem
        adata.obs["sample"] = sample

        sc.pp.calculate_qc_metrics(
            adata,
            inplace=True,
            percent_top=None,
        )

        counts = adata.obs["total_counts"].to_numpy()
        low = float(np.percentile(counts, LOW_PERCENTILE))
        high = float(np.percentile(counts, HIGH_PERCENTILE))
        keep = (counts >= low) & (counts <= high)

        n_before = adata.n_obs
        adata = adata[keep].copy()
        sc.pp.filter_genes(adata, min_cells=MIN_CELLS_PER_GENE)

        output_path = OUTPUT_DIR / f"{sample}.filtered.h5ad"
        adata.write_h5ad(output_path)

        summaries.append({
            "sample": sample,
            "n_nuclei_before_qc": n_before,
            "lower_total_count_cutoff": low,
            "upper_total_count_cutoff": high,
            "n_nuclei_after_qc": adata.n_obs,
            "n_nuclei_removed": n_before - adata.n_obs,
            "n_genes_after_filtering": adata.n_vars,
        })

    pd.DataFrame(summaries).to_csv(
        OUTPUT_DIR / "nuclei_qc_summary.csv",
        index=False,
    )


if __name__ == "__main__":
    main()
