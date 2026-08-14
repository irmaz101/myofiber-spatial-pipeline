"""
Assign transcripts to segmented objects and preprocess each sample.

Transcript coordinates must already be expressed in mask pixel
coordinates. Only genes listed in genes.txt are retained.
"""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import tifffile
from scipy import ndimage


TRANSCRIPT_DIR = Path("transcripts")
MASK_DIR = Path("masks_mf")
GENE_FILE = Path("genes.txt")
OUTPUT_DIR = Path("anndata_mf_density")

DENSITY_PERCENTILES = (5, 95)
MIN_CELLS_PER_GENE = 5
N_HVGS = 1000
RANDOM_STATE = 0

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def load_genes():
    """Load genes retained for downstream analysis."""
    if not GENE_FILE.exists():
        raise FileNotFoundError(
            f"Gene list not found: {GENE_FILE}"
        )

    genes = set(
        pd.read_csv(
            GENE_FILE,
            header=None,
        )[0]
        .dropna()
        .astype(str)
    )

    if not genes:
        raise ValueError(
            f"No genes found in {GENE_FILE}."
        )

    return genes


def get_file_pairs():
    """Match transcript CSV and mask TIFF files by sample name."""
    pairs = []

    for csv_path in sorted(
        TRANSCRIPT_DIR.glob("*.csv")
    ):
        mask_path = MASK_DIR / f"{csv_path.stem}.tif"

        if mask_path.exists():
            pairs.append((csv_path, mask_path))
        else:
            print(
                f"Warning: no mask found for "
                f"{csv_path.name}"
            )

    return pairs


def load_and_assign_transcripts(
    csv_path,
    mask_path,
    keep_genes,
):
    """Assign filtered transcripts to mask objects."""
    transcripts = pd.read_csv(csv_path)
    mask = tifffile.imread(mask_path)

    required_columns = {
        "feature_name",
        "x_location",
        "y_location",
    }

    missing_columns = required_columns.difference(
        transcripts.columns
    )

    if missing_columns:
        raise ValueError(
            f"{csv_path.name} is missing columns: "
            f"{sorted(missing_columns)}"
        )

    if mask.ndim != 2:
        raise ValueError(
            f"Expected a 2D mask for {mask_path.name}, "
            f"found shape {mask.shape}."
        )

    transcripts = transcripts[
        transcripts["feature_name"].astype(str).isin(
            keep_genes
        )
    ].copy()

    transcripts["x_location"] = pd.to_numeric(
        transcripts["x_location"],
        errors="coerce",
    )

    transcripts["y_location"] = pd.to_numeric(
        transcripts["y_location"],
        errors="coerce",
    )

    height, width = mask.shape

    valid_coordinates = (
        np.isfinite(transcripts["x_location"])
        & np.isfinite(transcripts["y_location"])
        & transcripts["x_location"].ge(0)
        & transcripts["x_location"].lt(width)
        & transcripts["y_location"].ge(0)
        & transcripts["y_location"].lt(height)
    )

    transcripts = transcripts.loc[
        valid_coordinates
    ].copy()

    transcripts["x_int"] = np.floor(
        transcripts["x_location"]
    ).astype(int)

    transcripts["y_int"] = np.floor(
        transcripts["y_location"]
    ).astype(int)

    transcripts["cell_id"] = mask[
        transcripts["y_int"].to_numpy(),
        transcripts["x_int"].to_numpy(),
    ]

    transcripts = transcripts[
        transcripts["cell_id"] != 0
    ].copy()

    if transcripts.empty:
        raise ValueError(
            f"No transcripts were assigned for "
            f"{csv_path.name}."
        )

    return transcripts, mask


def create_anndata(transcripts, mask):
    """Create an object-by-gene count matrix."""
    counts = pd.crosstab(
        transcripts["cell_id"],
        transcripts["feature_name"],
    )

    counts.index = counts.index.astype(str)

    cell_ids = np.unique(mask)
    cell_ids = cell_ids[cell_ids != 0]

    centroids = ndimage.center_of_mass(
        mask,
        labels=mask,
        index=cell_ids,
    )

    areas = np.bincount(
        mask.ravel().astype(np.int64)
    )[cell_ids]

    obs = pd.DataFrame({
        "cell_id": cell_ids.astype(str),
        "y_centroid": [
            centroid[0]
            for centroid in centroids
        ],
        "x_centroid": [
            centroid[1]
            for centroid in centroids
        ],
        "cell_area": areas,
    }).set_index("cell_id")

    obs = obs.loc[counts.index]

    adata = ad.AnnData(
        X=sp.csr_matrix(counts.to_numpy()),
        obs=obs,
    )

    adata.var_names = counts.columns.astype(str)

    # Napari-compatible row-column convention: [y, x].
    adata.obsm["spatial"] = obs[
        ["y_centroid", "x_centroid"]
    ].to_numpy()

    return adata


def apply_qc(adata):
    """Apply transcript-density and gene-detection filters."""
    sc.pp.calculate_qc_metrics(
        adata,
        inplace=True,
        percent_top=None,
    )

    adata.obs["transcript_density"] = (
        adata.obs["total_counts"]
        / adata.obs["cell_area"]
    )

    low, high = np.percentile(
        adata.obs["transcript_density"],
        DENSITY_PERCENTILES,
    )

    keep = adata.obs["transcript_density"].between(
        low,
        high,
    )

    adata = adata[keep].copy()

    sc.pp.filter_genes(
        adata,
        min_cells=MIN_CELLS_PER_GENE,
    )

    return adata


def preprocess_expression(adata):
    """Normalize expression and calculate PCA, neighbors and UMAP."""
    # Raw counts used by fibre-type annotation.
    adata.raw = adata.copy()

    # Raw counts used for seurat_v3 HVG selection.
    adata.layers["counts"] = adata.X.copy()

    sc.pp.normalize_total(
        adata,
        target_sum=10_000,
    )

    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(
        adata,
        layer="counts",
        n_top_genes=N_HVGS,
        flavor="seurat_v3",
        subset=False,
    )

    sc.pp.scale(
        adata,
        max_value=10,
    )

    sc.tl.pca(
        adata,
        n_comps=40,
        use_highly_variable=True,
        random_state=RANDOM_STATE,
    )

    sc.pp.neighbors(
        adata,
        n_neighbors=15,
        n_pcs=30,
        random_state=RANDOM_STATE,
    )

    sc.tl.umap(
        adata,
        random_state=RANDOM_STATE,
    )

    return adata


def preprocess_sample(
    csv_path,
    mask_path,
    keep_genes,
):
    """Preprocess one matched transcript and mask pair."""
    transcripts, mask = load_and_assign_transcripts(
        csv_path,
        mask_path,
        keep_genes,
    )

    adata = create_anndata(
        transcripts,
        mask,
    )

    n_objects_before_qc = adata.n_obs

    adata = apply_qc(adata)
    adata = preprocess_expression(adata)

    output_path = (
        OUTPUT_DIR
        / f"{csv_path.stem}.h5ad"
    )

    adata.write_h5ad(output_path)

    print(
        f"{csv_path.stem}: "
        f"{n_objects_before_qc} → {adata.n_obs} objects, "
        f"{adata.n_vars} genes"
    )
    print(f"Saved: {output_path}")


def main():
    keep_genes = load_genes()
    pairs = get_file_pairs()

    if not pairs:
        raise FileNotFoundError(
            "No matching transcript CSV and mask TIFF "
            "files were found."
        )

    for csv_path, mask_path in pairs:
        preprocess_sample(
            csv_path,
            mask_path,
            keep_genes,
        )


if __name__ == "__main__":
    main()
