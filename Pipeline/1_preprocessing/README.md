# Preprocessing

`preprocessing.py` assigns transcripts to myofibre masks and creates one
processed AnnData object per sample. Its output is used independently for
clustering, fibre-type annotation and SVG detection.

## Setup and inputs

Use the [core environment](../../README.md#environments). From the repository
root, run:

```bash
python Pipeline/1_preprocessing/preprocessing.py
```

Set these paths near the top of the script:

| Setting | Default | Contents |
| --- | --- | --- |
| `TRANSCRIPT_DIR` | `transcripts/` | Transcript CSV files |
| `MASK_DIR` | `masks_mf/` | Integer-labelled TIFF masks |
| `GENE_FILE` | `genes.txt` | One gene name per row, without a header |
| `OUTPUT_DIR` | `anndata_mf_density/` | AnnData files and QC comparisons |

Transcript and mask files must have matching basenames, such as
`sample_01.csv` and `sample_01.tif`. Each CSV needs `feature_name`,
`x_location` and `y_location` columns. Masks must be two-dimensional, with a
unique positive integer label per object and `0` for background.

Coordinates must already be in mask pixels, including any adjustments for
image cropping or resizing.

## Transcript assignment and count matrix

Only genes in `genes.txt` are retained. Invalid and out-of-bounds coordinates
are removed before converting the remaining coordinates to pixel indices:

```python
x_pixel = np.floor(x_location).astype(int)
y_pixel = np.floor(y_location).astype(int)
cell_id = mask[y_pixel, x_pixel]
```

Transcripts assigned to background are discarded. Counts are aggregated into
an object-by-gene matrix; objects with no assigned transcripts are omitted.

## Quality control

Transcript density is the number of assigned transcripts divided by mask
area. Within each sample, objects between the 5th and 95th density percentiles
are retained. Genes detected in fewer than five retained objects are removed.

The output comparison CSV records both the count-based and density-based
retain/remove decisions for the [QC benchmark](../../Benchmarking/README.md#preprocessing).
Only the density filter is applied to the saved myofibre AnnData.

## Expression processing

After QC, raw counts are copied to `adata.raw` and `adata.layers["counts"]`.
The expression matrix is normalized and reduced using these settings:

| Step | Setting |
| --- | --- |
| Normalization | 10,000 transcripts per object |
| Transformation | `log1p` |
| Highly variable genes | 1,000; `seurat_v3` on the raw-count layer |
| Scaling | Maximum scaled value of 10 |
| PCA | 40 components, using HVGs |
| Neighbour graph | 15 neighbours, using the first 30 PCs |
| UMAP | Calculated from the neighbour graph |
| Random state | 0 |

## Saved data

For `sample_01.csv`, the outputs are
`anndata_mf_density/sample_01.h5ad` and
`anndata_mf_density/sample_01_filtering_comparison.csv`.

| AnnData field | Contents |
| --- | --- |
| `obs_names` (index name `cell_id`) | Mask object IDs, stored as strings |
| `obs["x_centroid"]`, `obs["y_centroid"]` | Object centroids in pixels |
| `obs["cell_area"]` | Mask area in pixels² |
| `obs["total_counts"]` | Assigned transcript count used for object QC |
| `obs["transcript_density"]` | QC count divided by mask area |
| `obsm["spatial"]` | Centroids in `(y, x)` order |
| `raw`, `layers["counts"]` | Raw counts after object and gene filtering |
| `X` | Normalized, log-transformed and scaled expression |

PCA, the neighbour graph and UMAP are also saved. Use the explicitly named
centroid columns when coordinate order matters.

## Preparing nuclear data

For nuclei, reuse transcript assignment and count-matrix construction with
nuclear masks, then save the unfiltered counts before calling `apply_qc()` or
`preprocess_expression()`. Nuclear QC uses total counts, and normalization is
performed during BBKNN integration. See the
[nuclear preparation instructions](../3_annotation/README.md#step-1-nucleus-level-qc).
