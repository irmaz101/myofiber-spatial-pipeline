# Preprocessing

This directory contains the preprocessing workflow used to assign Xenium transcripts to segmented objects and create sample-level AnnData objects for downstream analysis.

The same workflow can be applied to myofibres and nuclei by changing the input mask directory and output directory.

## Files

| File               | Purpose                                                                                                                                                                    |
| ------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `preprocessing.py` | Assigns transcripts to segmentation masks, constructs object-by-gene count matrices, performs quality control, normalizes expression, and saves processed AnnData objects. |

## Installation

Create and activate the core analysis environment from the repository root:

```bash
conda create -n myofiber-core python=3.10
conda activate myofiber-core
python -m pip install -r Pipeline/1_preprocessing/requirements.txt
```

The same environment can be used for the provided Leiden clustering and
myofibre-type annotation scripts.

## Input files

The script requires:

1. Xenium transcript tables in CSV format.
2. Corresponding integer-labelled segmentation masks in TIFF format.

Transcript and mask files are matched by their basename:

```text
transcripts/sample_01.csv
masks_mf/sample_01.tif
```

The transcript CSV must contain the following columns:

| Column         | Description                                            |
| -------------- | ------------------------------------------------------ |
| `feature_name` | Gene assigned to the transcript.                       |
| `x_location`   | Transcript x coordinate in the mask coordinate system. |
| `y_location`   | Transcript y coordinate in the mask coordinate system. |

The segmentation mask must:

* have the same coordinate system as the transcript table;
* contain non-negative integer labels;
* use `0` for background;
* assign a unique positive integer label to each segmented object.

## Input and output directories

The default directories are defined in `preprocessing.py`:

```python
transcript_dir = "transcripts"
mask_dir = "masks_mf"
output_dir = "anndata_mf_density"
```

These values can be changed to process different object types.

For example, nucleus-level preprocessing can use:

```python
transcript_dir = "transcripts"
mask_dir = "masks_nuclei"
output_dir = "anndata_nuclei_density"
```

## Running the preprocessing pipeline

Run the script from the repository root:

```bash
python Pipeline/1_preprocessing/preprocessing.py
```

The script identifies transcript CSV files and matches each file to a TIFF mask with the same basename.

## Transcript assignment

Transcript coordinates are converted to integer pixel indices:

```python
x_pixel = x_location.astype(int)
y_pixel = y_location.astype(int)
```

Transcripts with invalid coordinates or coordinates outside the mask boundaries are excluded.

Each remaining transcript is assigned to the segmentation label at its pixel position:

```python
cell_id = mask[y_pixel, x_pixel]
```

Transcripts assigned to background label `0` are excluded.

## Count-matrix construction

An object-by-gene count matrix is constructed from the assigned transcripts.

* Rows represent segmented myofibres or nuclei.
* Columns represent genes.
* Values represent the number of transcripts assigned to each object.

Objects without any assigned transcripts are not included in the resulting AnnData object.

## Spatial information

The script calculates the centroid and area of each segmented object.

The following information is stored in `adata.obs`:

| Field                | Description                                            |
| -------------------- | ------------------------------------------------------ |
| `cell_id`            | Integer label of the segmented object.                 |
| `x_centroid`         | Object centroid along the x axis.                      |
| `y_centroid`         | Object centroid along the y axis.                      |
| `cell_area`          | Object area in pixels.                                 |
| `total_counts`       | Total number of assigned transcripts.                  |
| `transcript_density` | Number of assigned transcripts divided by object area. |

Spatial coordinates are also stored in:

```python
adata.obsm["spatial"]
```

The original analysis stores these coordinates in `(y, x)` order. The explicitly named `x_centroid` and `y_centroid` columns should therefore be used when coordinate order is important.

## Quality control

Transcript density is calculated as:

```text
transcript density = total assigned transcripts / object area
```

Two quality-control strategies are evaluated independently:

1. filtering based on total transcript counts;
2. filtering based on transcript density.

For each sample, the 5th and 95th percentiles are calculated separately for total counts and transcript density.

The script records whether each object would be retained by either approach:

```text
keep_count_filter
keep_density_filter
```

These object-level results are saved as:

```text
<sample>_filtering_comparison.csv
```

The transcript-density filter is used for downstream analysis. Objects between the sample-specific 5th and 95th percentiles of transcript density are retained.

Genes detected in fewer than five retained objects are removed.

## Normalization and dimensionality reduction

After quality control, the following preprocessing steps are performed:

1. Raw transcript counts are stored in:

   ```python
   adata.layers["counts"]
   ```

2. The 1,000 most highly variable genes are identified using the `seurat_v3` method.

3. Counts are normalized to 10,000 transcripts per object:

   ```python
   sc.pp.normalize_total(adata, target_sum=10_000)
   ```

4. Expression values are log-transformed:

   ```python
   sc.pp.log1p(adata)
   ```

5. Normalized, unscaled expression is preserved in:

   ```python
   adata.raw
   ```

6. The working expression matrix is scaled with a maximum value of 10.

7. Principal component analysis is calculated using 40 components.

8. A neighbourhood graph is constructed using 15 neighbours and the first 30 principal components.

9. A UMAP representation is calculated.

## Main study settings

| Parameter                | Setting                       |
| ------------------------ | ----------------------------- |
| Object-level QC variable | Transcript density            |
| Lower QC percentile      | 5th percentile                |
| Upper QC percentile      | 95th percentile               |
| Minimum gene detection   | 5 objects                     |
| Highly variable genes    | 1,000                         |
| HVG method               | `seurat_v3`                   |
| Normalization target     | 10,000 transcripts per object |
| Transformation           | `log1p`                       |
| Maximum scaled value     | 10                            |
| PCA components           | 40                            |
| Nearest neighbours       | 15                            |
| PCs used for neighbours  | 30                            |

## Outputs

For each sample, the script produces:

```text
anndata_mf_density/
├── sample_01.h5ad
└── sample_01_filtering_comparison.csv
```

The `.h5ad` file contains the processed expression matrix, spatial coordinates, object-level metadata, raw counts, normalized expression, PCA representation, neighbourhood graph, and UMAP coordinates.

The filtering-comparison CSV contains object-level measurements and the decisions produced by the count- and density-based QC strategies.
