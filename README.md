# Skeletal muscle Xenium analysis workflow

This repository contains the code used to analyse imaging-based spatial
transcriptomics data from skeletal muscle. The workflow treats multinucleated
myofibres and nuclei as separate analytical units and covers segmentation,
transcript assignment, preprocessing, clustering, annotation, spatially
variable gene detection and exploratory subcellular analysis.

The repository accompanies a manuscript currently in preparation. Code and
documentation may be updated during internal and peer review.

## Workflow

| Step                 | Description                                                                                                             |
| -------------------- | ----------------------------------------------------------------------------------------------------------------------- |
| Segmentation         | Segment myofibres and nuclei from Xenium morphology images using Cellpose-SAM.                                          |
| Preprocessing        | Assign transcripts to labelled masks, calculate object-level measurements and apply transcript-density quality control. |
| Clustering           | Compare Leiden, GraphST, SpaGCN and BANKSY for fibre-level clustering.                                                  |
| Annotation           | Annotate myofibre types using canonical *Myh* markers and nuclear populations using curated markers and decoupler.      |
| SVG detection        | Detect spatially variable genes using eight complementary methods.                                                      |
| Subcellular analysis | Explore subcellular transcript localization using Bento.                                                                |


## Input data

Raw Xenium data and large intermediate files are not included in this
repository. The main expected inputs are:

| Analysis                  | Input                                                                    |
| ------------------------- | ------------------------------------------------------------------------ |
| Segmentation              | Two-dimensional TIFF morphology images                                   |
| Preprocessing             | Transcript CSV files and matching integer-labelled TIFF masks            |
| Clustering and annotation | Sample-level AnnData (`.h5ad`) files                                     |
| SVG detection             | Preprocessed AnnData files containing expression and spatial coordinates |
| Segmentation benchmarking | Manual and predicted masks plus the corresponding transcript table       |
| Annotation validation     | Manual IHC annotations and annotated AnnData files                       |
| Bento                     | A prepared SpatialData Zarr store                                        |

Transcript tables must contain `feature_name`, `x_location` and `y_location`.
Transcript and mask files are matched by basename, for example:

```text
transcripts/sample_01.csv
masks_mf/sample_01.tif
```

## Running the workflow

The main stages are run in the following order:

```text
segmentation → preprocessing → clustering → annotation → SVG detection
```

The subcellular Bento analysis is an independent exploratory analysis. See the
README within each pipeline folder for the exact commands and parameters.

### 1. Segmentation

Scripts in `Pipeline/0_segmentation/` extract training patches, open them for
manual annotation in Napari, train custom Cellpose-SAM models and segment full
images. The same workflow can be configured for myofibres or nuclei.

### 2. Preprocessing

`Pipeline/1_preprocessing/preprocessing.py` assigns transcripts to mask labels,
constructs object-by-gene count matrices and stores object area and centroid
coordinates. Objects between the sample-specific 5th and 95th percentiles of
transcript density are retained, and genes detected in fewer than five objects
are removed.

Expression is normalized to 10,000 transcripts per object, log-transformed and
used to calculate highly variable genes, PCA, a neighbourhood graph and UMAP.

### 3. Clustering

Leiden, GraphST, SpaGCN and BANKSY were evaluated for fibre-level clustering.
GraphST, SpaGCN and BANKSY are external methods and should be installed and
cited from their original repositories. Study-specific settings are documented
in `Pipeline/2_clustering/README.md`.

### 4. Annotation

Myofibre types are assigned using expression of *Myh7*, *Myh2*, *Myh1* and
*Myh4*. Nucleus-level datasets are integrated using BBKNN and annotated using
curated marker genes. Decoupler and PanglaoDB markers provide an additional
annotation for comparison.

### 5. Spatially variable genes

The following SVG methods are included:

- Moran's I;
- SpatialDE;
- Hotspot;
- scGCO;
- SOMDE;
- SMASH;
- SpaGFT;
- SVGbit.

The scripts in `Pipeline/4_SVG/` process one AnnData file at a time. They can
save both the selected SVGs and a complete gene-level table containing an
`is_svg` column. The SVG-only files should be used for benchmarking.

### 6. Subcellular analysis

Bento was used for an exploratory proof-of-concept analysis of subcellular
transcript localization. Bento is an external package and should be installed
and cited from its original repository.

## Benchmarking

The `Benchmarking/` directory contains study-specific evaluation scripts for:

- segmentation performance and transcript-assignment error;
- comparison of count- and density-based preprocessing filters;
- validation of myofibre annotation against manual IHC labels;
- SVG counts, overlap between methods, and technical-replicate consistency.

The clustering evaluation and synthetic SVG-data generation are described in
the manuscript but are not provided as self-contained workflows in the current
repository version.


