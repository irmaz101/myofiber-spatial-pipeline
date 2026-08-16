# Skeletal muscle Xenium analysis workflow

This repository contains the code used to analyse imaging-based spatial
transcriptomics data from skeletal muscle. The workflow treats multinucleated
myofibres and nuclei as separate analytical units and covers segmentation,
transcript assignment, preprocessing, clustering, annotation, spatially
variable gene detection and exploratory subcellular analysis.

The repository accompanies a manuscript currently in preparation. Code and
documentation may be updated during internal and peer review.

## Workflow

| Step | Description | Documentation |
| --- | --- | --- |
| Segmentation | Segment myofibres and nuclei from Xenium morphology images using Cellpose-SAM. | [Segmentation README](Pipeline/0_segmentation/README.md) |
| Preprocessing | Assign transcripts to labelled masks, calculate object-level measurements and apply transcript-density quality control. | [Preprocessing README](Pipeline/1_preprocessing/README.md) |
| Clustering | Run the included Leiden workflow; GraphST, SpaGCN and BANKSY are documented as external methods. | [Clustering README](Pipeline/2_clustering/README.md) |
| Annotation | Annotate myofibre types using canonical *Myh* markers and nuclear populations using curated markers and decoupler. | [Annotation README](Pipeline/3_annotation/README.md) |
| SVG detection | Detect spatially variable genes using eight complementary methods. | [SVG README](Pipeline/4_SVG/README.md) |
| Subcellular analysis | Explore subcellular transcript localization using Bento. | [Bento documentation](Pipeline/5_subcellular_analysis/README.md) |
| Benchmarking | Evaluate segmentation, preprocessing, annotation and SVG-detection performance. | [Benchmarking README](Benchmarking/README.md) |

## Installation

The repository is organized as a sequence of independent stages. Each stage
must be run separately after updating the paths and settings in its scripts;
there is no single command that executes the complete workflow. Outputs from
one stage serve as inputs to the next stage.

Several SVG packages have incompatible dependency requirements. We therefore
recommend creating a separate conda environment for each major stage and for
each SVG method. Stage-specific requirements are provided in the corresponding
directories.

For example, create the preprocessing environment from the repository root:

```bash
conda create -n myofiber-preprocessing python=3.10
conda activate myofiber-preprocessing
python -m pip install -r Pipeline/1_preprocessing/requirements.txt
```

The provided Leiden clustering and myofibre-annotation scripts can use this
same core environment. Nuclear integration and decoupler annotation use:

```bash
conda create -n myofiber-annotation python=3.10
conda activate myofiber-annotation
python -m pip install -r Pipeline/3_annotation/requirements.txt
```

For SVG methods, use the method-specific files under
`Pipeline/4_SVG/requirements/`. For example:

```bash
conda create -n myofiber-morans python=3.10
conda activate myofiber-morans
python -m pip install -r Pipeline/4_SVG/requirements/morans.txt
```

The comments at the top of each requirements file state the Python version
used for that analysis. GPU-enabled PyTorch should be installed according to
the user's hardware and the official PyTorch instructions.


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
linked documentation for the exact commands and parameters.

### 1. Segmentation

Scripts in `Pipeline/0_segmentation/` extract training patches, open them for
manual annotation in Napari, train custom Cellpose-SAM models and segment full
images. The same workflow can be configured for myofibres or nuclei.

### 2. Preprocessing

[`Pipeline/1_preprocessing/preprocessing.py`](Pipeline/1_preprocessing/preprocessing.py)
assigns transcripts to mask labels,
constructs object-by-gene count matrices and stores object area and centroid
coordinates. Objects between the sample-specific 5th and 95th percentiles of
transcript density are retained, and genes detected in fewer than five objects
are removed.

Expression is normalized to 10,000 transcripts per object, log-transformed and
used to calculate highly variable genes, PCA, a neighbourhood graph and UMAP.

### 3. Clustering

Leiden, GraphST, SpaGCN and BANKSY were evaluated for fibre-level clustering.
The executable Leiden workflow is included. GraphST, SpaGCN and BANKSY are
external methods and should be installed and cited from their original
repositories. Study-specific settings and source links are documented in the
[clustering README](Pipeline/2_clustering/README.md).

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
transcript localization. Study-specific code is not redistributed; the
[subcellular-analysis README](Pipeline/5_subcellular_analysis/README.md) links
to the current Bento tutorials.

## Benchmarking

The [`Benchmarking/`](Benchmarking/README.md) directory contains
study-specific evaluation scripts for:

- segmentation performance and transcript-assignment error;
- comparison of count- and density-based preprocessing filters;
- validation of myofibre annotation against manual IHC labels;
- SVG counts, overlap between methods, and technical-replicate consistency.

The synthetic SVG-data generation and downstream benchmark workflows are
included. Clustering implementations for SpaGCN, GraphST and BANKSY are not
redistributed; their exact study settings and implementation references are
recorded in the [clustering README](Pipeline/2_clustering/README.md).
