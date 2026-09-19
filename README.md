# Skeletal muscle Xenium analysis workflow

Code for fibre- and nucleus-resolved analysis of Xenium spatial transcriptomics
in skeletal muscle. The analyses cover segmentation, transcript assignment,
quality control, clustering, annotation and spatially variable gene (SVG)
detection, with benchmarks for segmentation, QC, annotation and SVG methods.

This repository accompanies a manuscript in preparation and may be updated
during review.

## Repository guide

| Stage | Analysis | Documentation |
| --- | --- | --- |
| Segmentation | Cellpose-SAM models for myofibres and nuclei | [Segmentation](Pipeline/0_segmentation/README.md) |
| Preprocessing | Transcript assignment, count matrices and fibre-density QC | [Preprocessing](Pipeline/1_preprocessing/README.md) |
| Clustering | Leiden benchmark and settings for SpaGCN, GraphST and BANKSY | [Clustering](Pipeline/2_clustering/README.md) |
| Annotation | Myh-based fibre types and BBKNN-integrated nuclear populations | [Annotation](Pipeline/3_annotation/README.md) |
| SVG detection | Eight methods applied to individual samples | [SVG methods](Pipeline/4_SVG/README.md) |
| Subcellular analysis | Notes on the exploratory Bento analysis | [Bento](Pipeline/5_subcellular_analysis/README.md) |
| Benchmarking | Segmentation, QC, IHC validation and SVG benchmarks | [Benchmarking](Benchmarking/README.md) |

## Running the analyses

Run commands from the repository root. Set the paths and parameters near the
top of each script, or use its command-line arguments where available. Stages
are run separately; there is no single command for the whole analysis.

For myofibres, preprocessing produces the sample-level AnnData files used
independently for Leiden clustering, fibre-type annotation and SVG detection.
Annotation and SVG detection use the preprocessed files, without passing
through the Leiden output.

Nuclei use the same transcript-assignment approach with nuclear masks. Their
count matrices are filtered by total counts and then integrated with BBKNN;
see [nuclear QC](Pipeline/3_annotation/README.md#step-1-nucleus-level-qc).

## Environments

Create separate environments for segmentation, core analysis, nuclear
annotation and each SVG method. The Python versions and dependencies are
listed in the corresponding requirements files. Some SVG packages require
incompatible dependency versions.

For preprocessing, Leiden clustering and fibre-type annotation:

```bash
conda create -n myofiber-core python=3.10
conda activate myofiber-core
python -m pip install -r Pipeline/1_preprocessing/requirements.txt
```

Installation commands for the other stages are in their linked READMEs.
Install GPU-enabled PyTorch for your hardware before setting up Cellpose.

## Inputs

Data and large intermediate files are not included.

| Analysis | Required input |
| --- | --- |
| Segmentation | Two-dimensional TIFF morphology images |
| Preprocessing | Transcript CSV files, matching labelled TIFF masks and a one-column `genes.txt` |
| Fibre clustering and annotation | Preprocessed sample-level `.h5ad` files |
| Nuclear QC and integration | Nuclear `.h5ad` files with raw counts in `.X` |
| SVG detection | Preprocessed `.h5ad` files with raw counts and centroid coordinates |
| Segmentation benchmark | Manual and predicted masks, plus a transcript table |
| Annotation validation | IHC labels and annotated AnnData files |
| Bento | A prepared SpatialData Zarr store |

Transcript tables need `feature_name`, `x_location` and `y_location` columns.
Coordinates must already match the mask pixel coordinate system. Transcript
and mask filenames are matched by basename, for example
`transcripts/sample_01.csv` and `masks_mf/sample_01.tif`.

## Included and external analyses

The repository includes the Leiden implementation and all eight SVG wrappers,
as well as synthetic-data generation and evaluation. The wrappers call the
original SVG packages; installation details and source links are in the
[SVG README](Pipeline/4_SVG/README.md).

SpaGCN, GraphST and BANKSY were run with external implementations. Their study
settings and source links are in the [clustering README](Pipeline/2_clustering/README.md).
The exploratory Bento analysis is documented separately; its study-specific
code is not included.
