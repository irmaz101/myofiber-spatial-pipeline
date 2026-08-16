# Benchmarking

This directory contains scripts used to evaluate four parts of the workflow:
segmentation, preprocessing filters, fibre-type annotation, and spatially
variable gene (SVG) detection.

The paths and main parameters are defined in the **Settings** section at the
top of each script. Update these values before running a script. Commands below
assume they are run from the repository root.

## Installation

The benchmarking scripts reuse the core analysis dependencies:

```bash
conda create -n myofiber-benchmarking python=3.10
conda activate myofiber-benchmarking
python -m pip install -r Benchmarking/requirements.txt
```

The synthetic-data notebook additionally uses R with `scDesign3`,
`SingleCellExperiment` and `reticulate`. These R dependencies are not installed
by the Python requirements file.

## Contents

| Folder | Purpose | Main outputs |
|---|---|---|
| `0_segmentation` | Compare a predicted instance mask with a manual mask and evaluate transcript assignment | Object-level metrics and a summary table |
| `1_preprocessing` | Compare total-count and transcript-density QC filters | Per-sample metrics, summary statistics, and diagnostic plots |
| `2_annotation` | Compare transcriptomic fibre labels with manual IHC labels | Class-wise precision, recall, F1, and matching tables |
| `3_SVG` | Compare SVG counts, method overlap, technical-replicate consistency, and performance on synthetic data | Summary tables, bootstrap rankings, and plots |

## 0. Segmentation benchmark

Script: `0_segmentation/benchmark_segmentation.py`

The script compares a predicted labelled mask with a manually annotated
reference mask. Manual and predicted objects are matched one-to-one based on
intersection over union (IoU), so each object can contribute to at most one
matched pair. Because the manual mask may contain only an annotated subset of
objects, unmatched predictions outside that subset are not treated as false
positives.

The following metrics are calculated:

- IoU and Dice coefficient;

- absolute and relative area error;

- detection recall at the selected IoU threshold;

- absolute and relative error in the number of assigned transcripts.

IoU, Dice, area-error and transcript-assignment summaries are calculated only
for one-to-one matched pairs meeting the IoU threshold. Precision and F1 are
not reported when the manual mask is not an exhaustive annotation.



## 1. Preprocessing benchmark

Script: `1_preprocessing/benchmark_qc_filters.py`

This script compares the objects removed by a total-transcript-count filter and
a transcript-density filter. It reports the overlap between the filters and
the correlations among object area, total counts, and transcript density.



## 2. Annotation benchmark

Script: `2_annotation/validate_fibre_annotation.py`

This script compares transcriptomic fibre-type predictions with manually
assigned IHC labels for type IIa and type IIb fibres. It calculates one-vs-rest
precision, recall, and F1 for both classes, as well as macro F1.



## 3. SVG benchmarks

The SVG benchmark contains two complementary analyses:

1. comparison of SVG sets detected in experimental datasets;
2. evaluation of method performance using synthetic spatial-expression data.

### Files

| File | Purpose |
|---|---|
| `0_summarize_svg_counts.py` | Count SVGs detected by each method in experimental samples |
| `1_calculate_svg_method_overlap.py` | Calculate pairwise overlap between methods |
| `2_benchmark_svg_technical_replicates.py` | Compare SVG sets between technical replicates |
| `3a_synthetic_datasets.ipynb` | Generate synthetic count matrices with scDesign3 |
| `3b_convert_simulation_to_anndata.py` | Convert one synthetic count matrix and its coordinates to AnnData |
| `4_compute_correlation.ipynb` | Calculate Kendall correlation with simulated spatial-signal strength |
| `5_compute_aupr.ipynb` | Calculate area under the precision-recall curve (AUPRC) |
| `6_bootstrap_svg_benchmark.py` | Bootstrap method rankings across synthetic datasets |

### Experimental samples and technical replicates

Place selected-gene CSV files in one directory per method, using identical
sample basenames across methods. The expected method directories are
`hotspot`, `moran_i`, `scgco`, `smash`, `somde`, `spagft`, `spatialde` and
`svgbit`.

Technical-replicate pairs are supplied in a user-created
`technical_replicates.csv`. Copy the neutral template before use:

```bash
cp Benchmarking/3_SVG/technical_replicates.example.csv \
   Benchmarking/3_SVG/technical_replicates.csv
```

The file has this structure:

```csv
sample_1,sample_2
reference_sample,replicate_sample
```

When this file is present, the count and method-overlap scripts exclude the
second sample in each pair from the independent-sample summaries. The technical
replicate script uses both samples to calculate overlap and Jaccard indices.

### Synthetic benchmark

`3a_synthetic_datasets.ipynb` processes one reference `.h5ad` file per run. Set
the neutral `input_file`, `output_id` and `n_genes` values at the top of the
notebook and rerun it for each required reference dataset. The simulation uses
raw counts from `adata.raw`, a Gaussian-process spatial mean model and signal
strengths from 0 to 1.

Convert each simulated count matrix to AnnData before running the eight SVG
methods:

```bash
python Benchmarking/3_SVG/3b_convert_simulation_to_anndata.py \
    --counts results_simulation/counts_100_sample.csv \
    --locations results_simulation/location_100_sample.csv \
    --output results_simulation/simulated_sample.h5ad
```

The transfer script preserves counts in both `adata.raw` and
`adata.layers["counts"]`, adds centroid coordinates, and retains the simulated
signal strength in the feature identifiers. After running every SVG method on
the synthetic AnnData file, place the complete method result tables in the
corresponding method directories. The Kendall and AUPRC notebooks recover the
base gene and signal strength from those identifiers and use only datasets
available for all eight methods.

The synthetic results characterize performance under the implemented
Gaussian-process signal-gradient scenario and should be interpreted within
that scope.

For the non-spatial component, the fitted mean matrix is permuted across
spatial objects while gene columns remain aligned. This preserves each gene's
fitted value distribution but breaks its association with location; it follows
the reference scDesign3 simulation strategy used for this benchmark.
