# Benchmarking

This directory contains scripts used to evaluate four parts of the workflow:
segmentation, preprocessing filters, fibre-type annotation, and spatially
variable gene (SVG) detection.

The paths and main parameters are defined in the **Settings** section at the
top of each script. Update these values before running a script. Commands below
assume they are run from the repository root.

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
reference mask. Each manual object is paired with the overlapping predicted
object that has the highest intersection over union (IoU). This is a
best-overlap match for each manual object, not a global one-to-one assignment.

The following metrics are calculated:

- IoU and Dice coefficient;
- absolute and relative area error;
- precision, recall, and F1 at the selected IoU threshold;
- oversegmentation and undersegmentation rates;
- absolute and relative error in the number of assigned transcripts.

.

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
| `3_synthetic_dataset.ipynb` | Generate synthetic count matrices with scDesign3 |
| `4_compute_correlation.ipynb` | Calculate Kendall correlation with simulated spatial-signal strength |
| `5_compute_aupr.ipynb` | Calculate area under the precision-recall curve (AUPRC) |
| `6_bootstrap_svg_benchmark.py` | Bootstrap method rankings across synthetic datasets |




