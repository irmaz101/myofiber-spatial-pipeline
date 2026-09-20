# Benchmarking

Benchmarks cover segmentation, count- versus density-based QC, IHC validation
of fibre types, and SVG detection in experimental and synthetic datasets.

## Environment

```bash
conda create -n myofiber-benchmarking python=3.10
conda activate myofiber-benchmarking
python -m pip install -r Benchmarking/requirements.txt
```

The synthetic-data notebook also requires R with `scDesign3`,
`SingleCellExperiment` and `reticulate`; these are not installed by pip.
Run commands from the repository root and set input/output paths in each
script before running.

## Segmentation

`0_segmentation/benchmark_segmentation.py` matches predicted objects to a
manual reference one-to-one, using an IoU threshold of 0.5 by default.
It reports detection recall and, for matched objects, IoU, Dice, absolute and
relative area error, and differences in assigned transcript counts.

The manual reference may cover only a subset of objects. Unmatched predictions
are therefore not treated as false positives, and precision/F1 are not
reported. Overlap, area and transcript-count summaries include only matched
pairs meeting the IoU threshold.

## Preprocessing

`1_preprocessing/benchmark_qc_filters.py` reads the `*_filtering_comparison.csv`
files from preprocessing. It compares the objects removed by count and density
filters, reports their overlap and the correlations between area, counts and
density, and saves per-sample diagnostic plots.

## Fibre annotation

`2_annotation/validate_fibre_annotation.py` compares `fiber_type_from_myh`
labels with IHC annotations. The Excel input has one sheet per sample, with
object IDs in the first column and `2a`/`2b` columns. A `1` denotes a positive
label; `0` or a blank denotes a negative label for that marker.

The script matches fibre IDs and calculates one-vs-rest precision, recall and
F1 for IIa and IIb, plus macro F1. Double-negative fibres remain negative
observations in both evaluations. Outputs include matched fibres, a confusion
matrix, reference counts and per-sample and summary metrics.

## SVG benchmarks

| File in `3_SVG/` | Analysis |
| --- | --- |
| `0_summarize_svg_counts.py` | SVG counts per method and sample |
| `1_calculate_svg_method_overlap.py` | Pairwise overlap coefficient between methods |
| `2_benchmark_svg_technical_replicates.py` | Overlap coefficient and Jaccard index between adjacent-section pairs |
| `3a_synthetic_datasets.ipynb` | scDesign3 simulations |
| `3b_convert_simulation_to_anndata.py` | Convert simulated counts and coordinates to AnnData |
| `4_compute_correlation.ipynb` | Mean gene-wise Kendall correlation with signal strength |
| `5_compute_aupr.ipynb` | Average precision for the simulated signal classes |
| `6_bootstrap_svg_benchmark.py` | Bootstrap method rankings across datasets |

### Experimental data

Use the selected-gene CSVs from each SVG wrapper's `--output`. Set
`results_dir` to their parent directory, containing `hotspot`, `moran_i`,
`scgco`, `smash`, `somde`, `spagft`, `spatialde` and `svgbit`. Sample basenames
must match across methods.

Create the technical-replicate table from the example:

```bash
cp Benchmarking/3_SVG/technical_replicates.example.csv \
   Benchmarking/3_SVG/technical_replicates.csv
```

Replace the example names with the actual pairs:

```csv
sample_1,sample_2
reference_sample,replicate_sample
```

Set `replicate_file` in each relevant script to this file's path. When it is
present, the count and method-overlap summaries exclude `sample_2` from each
pair. The technical-replicate script requires the table and compares both
members of each pair.

### Synthetic data

In `3a_synthetic_datasets.ipynb`, set `input_dir`, `input_file`, `output_id` and
`n_genes` for one reference sample. The reference AnnData must contain raw
counts in `.raw`. The notebook fits Gaussian-process spatial mean models and
simulates a signal gradient from 0 to 1.

For the non-spatial component, fitted means are shuffled across spatial
objects while gene columns remain aligned. This preserves each gene's fitted
value distribution while breaking its spatial association. The benchmark
therefore evaluates this Gaussian-process signal-gradient scenario.

Convert the count and location files before running the SVG wrappers:

```bash
python Benchmarking/3_SVG/3b_convert_simulation_to_anndata.py \
    --counts results_simulation/counts_100_sample.csv \
    --locations results_simulation/location_100_sample.csv \
    --output results_simulation/simulated_sample.h5ad
```

The conversion preserves counts in `.raw` and `layers["counts"]`, adds centroid
coordinates, and retains the gene and signal strength in feature identifiers.

Run all eight SVG methods with `--all-output`. Place these complete result
tables in method directories under a common parent and set `input_dir` in the
Kendall and AUPRC notebooks to that parent. Both notebooks use datasets
available for all eight methods and parse signal strength from feature names.
For average precision (reported as AUPRC), positives have signal strength
**greater than 0.4**.

The notebooks save `results/kendall_results.csv` and `results/aupr.csv`.
`6_bootstrap_svg_benchmark.py` resamples complete datasets with replacement
(100 iterations, seed 42), then records the best and second-best mean method
scores. Outputs include bootstrap scores, winning proportions and excluded
datasets.
