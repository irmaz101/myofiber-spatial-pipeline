# Annotation

This directory contains the workflows used to annotate myofibre types and nucleus-level cell populations.

Myofibres and nuclei are analysed separately:

* Myofibre types are assigned using canonical myosin heavy-chain genes.
* Nuclear populations are identified using Leiden clustering, curated marker genes, and an independent decoupler-based annotation.

## Files

| File                              | Purpose                                                                                  |
| --------------------------------- | ---------------------------------------------------------------------------------------- |
| `0_fibre_annotation.py`             | Annotates myofibre types using canonical *Myh* marker expression and k-means clustering. |
| `1_integrate_nuclei_bbknn.py`       | Integrates sample-level nuclear datasets using BBKNN and performs Leiden clustering.     |
| `3_annotate_nuclear_populations.py` | Assigns manually curated biological labels to nuclear Leiden clusters.                   |
| `4_annotate_nuclei_decoupler.py`    | Performs an independent marker-based nuclear annotation using decoupler and PanglaoDB.   |

## Installation

The myofibre annotation script can use the core preprocessing environment.
For the complete nuclear workflow, create the annotation environment:

```bash
conda create -n myofiber-annotation python=3.10
conda activate myofiber-annotation
python -m pip install -r Pipeline/3_annotation/requirements.txt
```

## Myofibre-type annotation

### Input

`0_fibre_annotation.py` expects sample-level AnnData files produced by the preprocessing pipeline:

```text
anndata_mf_density/
├── sample_01.h5ad
├── sample_02.h5ad
└── ...
```

Normalized, unscaled expression is read from `adata.raw` when available.

### Marker genes

Myofibre types are identified using canonical myosin heavy-chain markers:

| Fibre type | Marker |
| ---------- | ------ |
| Type I     | *Myh7* |
| Type IIa   | *Myh2* |
| Type IIx   | *Myh1* |
| Type IIb   | *Myh4* |

Gene matching is case-insensitive. Missing markers are reported in the terminal.

### Method

For each sample:

1. Expression of *Myh7*, *Myh2*, *Myh1*, and *Myh4* is extracted.
2. Marker expression is standardized independently within the sample.
3. Fibres are grouped into four clusters using k-means.
4. A dotplot of marker expression across the clusters is generated.
5. The user inspects the marker profiles and assigns each cluster a fibre-type label.
6. The final fibre-type labels are stored in `adata.obs["fiber_type"]`.

The allowed labels are:

```text
type_1
type_2a
type_2x
type_2b
other
```

The `other` label should be used when a cluster has mixed or ambiguous marker expression.

### Main settings

| Parameter                         | Setting    |
| --------------------------------- | ---------- |
| Number of k-means clusters        | 4          |
| Random state                      | 0          |
| Number of k-means initializations | 50         |
| Marker scaling                    | Per sample |
| Ambiguous label                   | `other`    |



## Nuclear-population annotation

Nucleus-level annotation consists of three steps:

1. integration and clustering with BBKNN;
2. manual annotation using curated markers;
3. independent marker-based annotation using decoupler.

### Step 1: BBKNN integration and clustering

`integrate_nuclei_bbknn.py` expects sample-level nuclear AnnData objects:

```text
nuclei_anndata/
├── sample_01.h5ad
├── sample_02.h5ad
└── ...
```

The input expression matrix should contain unnormalized transcript counts.


### Nuclear integration settings

| Parameter                          | Setting             |
| ---------------------------------- | ------------------- |
| Batch variable                     | `sample`            |
| Highly variable genes              | 2,000               |
| HVG method                         | `seurat_v3`         |
| Normalization target               | 10,000              |
| PCA components                     | 30                  |
| BBKNN neighbours within each batch | 3                   |
| Leiden resolution                  | 0.2                 |
| Cluster column                     | `leiden_bbknn_r0.2` |


## Manual annotation of nuclear populations



### Nuclear marker sets

| Population                   | Example markers                                   |
| ---------------------------- | ------------------------------------------------- |
| Myonuclei                    | *Ttn*, *Mybpc1*, *Mybpc2*, *Myh1*, *Myh2*, *Myh4* |
| Endothelial cells            | *Pecam1*, *Cdh5*, *Kdr*, *Esam*                   |
| Fibro-adipogenic progenitors | *Dcn*, *Col1a1*, *Col3a1*, *Pdgfra*               |
| Vascular mural cells         | *Myh11*, *Notch3*, *Pdgfrb*                       |
| Macrophages                  | *Mrc1*, *Csf1r*, *Adgre1*, *Cd163*                |
| Schwann cells                | *Cnp*, *Itgb4*, *Cd9*                             |
| Satellite/myogenic cells     | *Pax7*, *Myf5*, *Fgfr4*                           |
| Lymphatic endothelial cells  | *Lyve1*, *Ccl21a*, *Flt4*                         |

The curated labels are stored in:

```python
adata.obs["nuclear_population"]
```

### Important note about cluster labels

The `cluster_labels` dictionary in the script is specific to the integrated dataset used in the study.

Leiden cluster numbers are not biologically fixed and may change if the input data, software versions, random seeds, or clustering parameters are changed. The cluster-specific marker expression and differentially expressed genes must therefore be inspected before applying these labels to a new analysis.

### Manual annotation outputs

The script saves:

```text
results/
├── nuclear_cluster_de_genes.csv
├── nuclear_marker_dotplot.png
├── nuclear_population_summary.csv
└── annotated_nuclei.h5ad
```

The population summary contains the number and percentage of nuclei assigned to each population.

## Decoupler-based nuclear annotation

An independent marker-based annotation is performed using decoupler and canonical PanglaoDB markers.

Run:

```bash
python Pipeline/3_annotation/4_annotate_nuclei_decoupler.py
```

The script:

1. retrieves canonical PanglaoDB markers;
2. retains markers supported for mouse;
3. converts human-style gene symbols to mouse-style capitalization;
4. keeps markers present in the nuclear dataset;
5. calculates per-nucleus cell-type scores using the univariate linear model method;
6. averages the scores within each Leiden cluster;
7. reports the five highest-scoring predictions for each cluster;
8. assigns the highest-scoring cell type as the decoupler label.

The evaluated populations include:

* endothelial cells;
* fibroblasts;
* macrophages;
* pericytes;
* smooth muscle cells;
* myoblasts;
* Schwann cells.

At least three marker genes must be available for a cell type to be evaluated.

### Decoupler outputs

The script saves:

```text
results/
├── panglaodb_markers_used.csv
├── decoupler_scores_by_cluster.csv
├── decoupler_top5_predictions.csv
└── nuclei_decoupler_annotated.h5ad
```

The automated label is stored in:

```python
adata.obs["decoupler_label"]
```

Decoupler predictions are used as supporting evidence for the manually curated annotations rather than as an independent ground truth.


