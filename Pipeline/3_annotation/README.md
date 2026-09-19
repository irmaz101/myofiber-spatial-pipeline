# Annotation

This directory contains the workflows used to annotate myofibre types and nucleus-level cell populations.

Myofibres and nuclei are analysed separately:

* Myofibre types are assigned using canonical myosin heavy-chain genes.
* Nuclear populations are identified using Leiden clustering and curated marker
  genes; decoupler/PanglaoDB signatures are used as supporting evidence.

## Files

| File                              | Purpose                                                                                  |
| --------------------------------- | ---------------------------------------------------------------------------------------- |
| `0_fibre_annotation.py`             | Annotates myofibre types using canonical *Myh* marker expression and k-means clustering. |
| `1_nuclei_qc.py`                    | Applies per-sample total-count QC to nucleus-level AnnData objects.                       |
| `1_integrate_nuclei_bbknn.py`       | Integrates sample-level nuclear datasets using BBKNN and performs Leiden clustering.     |
| `2_evaluate_nuclei_resolution.py`   | Compares Leiden resolutions on the saved BBKNN graph.                                    |
| `3_annotate_nuclear_populations.py` | Assigns manually curated biological labels to nuclear Leiden clusters.                   |
| `4_annotate_nuclei_decoupler.py`    | Calculates complementary PanglaoDB/decoupler signatures for nuclear clusters.            |

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

`0_fibre_annotation.py` expects sample-level AnnData files in the configured
`ANNDATA_DIR` (default: `anndata/`):

```text
anndata/
├── sample_01.h5ad
├── sample_02.h5ad
└── ...
```

Marker values are read from `adata.raw` when available; with the supplied
preprocessing script this contains filtered raw counts. If `adata.raw` is
absent, the script falls back to `adata.X` and prints a warning. Confirm that
the fallback matrix contains the intended expression representation before
continuing.

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

1. Raw-count expression of *Myh7*, *Myh2*, *Myh1*, and *Myh4* is extracted.
2. Marker expression is standardized independently within the sample.
3. Fibres are grouped into four clusters using k-means.
4. A dotplot of marker expression across the clusters is generated.
5. The user inspects the marker profiles and assigns each cluster a fibre-type label.
6. The final fibre-type labels are stored in
   `adata.obs["fiber_type_from_myh"]`.

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

Nucleus-level annotation consists of five separately executed steps:

1. per-sample total-count QC;
2. integration and clustering with BBKNN;
3. optional evaluation of Leiden resolutions;
4. manual annotation using curated markers;
5. complementary signature scoring using decoupler.

### Step 1: nucleus-level QC

Nuclear count matrices can be constructed using the same transcript-to-mask
assignment and count-matrix construction steps as for myofibres, with
nuclear masks as input. For this workflow, save the unfiltered nuclear
AnnData objects in `nuclei_anndata_unfiltered/` with raw transcript counts
in `adata.X`, before applying object-level QC, gene filtering, normalization
or scaling. Nuclear QC uses total transcript counts, whereas myofibre QC
uses transcript density (counts divided by object area). After nuclear QC,
retain raw counts in `adata.X`; normalization and log transformation are
performed later by `1_integrate_nuclei_bbknn.py`.

`1_nuclei_qc.py` reads unfiltered count matrices from
`nuclei_anndata_unfiltered/`. Within each sample, it retains nuclei between the
5th and 95th percentiles of `total_counts`, removes genes detected in fewer
than five retained nuclei, and writes `*.filtered.h5ad` files plus
`nuclei_qc_summary.csv` to `nuclei_anndata/`.

```bash
python Pipeline/3_annotation/1_nuclei_qc.py
```


### Step 2: BBKNN integration and clustering

`1_integrate_nuclei_bbknn.py` expects sample-level nuclear AnnData objects:

```text
nuclei_anndata/
├── sample_01.filtered.h5ad
├── sample_02.filtered.h5ad
└── ...
```

The input expression matrix must contain unnormalized transcript counts. Paths
are settings inside each script; the scripts do not form a command-line
pipeline. Run the integration script first, then point both annotation scripts
to its `results/combined_nuclei_bbknn.h5ad` output (or copy that file to the
configured input location).


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

Run from the repository root after configuring paths:

```bash
python Pipeline/3_annotation/1_integrate_nuclei_bbknn.py
```

The integration script writes `results/combined_nuclei_bbknn.h5ad`,
`results/cluster_counts_by_sample.csv`, and BBKNN UMAP plots. The integrated
embedding is stored in `adata.obsm["X_umap_bbknn"]`; the BBKNN graph is retained
in `adata.obsp` for subsequent resolution evaluation.

### Step 3: resolution evaluation

```bash
python Pipeline/3_annotation/2_evaluate_nuclei_resolution.py
```

The script evaluates resolutions 0.1--0.5, computes silhouette scores on a
reproducible sample of up to 10,000 PCA observations, and saves
`results/leiden_resolution_comparison.csv`. It does not overwrite the saved
integrated AnnData object.


## Step 4: manual annotation of nuclear populations

After setting `input_file` to the integrated object, run:

```bash
python Pipeline/3_annotation/3_annotate_nuclear_populations.py
```



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

## Step 5: decoupler-based nuclear annotation

A complementary marker-signature analysis is performed using decoupler and
canonical PanglaoDB markers. It is not an independent validation because it
uses the same expression data.

Run:

```bash
python Pipeline/3_annotation/4_annotate_nuclei_decoupler.py
```

The script:

1. retrieves canonical PanglaoDB markers;
2. retains markers supported for mouse;
3. matches mouse-supported PanglaoDB entries case-insensitively to the actual
   Xenium gene symbols (this is symbol matching, not complete orthologue
   mapping);
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
* myocytes;
* myoblasts;
* satellite cells;
* Schwann cells.

At least three marker genes must be available for a cell type to be evaluated.

### Decoupler outputs

The script saves:

```text
results/decoupler/
├── panglaodb_marker_counts_in_xenium.csv
├── panglaodb_mouse_supported_markers_used.csv
├── decoupler_scores_by_cluster.csv
├── decoupler_top5_predictions.csv
└── nuclei_panglaodb_support.h5ad
```

The automated label is stored in:

```python
adata.obs["panglaodb_top_signature"]
```

Decoupler predictions are used as supporting evidence for the manually curated annotations rather than as an independent ground truth.
