# Annotation

Myofibre types are assigned from Myh expression. Nuclear populations are
identified after BBKNN integration using cluster markers, with PanglaoDB
signature scores as supporting evidence.

## Environments

Fibre annotation uses the [core environment](../../README.md#environments).
For the nuclear workflow:

```bash
conda create -n myofiber-annotation python=3.10
conda activate myofiber-annotation
python -m pip install -r Pipeline/3_annotation/requirements.txt
```

Run commands from the repository root. Paths are set near the top of each
script.

## Myofibre annotation

Set `ANNDATA_DIR` in `0_fibre_annotation.py` to the preprocessed fibre files
(e.g. `anndata_mf_density/`; the script default is `anndata/`), then run:

```bash
python Pipeline/3_annotation/0_fibre_annotation.py
```

Marker expression is read from `adata.raw`, which contains filtered raw counts
in the preprocessing output. If `.raw` is absent, the script warns and uses
`.X`; check the expression representation before using that fallback.

| Fibre type | Marker |
| --- | --- |
| Type I | *Myh7* |
| Type IIa | *Myh2* |
| Type IIx | *Myh1* |
| Type IIb | *Myh4* |

Gene names are matched case-insensitively. Marker expression is standardized
within each sample, then grouped by k-means (`k = 4`, `random_state = 0`,
`n_init = 50`). Inspect the saved Myh dotplot and enter a label for each cluster:
`type_1`, `type_2a`, `type_2x`, `type_2b` or `other`. Use `other` for mixed or
ambiguous marker profiles.

The default output directory is `annotated_anndata/`. Each sample produces a
`*_myh_dotplot.png` and `*_with_myh_fiber_types.h5ad`. The AnnData contains
`myh_cluster`, `fiber_type_from_myh` and the standardized marker matrix in
`obsm["X_myh_scaled"]`.

## Nuclear annotation

### Step 1: nucleus-level QC

Nuclear count matrices can be constructed using the same transcript-to-mask
assignment and count-matrix construction steps as for myofibres, with nuclear
masks as input. Save them in `nuclei_anndata_unfiltered/` with raw counts in
`adata.X`, before object QC, gene filtering, normalization or scaling.

Nuclear QC uses **total transcript counts**; myofibre QC uses **transcript
density** (counts divided by mask area). After nuclear QC, keep raw counts in
`.X`: normalization and log transformation take place during BBKNN integration.

```bash
python Pipeline/3_annotation/1_nuclei_qc.py
```

Within each sample, nuclei between the 5th and 95th count percentiles are
retained, followed by removal of genes detected in fewer than five nuclei.
Outputs are `nuclei_anndata/<sample>.filtered.h5ad` and
`nuclei_anndata/nuclei_qc_summary.csv`.

### Step 2: BBKNN integration

```bash
python Pipeline/3_annotation/1_integrate_nuclei_bbknn.py
```

Inputs are the filtered nuclear files in `nuclei_anndata/`, with unnormalized
counts in `.X`.

| Setting | Value |
| --- | --- |
| Batch variable | `sample` |
| HVGs | 2,000; `seurat_v3` on raw counts, accounting for sample |
| Normalization | 10,000 transcripts per nucleus, followed by `log1p` |
| PCA | 30 components from scaled HVGs |
| BBKNN neighbours per batch | 3 |
| Leiden resolution | 0.2 |
| Random state | 0 |
| Cluster key | `leiden_bbknn_r0.2` |

Results are saved in `results/combined_nuclei_bbknn.h5ad`, with cluster counts
in `results/cluster_counts_by_sample.csv` and UMAP plots by cluster and sample.
The integrated object keeps raw counts in `layers["counts"]`, normalized and
log-transformed full-gene expression in `.X` and `.raw`, the BBKNN graph in
`.obsp`, and UMAP coordinates in `obsm["X_umap_bbknn"]`.

### Step 3: compare resolutions (optional)

```bash
python Pipeline/3_annotation/2_evaluate_nuclei_resolution.py
```

Leiden resolutions 0.1–0.5 are evaluated on the saved BBKNN graph. Silhouette
scores use a reproducible sample of up to 10,000 nuclei in PCA space. Results
go to `results/leiden_resolution_comparison.csv`; the integrated AnnData is
not overwritten.

### Step 4: assign population labels

```bash
python Pipeline/3_annotation/3_annotate_nuclear_populations.py
```

Inspect cluster marker expression and differential expression before applying
the `cluster_labels` dictionary. The included mapping belongs to the study
dataset: cluster numbers can change with data, parameters or software versions.

| Population | Marker genes |
| --- | --- |
| Myonuclei | *Ttn*, *Mybpc1*, *Mybpc2*, *Myh1*, *Myh2*, *Myh4* |
| Endothelial | *Pecam1*, *Cdh5*, *Kdr*, *Esam* |
| FAPs | *Dcn*, *Col1a1*, *Col3a1*, *Pdgfra* |
| Vascular mural | *Myh11*, *Notch3*, *Pdgfrb* |
| Macrophages | *Mrc1*, *Csf1r*, *Adgre1*, *Cd163* |
| Schwann | *Cnp*, *Itgb4*, *Cd9* |
| Satellite / myogenic | *Pax7*, *Myf5*, *Fgfr4* |
| Lymphatic endothelial | *Lyve1*, *Ccl21a*, *Flt4* |

Labels are stored in `obs["nuclear_population"]` in
`results/annotated_nuclei.h5ad`. The script also saves cluster DE genes, a marker
dotplot, a population UMAP and a table of population counts and percentages.

### Step 5: PanglaoDB signatures

```bash
python Pipeline/3_annotation/4_annotate_nuclei_decoupler.py
```

Canonical PanglaoDB markers supported for mouse are matched case-insensitively
to Xenium gene symbols; this is symbol matching, not full orthologue mapping.
Cell types with at least three measured markers are scored by decoupler ULM
using the normalized `.raw` from the integrated analysis. Scores are averaged
within Leiden clusters, and the top five signatures are reported.

The candidate signatures cover endothelial cells, fibroblasts, macrophages,
pericytes, smooth muscle cells, myocytes, myoblasts, satellite cells and Schwann
cells. The highest-scoring signature is stored in
`obs["panglaodb_top_signature"]`.

Outputs in `results/decoupler/` include marker coverage, the marker list used,
cluster scores, top-five predictions and `nuclei_panglaodb_support.h5ad`.
These signatures support manual annotation; they use the same expression data
and are not an independent validation.
