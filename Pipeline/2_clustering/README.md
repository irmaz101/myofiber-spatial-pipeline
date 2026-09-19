# Clustering

`0_leiden.py` evaluates Leiden resolutions separately for each myofibre sample.
SpaGCN, GraphST and BANKSY were evaluated with external implementations; their
settings and sources are listed below.

## Run Leiden

Use the [core environment](../../README.md#environments), set `INPUT_DIR`,
`GROUND_TRUTH_FILE` and `OUTPUT_DIR`, then run from the repository root:

```bash
python Pipeline/2_clustering/0_leiden.py
```

Inputs are the preprocessed sample-level AnnData files with raw counts in
`adata.raw`. Preprocessing is reconstructed for this benchmark. Resolutions
0.1–0.5 are evaluated, and resolution 0.2 is saved as `adata.obs["leiden"]`.
Silhouette scores use the first 30 PCs. ARI is calculated when the manual
reference file is available.

Each sample produces a `*_leiden_resolutions.h5ad` file. Summary CSVs contain
scores by sample and resolution, cluster sizes and reference-matching counts.
The saved `.raw` contains normalized, log-transformed expression. These files
are clustering results; fibre annotation and SVG detection use the original
preprocessed files.

## Method settings

All methods were run independently on each sample.

| Method | Study settings | Source |
| --- | --- | --- |
| Leiden | 1,000 Seurat-v3 HVGs; 40 PCs; 15 neighbours using the first 30 PCs; resolutions 0.1–0.5; selected resolution 0.2 | `0_leiden.py` |
| SpaGCN | Counts normalized to 10,000 and log1p-transformed; no histology; `p = 0.5`; three domains; spatial refinement with `shape="square"` | [Benchmark implementation](https://github.com/BMILAB/Benchmark_ST_analysis/tree/master/1.Benchmark_SRT-main), [SpaGCN](https://github.com/jianhuupenn/SpaGCN) |
| GraphST | 3,000 Seurat-v3 HVGs; normalization to 10,000, log1p and scaling; 600 epochs; Leiden used to obtain three domains; no spatial label refinement | [Benchmark implementation](https://github.com/BMILAB/Benchmark_ST_analysis/tree/master/1.Benchmark_SRT-main), [GraphST](https://github.com/JinmiaoChenLab/GraphST) |
| BANKSY | 3,000 Seurat-v3 HVGs; normalization to 10,000, log1p and scaling; `lambda = 0.2`; `k_geom = 15`; scaled-Gaussian neighbour weights; `max_m = 1`; 20 PCs; Leiden used to obtain three domains | [Benchmark implementation](https://github.com/ZJUFanLab/SRTBenchmark/tree/main), [BANKSY](https://prabhakarlab.github.io/Banksy/) |

Install external methods in separate environments following their source
repositories, particularly where PyTorch or graph-library requirements differ.

## Manual reference

The reference CSV needs `sample`, `cell_id` and `fiber_type` columns. The three
categories are type IIa, type IIb and SC71/BF-F3 double-negative fibres.
Double-negative fibres form a separate reference category; they are not
independently validated type IIx fibres.

Only object IDs shared by the AnnData and reference table enter the ARI
calculation. The study reference table is not included in this repository.
