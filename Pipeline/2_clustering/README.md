# Clustering

This folder contains the executable Leiden workflow used for sample-level
myofibre clustering. SpaGCN, GraphST and BANKSY were evaluated with their
published implementations, but their source code is not redistributed here.
The exact study settings and the benchmark repositories used as implementation
references are documented below.

## Included file

| File | Description |
| --- | --- |
| `0_leiden.py` | Runs Leiden clustering at resolutions 0.1-0.5, retains resolution 0.2 as the study setting, calculates silhouette scores and optional ARI, and saves clustered AnnData files. |

## Installation and usage

The included Leiden script uses the core dependencies from the preprocessing
stage:

```bash
conda create -n myofiber-core python=3.10
conda activate myofiber-core
python -m pip install -r Pipeline/1_preprocessing/requirements.txt
python Pipeline/2_clustering/0_leiden.py
```

GraphST, SpaGCN and BANKSY are not redistributed here. Install them according
to the linked original sources.

## Methods and study settings

All methods were applied independently to each sample-level AnnData object.

| Method | Study settings | Source |
| --- | --- | --- |
| Leiden | Scanpy; 1,000 Seurat-v3 HVGs; 40 PCs; 15 neighbours using the first 30 PCs; resolutions 0.1-0.5; selected resolution 0.2 | Included in this folder |
| SpaGCN | Raw counts normalized to 10,000 and log1p-transformed; no histology; `p = 0.5`; three domains; spatial refinement with `shape="square"` | [Benchmark implementation](https://github.com/BMILAB/Benchmark_ST_analysis/tree/master/1.Benchmark_SRT-main) and [SpaGCN](https://github.com/jianhuupenn/SpaGCN) |
| GraphST | 3,000 Seurat-v3 HVGs; normalization to 10,000, log1p and scaling; 600 epochs; Leiden used to obtain three domains; no spatial label refinement | [Benchmark implementation](https://github.com/BMILAB/Benchmark_ST_analysis/tree/master/1.Benchmark_SRT-main) and [GraphST](https://github.com/JinmiaoChenLab/GraphST) |
| BANKSY | 3,000 Seurat-v3 HVGs; normalization to 10,000, log1p and scaling; `lambda = 0.2`; `k_geom = 15`; scaled-Gaussian neighbour weights; `max_m = 1`; 20 PCs; Leiden used to obtain three domains | [Benchmark implementation](https://github.com/ZJUFanLab/SRTBenchmark/tree/main) and [BANKSY](https://prabhakarlab.github.io/Banksy/) |

The external methods may require separate environments because their
dependencies, particularly PyTorch and graph-learning libraries, can conflict.

## External-reference evaluation

The study compared clusters with three manually annotated reference categories:
type IIa, type IIb and SC71/BF-F3 double-negative fibres. Double-negative fibres
were retained as a third category but were not interpreted as independently
validated type IIx fibres. ARI was calculated only for object identifiers shared
by an AnnData file and the reference table. The reference table is not included
because its filenames and object identifiers are dataset-specific.
