# Clustering

This folder contains the Leiden clustering workflow used for sample-level
myofibre analysis. GraphST, SpaGCN and BANKSY were also evaluated in the study,
but their implementations are not redistributed in this repository. Those
methods should be installed and cited from their original sources.

## Included file

| File | Description |
| --- | --- |
| `leiden.py` | Runs Leiden clustering, calculates silhouette scores and adjusted Rand index (ARI), creates spatial and UMAP plots, and saves clustered AnnData files. |

## Methods and study settings

All methods were applied independently to each sample-level AnnData object.

| Method | Study settings | Source |
| --- | --- | --- |
| Leiden | Scanpy implementation; resolution = 0.2 | Included in this folder |
| GraphST | 600 training epochs; clustering with `k = 2` or `k = 3` | [GraphST](https://github.com/JinmiaoChenLab/GraphST) |
| SpaGCN | Clustering with `k = 3` | [SpaGCN](https://github.com/jianhuupenn/SpaGCN) |
| BANKSY | `lambda = 0.2` | [BANKSY](https://prabhakarlab.github.io/Banksy/) |

The original software version or Git commit used for each external method
should be recorded for the final analysis.


