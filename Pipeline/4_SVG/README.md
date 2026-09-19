# Spatially variable gene detection

Eight SVG methods are run separately on each preprocessed myofibre AnnData
file. Each wrapper prepares the method input and saves its gene-level results.

## Inputs and environments

Use the files from [preprocessing](../1_preprocessing/README.md). Required
fields are:

| Field | Used by |
| --- | --- |
| `adata.raw` with raw counts | Moran's I, Hotspot, SMASH, SpaGFT, SVGbit |
| `adata.layers["counts"]` | SpatialDE, scGCO, SOMDE |
| `obs["x_centroid"]`, `obs["y_centroid"]` | All methods |

Create a separate environment per method, using the Python version noted in
its requirements file. For example:

```bash
conda create -n myofiber-morans python=3.10
conda activate myofiber-morans
python -m pip install -r Pipeline/4_SVG/requirements/morans.txt
```

Requirements filenames follow the method names (`hotspot.txt`, `somde.txt`,
etc.), even where the wrapper has a `run_` prefix. SVGbit dependencies are
unpinned because the exact study versions were not recorded.

scGCO and SMASH also need local source directories, supplied with `--source`.
The scGCO wrapper expects `scGCO_simple.py` (default directory: `scgco_utils/`);
that helper is not included here. SMASH defaults to `SMASH/SMASH/`.

## Method settings

| Method | Wrapper | Expression and parameters | SVG selection |
| --- | --- | --- | --- |
| [Moran's I / Squidpy](https://github.com/scverse/squidpy) | `morans.py` | Counts normalized to 10,000 and log1p-transformed; 6 spatial neighbours | BH-FDR < 0.05 |
| [SpatialDE](https://github.com/teichlab/spatialde) | `spatial_de.py` | Raw counts; NaiveDE stabilization and regression of `log(total_counts)` | q-value < 0.05 |
| [Hotspot](https://github.com/YosefLab/Hotspot) | `run_hotspot.py` | Raw counts; Bernoulli model; 20 unweighted neighbours | FDR < 0.05 |
| [scGCO](https://github.com/WangPeng-Lab/scGCO) | `scgco.py` | Raw counts; Cell Ranger normalization; detection in at least 2 objects | FDR < 0.05 |
| [SOMDE](https://github.com/WhirlFirst/somde) | `run_somde.py` | Raw counts; `k = 10` | q-value < 0.05 |
| [SMASH](https://github.com/sealx017/SMASH-) | `smash.py` | Counts normalized to 10,000 and log1p-transformed; all covariance kernels; `mean_only=False`; `forcePD=False` | BY-adjusted P ≤ 0.05 |
| [SpaGFT](https://github.com/jxLiu-bio/SpaGFT) | `spagft.py` | Counts normalized to 10,000 and log1p-transformed; `ratio_neighbors=1`; `filter_peaks=True`; `S=6` | Score cutoff and FDR ≤ 0.05 |
| [SVGbit](https://github.com/CPenglab/svgbit) | `run_svgbit.py` | Raw counts with log-CPM normalization; variance filter = 0; no upper-expression quantile filter; 6 neighbours | Top 1,000 genes by AI score |

## Run a method

Every wrapper accepts `--input`, `--output` and optional `--all-output`.
Use `--help` for method-specific arguments. From the repository root:

```bash
python Pipeline/4_SVG/morans.py \
    --input anndata_mf_density/sample_01.h5ad \
    --output results/moran_i/sample_01.csv \
    --all-output results_all/moran_i/sample_01.csv
```

Substitute the wrapper and output directory for another method, for example:

```bash
python Pipeline/4_SVG/run_hotspot.py \
    --input anndata_mf_density/sample_01.h5ad \
    --output results/hotspot/sample_01.csv \
    --all-output results_all/hotspot/sample_01.csv
```

## Result files

| Argument | Contents | Benchmark use |
| --- | --- | --- |
| `--output` | Genes meeting the method's selection rule | SVG counts, method overlap and technical-replicate overlap |
| `--all-output` | Every returned gene, its statistics and an `is_svg` flag | Synthetic Kendall and AUPRC benchmarks |

Use identical sample basenames across methods. The benchmark scripts expect
method directories named `hotspot`, `moran_i`, `scgco`, `smash`, `somde`,
`spagft`, `spatialde` and `svgbit`. Keep selected-gene results and complete
results in separate parent directories, as in the examples above.
