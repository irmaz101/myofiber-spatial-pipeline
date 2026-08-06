# Spatially variable gene detection

This folder contains the wrappers used to run eight spatially variable gene
(SVG) methods on sample-level myofibre AnnData objects. The statistical methods
are provided by the original packages; these scripts only prepare the input,
set the study parameters and save the selected genes.

## Input

Each script processes one `.h5ad` file. The AnnData object must contain:

- raw counts in `adata.layers["counts"]`;
- normalized, unscaled expression in `adata.raw`;
- `x_centroid` and `y_centroid` in `adata.obs`.

The first two fields are created by `Pipeline/1_preprocessing/preprocessing.py`.

## Methods and parameters

| Method | Expression | Parameters | SVG selection |
| --- | --- | --- | --- |
| Moran's I | `adata.X` | Squidpy; 6 spatial neighbours | BH-FDR < 0.05 |
| SpatialDE | Raw counts | NaiveDE stabilization; regression of `log(total_counts)` | q-value < 0.05 |
| Hotspot | Raw counts | Bernoulli model; 20 unweighted neighbours | FDR < 0.05 |
| scGCO | Raw counts | Cell Ranger normalization; genes detected in at least 2 objects | FDR < 0.05 |
| SOMDE | Raw counts | SOM grid size = 10 | q-value < 0.05 |
| SMASH | `adata.X` | All covariance kernels; `mean_only=False`; `forcePD=False` | BY-adjusted P < 0.01 |
| SpaGFT | AnnData expression | `ratio_neighbors=1`; `filter_peaks=True`; `S=6` | SpaGFT score cutoff and FDR < 0.05 |
| SVGbit | Normalized expression from `adata.raw` | Variance filter = 0; quantile filter = 0.99; 6 neighbours | Top 1,000 genes by AI score |

The required `--output` CSV contains only genes meeting the stated selection
rule. This is the file used by the benchmarking scripts, which treat each row
as one detected SVG.

An optional `--all-output` CSV can also be saved. It contains every gene
returned by the method and an `is_svg` column indicating whether the gene met
the selection rule.

## Usage

All scripts use the same input and output arguments. For example:

```bash
python Pipeline/4_SVG/morans.py \
    --input anndata_mf_density/sample_01.h5ad \
    --output results/moran_i/sample_01.csv \
    --all-output results_all/moran_i/sample_01.csv
```

Omit `--all-output` if only the selected SVG table is needed.

Equivalent commands can be used for the other methods:

```bash
python Pipeline/4_SVG/spatial_de.py -i INPUT.h5ad -o OUTPUT.csv
python Pipeline/4_SVG/hotspot.py    -i INPUT.h5ad -o OUTPUT.csv
python Pipeline/4_SVG/scgco.py      -i INPUT.h5ad -o OUTPUT.csv
python Pipeline/4_SVG/somde.py      -i INPUT.h5ad -o OUTPUT.csv
python Pipeline/4_SVG/smash.py      -i INPUT.h5ad -o OUTPUT.csv
python Pipeline/4_SVG/spagft.py     -i INPUT.h5ad -o OUTPUT.csv
python Pipeline/4_SVG/svgbit.py     -i INPUT.h5ad -o OUTPUT.csv
```

Use one output directory per method, with identical sample basenames:

```text
results/
├── hotspot/sample_01.csv
├── moran_i/sample_01.csv
├── scgco/sample_01.csv
├── smash/sample_01.csv
├── somde/sample_01.csv
├── spagft/sample_01.csv
├── spatialde/sample_01.csv
└── svgbit/sample_01.csv
```

The same layout can be used under `results_all/` for the optional complete
gene-level tables. Use files under `results/`, not `results_all/`, as input to
the SVG benchmarking scripts.

scGCO and SMASH use local source directories by default. Alternative locations
can be supplied with `--source`.

## Software

- [Squidpy](https://github.com/scverse/squidpy)
- [SpatialDE](https://github.com/teichlab/spatialde)
- [Hotspot](https://github.com/YosefLab/Hotspot)
- [scGCO](https://github.com/WangPeng-Lab/scGCO)
- [SOMDE](https://github.com/WhirlFirst/somde)
- [SMASH](https://github.com/sealx017/SMASH-)
- [SpaGFT](https://github.com/jxLiu-bio/SpaGFT)
- [SVGbit](https://github.com/CPenglab/svgbit)

Some packages have incompatible dependency requirements. Separate environments
may therefore be needed. Record the version or Git commit used for each final
analysis and cite the corresponding method publication.
