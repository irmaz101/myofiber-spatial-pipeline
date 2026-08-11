# Subcellular analysis

This directory contains the notebook used for the exploratory Bento
proof-of-concept analysis of subcellular transcript localization.

## Included file

| File | Description |
| --- | --- |
| `run_bento.ipynb` | Runs the study-specific Bento, RNAforest and RNAcoloc workflow on a prepared SpatialData object. |

## Installation

Bento is an external workflow and is not redistributed or pinned in this
repository. Install Bento following the instructions in its
[original repository](https://github.com/ckmah/bento-tools). Because Bento and
its SpatialData dependencies evolve rapidly, a separate environment is
recommended.

## Input

The notebook expects a prepared SpatialData Zarr store containing transcript
locations and matching cell and nuclear segmentation shapes. Input paths and
analysis settings must be updated inside the notebook before execution.

## Scope

This notebook documents the study-specific exploratory analysis and is not
required for the main segmentation, preprocessing, clustering, annotation or
SVG-detection workflow.
