# Subcellular analysis

This directory documents the exploratory Bento proof-of-concept analysis of
subcellular transcript localization. Study-specific Bento code is not
redistributed because the analysis depends on a separately prepared
SpatialData object and version-specific Bento interfaces.

## Installation

Bento is an external workflow and is not redistributed or pinned in this
repository. Install Bento following the instructions in its
[original repository](https://github.com/ckmah/bento-tools) and reproduce the
workflow from the project's tutorials. Because Bento and its SpatialData
dependencies evolve rapidly, a separate environment is recommended.

## Input

The workflow requires a prepared SpatialData Zarr store containing transcript
locations and matching cell and nuclear segmentation shapes. Follow the
[official Bento tutorials](https://bento-tools.readthedocs.io/) for current
data preparation, RNAforest and RNAcoloc usage.

## Scope

This exploratory analysis is not required for the main segmentation,
preprocessing, clustering, annotation or SVG-detection workflow.
