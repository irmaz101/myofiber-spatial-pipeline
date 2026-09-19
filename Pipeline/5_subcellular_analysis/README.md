# Subcellular analysis

Bento was used to explore subcellular transcript localization. This
proof-of-concept analysis is separate from the main fibre and nuclear workflow.

It requires a prepared SpatialData Zarr store containing transcript locations
and matching cell and nuclear segmentation shapes. The study-specific code,
prepared input and software versions are not included in this repository.

For setup, data preparation, RNAforest and RNAcoloc, see the
[Bento repository](https://github.com/ckmah/bento-tools) and
[documentation](https://bento-tools.readthedocs.io/). Use a separate environment
for Bento and its SpatialData dependencies.
