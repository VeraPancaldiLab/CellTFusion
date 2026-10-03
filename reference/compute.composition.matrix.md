# Compute a cell-type composition matrix from deconvolution subgroups

Builds a binary presence matrix indicating which original cell types are
present in higher-level cell groups.

## Usage

``` r
compute.composition.matrix(
  deconvolution.subgroupped,
  cell.groups,
  cells_extra = NULL
)
```

## Arguments

- deconvolution.subgroupped:

  Deconvolution subgroups as returned by
  [`multideconv::compute.deconvolution.analysis()`](https://verapancaldilab.github.io/multideconv/reference/compute.deconvolution.analysis.html)
  (the elements `Deconvolution matrix` and
  `Deconvolution subgroups composition` are used).

- cell.groups:

  A list containing cell group definitions, where the second element
  holds the groupings (e.g. the output of
  [`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md)).

- cells_extra:

  Optional vector of additional cell identifiers to consider during
  extraction.

## Value

A binary matrix (data.frame) where rows are cell groups and columns are
cell types (1 = present, 0 = absent). Each cell group appears twice, as
`<group>_pos` and `<group>_neg` (same composition), to match the feature
names used by
[`compute.latent_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.latent_factors.md).
Cell types not present in any cell group are dropped (with a warning).
