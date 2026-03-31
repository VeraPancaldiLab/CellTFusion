# Compute a cell-type composition matrix from deconvolution subgroups

Builds a binary presence matrix indicating which original cell types are
present in higher-level cell groups.

## Usage

``` r
# S3 method for class 'composition.matrix'
compute(deconvolution.subgroupped, cell.groups, cells_extra = NULL)
```

## Arguments

- deconvolution.subgroupped:

  A list containing "Deconvolution subgroups composition" per model.

- cell.groups:

  A list containing cell group definitions, where the second element
  holds the groupings.

- cells_extra:

  Optional vector of additional cell identifiers to consider during
  extraction.

## Value

A binary matrix (data.frame) where rows are cell groups and columns are
cell types (1 = present, 0 = absent).
