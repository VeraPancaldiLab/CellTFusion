# Remove cell groups composed of a single feature

Filters out cell groups whose composition contains only one
deconvolution feature, as these groups lack multi-cellular context.

## Usage

``` r
remove_single_groups(cell.values, cell.composition, cell.loadings)
```

## Arguments

- cell.values:

  A list of numeric vectors of cell group scores.

- cell.composition:

  A list of character vectors with the deconvolution features of each
  group.

- cell.loadings:

  A list of loading vectors corresponding to each cell group.

## Value

A list of three elements (scores, compositions, loadings) with singleton
groups removed, or `NULL` if all groups are removed.
