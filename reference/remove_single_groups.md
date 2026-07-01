# Remove cell groups composed of a single cell type

Filters out cell groups whose composition contains only one cell type,
as these groups lack multi-cellular context.

## Usage

``` r
remove_single_groups(cell.values, cell.composition, cell.loadings)
```

## Arguments

- cell.values:

  A list of numeric vectors of cell group scores.

- cell.composition:

  A list of character vectors describing cell-type membership per group.

- cell.loadings:

  A list of loading vectors corresponding to each cell group.

## Value

A list of three elements (scores, compositions, loadings) with singleton
groups removed, or `NULL` if all groups are removed.
