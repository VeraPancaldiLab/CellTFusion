# Extract cells from cell type groups

Extract cells from cell type groups

## Usage

``` r
extract_cells(groups, cells_extra = NULL)
```

## Arguments

- groups:

  A character vector of cell type group names, typically from
  cell.groups.computation()

- cells_extra:

  Optional character vector of additional cell type names to include in
  the extraction

## Value

A character vector of unique cell types found in the groups, ignoring
method and signature suffixes
