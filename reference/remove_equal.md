# Remove cell groups with duplicate composition

Identifies cell groups whose cell-type composition is identical to
another group and removes the duplicates, keeping only the first
occurrence.

## Usage

``` r
remove_equal(cell.values, cell.composition, cell.loadings)
```

## Arguments

- cell.values:

  A list of numeric vectors of cell group scores.

- cell.composition:

  A list of character vectors describing cell-type membership per group.

- cell.loadings:

  A list of loading vectors corresponding to each cell group.

## Value

A list of three elements:

- `[[1]]`: Deduplicated cell group scores.

- `[[2]]`: Deduplicated cell group compositions.

- `[[3]]`: Deduplicated cell group loadings.
