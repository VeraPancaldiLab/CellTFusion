# Remove highly correlated cell groups

Computes pairwise Spearman correlations among cell group score vectors
and removes one member of each pair whose absolute correlation exceeds
`threshold`.

## Usage

``` r
remove.cell.groups.corr(data, threshold = 0.95)
```

## Arguments

- data:

  A list of three elements:

  scores

  :   A numeric data frame or matrix of cell group scores (samples x
      groups).

  compositions

  :   A named list of cell-type vectors describing group membership.

  loadings

  :   A named list of loading vectors corresponding to each cell group.

- threshold:

  Numeric. Correlation threshold above which one of a correlated pair is
  removed. Default is 0.95.

## Value

A list of three elements (scores, compositions, loadings) with redundant
cell groups removed.
