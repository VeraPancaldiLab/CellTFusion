# Remove highly correlated cell groups

Remove highly correlated cell groups

## Usage

``` r
remove.cell.groups.corr(data, threshold = 0.95)
```

## Arguments

- data:

  A named list with three elements:

  scores

  :   A data frame or matrix of cell group scores.

  compositions

  :   A named list of cell group compositions.

  loadings

  :   A named list of loadings corresponding to the cell groups.

- threshold:

  Numeric value for correlation threshold above which features are
  considered highly correlated (default 0.95)

## Value

A list containing:

- Cell group scores after removal/combination of highly correlated
  features

- Cell group compositions updated accordingly

- Loadings updated accordingly
