# Perform pairwise Spearman correlation across all features

Perform pairwise Spearman correlation across all features

## Usage

``` r
correlation(data)
```

## Arguments

- data:

  A numeric matrix or data frame where columns are features.

## Value

A data frame of pairwise significant correlations (p \< 0.05), with
columns `measure1`, `measure2`, `r`, `p`, `sig_p`, `p_if_sig`,
`r_if_sig`, and `AbsR`, ordered by descending `r`.
