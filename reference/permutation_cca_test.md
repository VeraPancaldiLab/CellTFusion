# Permutation test for the first canonical correlation

Estimates the null distribution of the first canonical correlation by
row-permuting X and recomputing CCA `n_perm` times, then returns the
empirical p-value (proportion of null correlations \>= observed).

## Usage

``` r
permutation_cca_test(X, Y, n_perm = 999)
```

## Arguments

- X:

  Numeric matrix (samples x features), already scaled.

- Y:

  Numeric matrix (samples x features), already scaled.

- n_perm:

  Integer. Number of permutations. Default 999.

## Value

A list with `r` (observed first canonical correlation) and `p_value`
(empirical one-sided p-value).
