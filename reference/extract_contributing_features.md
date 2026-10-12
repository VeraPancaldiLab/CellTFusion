# Extract top-contributing features from NMF latent factors

For each NMF factor (column of the W basis matrix), selects the features
(cell group parts `<group>_pos` and `<group>_neg`) whose weight is at or
above the specified quantile of the factor's positive weights (zero
weights are left out before computing the quantile).

## Usage

``` r
extract_contributing_features(latent_factors, quantile_cutoff = 0.7)
```

## Arguments

- latent_factors:

  A list returned by
  [`compute.latent_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.latent_factors.md),
  containing at minimum a matrix `W` (features x factors) of NMF basis
  weights.

- quantile_cutoff:

  Numeric between 0 and 1. Features with weight at or above this
  quantile of the factor's positive weights are retained. Default 0.7.

## Value

A named list (one element per factor) of named numeric vectors, where
names are feature names (with their `_pos` / `_neg` suffix) and values
are their NMF basis weights, sorted descending.
