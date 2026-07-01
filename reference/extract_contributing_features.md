# Extract top-contributing features from NMF latent factors

For each NMF factor (column of the W basis matrix), selects the features
whose weight exceeds the specified quantile threshold.

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

  Numeric between 0 and 1. Features with weight above this quantile of
  the factor's positive weights are retained. Default 0.7.

## Value

A named list (one element per factor) of named numeric vectors, where
names are feature names and values are their NMF basis weights, sorted
descending.
