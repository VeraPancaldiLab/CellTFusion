# Extract top-contributing features from NMF latent factors

For each NMF factor (column of the W basis matrix), selects the `top_n`
cell groups with the highest weight (cell groups with zero weight are
never selected).

## Usage

``` r
extract_contributing_features(latent_factors, top_n = 5)
```

## Arguments

- latent_factors:

  A list returned by
  [`compute.latent_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.latent_factors.md),
  containing at minimum a matrix `W` (features x factors) of NMF basis
  weights.

- top_n:

  Integer. Number of cell groups selected per factor. Default 5.

## Value

A named list (one element per factor) of named numeric vectors, where
names are cell group names and values are their NMF basis weights,
sorted descending.
