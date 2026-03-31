# Compute composite score for cell groups

Computes a composite score by performing Canonical Correlation Analysis
(CCA) between cell group features and corresponding TF module scores.

## Usage

``` r
compute_composite_score(
  cell_group,
  module_group,
  tfs.module.network,
  batch = NULL,
  discard = T
)
```

## Arguments

- cell_group:

  A numeric matrix of cell deconvolution features for a cell group
  (samples x features).

- module_group:

  A character vector indicating TF module group colors corresponding to
  the cell group (can be obtained via
  [`extract_colors()`](https://verapancaldilab.github.io/CellTFusion/reference/extract_colors.md)).

- tfs.module.network:

  Output of compute.WTCNA().

- batch:

  Optional vector indicating batch assignment for samples.

- discard:

  Logical; whether to discard cell groups whose canonical correlation is
  below 0.6 (default TRUE).

## Value

A list with:

- `selected_components`: Numeric matrix of the first canonical component
  scores across samples.

- `xcoef`: The canonical weights (coefficients) for the cell group
  features.

If discarded due to low correlation, returns `list("NA", "NA")`.
