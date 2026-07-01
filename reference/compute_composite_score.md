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
  discard = T,
  pval = 0.05,
  n_perm = 999
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

  Logical; whether to discard cell groups that do not pass the
  permutation test for the first canonical correlation (default TRUE).

- pval:

  Numeric. Significance threshold for the permutation test (default
  0.05).

- n_perm:

  Integer. Number of permutations used to build the null distribution
  (default 999).

## Value

A list with:

- `selected_components`: Numeric matrix of the first canonical component
  scores across samples.

- `xcoef`: The canonical weights (coefficients) for the cell group
  features.

If discarded due to low correlation, returns `list("NA", "NA")`.
