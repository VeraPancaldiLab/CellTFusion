# Construct cell groups based on TF networks and deconvolution

Identifies and projects cell groups using module relationships derived
from TF networks and deconvolution outputs. If a binary trait is
specified, the function splits the data and constructs cell groups for
both classes (supervised analysis).

## Usage

``` r
construct_cell_groups(
  network,
  dt,
  batch = NULL,
  pval = 0.05,
  high_corr_groups = 0.8,
  clustering.method = "ward.D2"
)
```

## Arguments

- network:

  A list containing TF networks for cell types.

- dt:

  A list containing deconvolution subgroup structures.

- batch:

  Optional vector indicating batch assignment for samples.

- pval:

  P-value threshold used to filter module relationships. Default: 0.05.

- high_corr_groups:

  Correlation threshold to merge or remove redundant cell groups.
  Default: 0.9.

- clustering.method:

  Clustering method for hierarchical clustering. Default: "ward.D2".

## Value

A list of 3 elements:

- scores:

  A data frame or matrix with the projected cell group scores (samples x
  groups).

- composition:

  A named list where each element is a character vector of original cell
  types per group.

- loadings:

  A list of numeric vectors indicating the loadings (feature
  contributions) for each group.
