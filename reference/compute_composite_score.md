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

  Character. Name (color) of the TF module the cell group was built
  from; must match a column name of the module matrix in
  `tfs.module.network` exactly.

- tfs.module.network:

  Output of compute.WTCNA().

- batch:

  Optional vector indicating batch assignment for samples. It is treated
  as categorical: per-batch means are regressed out of the cell group
  features, the module eigengene and the module TFs before the CCA.

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

An unnamed list of two elements:

- `[[1]]`: Numeric matrix (samples x 1) with the composite score, i.e.
  the scaled cell group features projected onto the first canonical
  component.

- `[[2]]`: Projection parameters used to score new samples: `xcoef`
  (canonical weights of the first component), `train_means` and
  `train_sds` (column means/SDs used for scaling). `train_means` is
  `NULL` when batch correction was applied: new samples are then centred
  on their own means, as each training cohort was.

If the permutation test is not significant (and `discard = TRUE`),
returns `list("NA", "NA")`.
