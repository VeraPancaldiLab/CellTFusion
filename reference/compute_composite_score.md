# Compute composite score for cell groups

Computes a composite score by performing Canonical Correlation Analysis
(CCA) between the deconvolution features of a cell group and the TFs of
its module. The score is the weighted sum of the (scaled) features of
the cell group, with the weights of the first canonical component. As
the sign of CCA weights is arbitrary, the score is oriented so that it
is aligned with the TF module: a positive score means that the sample is
on the active side of the module for that cell group. The score
therefore increases with the features of the cell group when they are
positively related to the module, and decreases with them when they are
negatively related. The relation of the cell group to its module (mean
correlation of its features with the module score) is returned to
translate the score into cell abundance.

## Usage

``` r
compute_composite_score(
  cell_group,
  module_group,
  tfs.module.network,
  batch = NULL
)
```

## Arguments

- cell_group:

  A numeric matrix of cell deconvolution features for a cell group
  (samples x features), with sample names as row names. The module score
  and TF activities of `tfs.module.network` are matched to these samples
  by name; every sample must be in the network.

- module_group:

  Character. Name (color) of the TF module the cell group was built
  from; must match a column name of the module matrix in
  `tfs.module.network` exactly.

- tfs.module.network:

  Output of compute.WTCNA() (single- or multi-cohort).

- batch:

  Optional vector with the cohort of each sample, in the row order of
  `cell_group`. It is treated as categorical: per-batch means are
  regressed out of the cell group features, the module score and the
  module TFs before the CCA.

## Value

An unnamed list of three elements:

- `[[1]]`: Numeric matrix (samples x 1) with the composite score, i.e.
  the scaled cell group features projected onto the first canonical
  component, oriented to be aligned with the TF module.

- `[[2]]`: Projection parameters used to score new samples: `xcoef`
  (canonical weights of the first component), `train_means` and
  `train_sds` (column means/SDs used for scaling). The sign of `xcoef`
  follows the orientation of the score. `train_means` is `NULL` when
  batch correction was applied: new samples are then centred on their
  own means, as each training cohort was.

- `[[3]]`: Relation of the cell group to its TF module: the mean
  correlation of its features with the module score (positive or
  negative).

## Details

The CCA is reliable only when the number of variables (features of the
cell group plus TFs of the module) is well below the number of samples.
When it reaches the number of samples, the first canonical correlation
is 1 whatever the data, the weights are not unique and the score is
decided by numerical noise (it can change completely with the sample
order or a tiny change in the data). Large modules with few samples
(e.g. cross-validation folds) are the most affected.

## Examples

``` r
data("network.tuto")
data("deconv_subgroups.tuto")
cell_groups <- construct_cell_groups(network.tuto, deconv_subgroups.tuto)
#> 20 of 764 features were not assigned to any cluster by the tree cut and are not part of any cell group (blue: 6, pink: 2, red: 2, yellow: 8, black: 2)
group <- cell_groups$Module_relation[1, ]
features <- deconv_subgroups.tuto[[1]][, cell_groups$Composition[[group$cell_group]], drop = FALSE]
score <- compute_composite_score(features, group$module, network.tuto)
head(score[[1]])
#>                         [,1]
#> SAM7f0d9cc7f001  0.039984111
#> SAM4305ab968b90 -0.133680921
#> SAMcf018fee2acd -0.026703505
#> SAMcc4675f394a1  0.061549350
#> SAM49f9b2e57aa5 -0.082696177
#> SAM2e7aa8fa0ab3 -0.004501591
```
