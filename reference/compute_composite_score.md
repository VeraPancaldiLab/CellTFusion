# Compute composite score for cell groups

Computes a composite score by performing Canonical Correlation Analysis
(CCA) between the deconvolution features of a cell group and the score
of its TF module (module eigengene). The score is the weighted sum of
the (scaled) features of the cell group, with the weights of the first
canonical component: the combination of the features that correlates
most with the TF module score (with a single module score, this is the
linear regression of the module score on the features, up to a scaling
factor). As the sign of CCA weights is arbitrary, the score is oriented
so that it is aligned with the TF module: a positive score means that
the sample is on the active side of the module for that cell group. The
score therefore increases with the features of the cell group when they
are positively related to the module, and decreases with them when they
are negatively related.

## Usage

``` r
compute_composite_score(cell_group, module, network, batch = NULL)
```

## Arguments

- cell_group:

  A numeric matrix of cell deconvolution features for a cell group
  (samples x features), with sample names as row names. The module score
  of `network` is matched to these samples by name; every sample must be
  in the network.

- module:

  Character. Name (color) of the TF module the cell group was built
  from; must match a column name of the module matrix in `network`
  exactly.

- network:

  Output of compute.WTCNA() (single- or multi-cohort).

- batch:

  Optional vector with the cohort of each sample, in the row order of
  `cell_group`. It is treated as categorical: the cell group features
  are z-scored within each cohort (mean 0, SD 1 per cohort; a feature
  constant in a cohort is 0 there) before the CCA. The module score is
  used as it is:
  [`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md)
  already scales it within each cohort.

## Value

An unnamed list of two elements:

- `[[1]]`: Numeric matrix (samples x 1) with the composite score, i.e.
  the scaled cell group features projected onto the first canonical
  component, oriented to be aligned with the TF module. With `batch`, it
  has mean 0 within each cohort.

- `[[2]]`: Projection parameters used to score new samples: `xcoef`
  (canonical weights of the first component), `train_means` and
  `train_sds` (column means/SDs used for scaling). The sign of `xcoef`
  follows the orientation of the score. When batch correction was
  applied, only `xcoef` is returned: new samples are then z-scored on
  their own means and SDs, as each training cohort was.

## Details

The CCA is reliable only when the number of features of the cell group
is well below the number of samples. When it reaches the number of
samples, the first canonical correlation is 1 whatever the data, the
weights are not unique and the score is decided by numerical noise. Cell
groups with many features and few samples (e.g. cross-validation folds)
are the most affected.

## Examples

``` r
data("network.tuto")
data("deconv_subgroups.tuto")
cell_groups <- construct_cell_groups(network.tuto, deconv_subgroups.tuto)
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Dendrogram plots saved to: Results/Dendrogram_color_clusters.pdf
#> 44 of 897 features were not assigned to any cluster by the tree cut and are not part of any cell group (pink: 2, yellow: 8, red: 2, blue: 6, brown: 4, turquoise: 12, black: 4, green: 6)
group <- names(cell_groups$Composition)[1]
module <- sub("^Dendrogram_(.*)_group_.*$", "\\1", group) # TF module of the cell group
features <- deconv_subgroups.tuto[[1]][, cell_groups$Composition[[group]], drop = FALSE]
score <- compute_composite_score(features, module, network.tuto)
head(score[[1]])
#>                        [,1]
#> SAM7f0d9cc7f001 -0.14044970
#> SAM4305ab968b90  0.05344901
#> SAMcf018fee2acd -0.10579244
#> SAMcc4675f394a1 -0.26856790
#> SAM49f9b2e57aa5  0.06540482
#> SAM2e7aa8fa0ab3  0.01583567
```
