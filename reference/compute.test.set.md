# Compute composite scores on test set based on previous cell groups

Scores the cell groups learned on a training set in an independent
(test) set, on the same scale as the training scores. Three steps:

1.  the training deconvolution subgroups are rebuilt in the test data
    with
    [`multideconv::replicate_deconvolution_subgroups()`](https://verapancaldilab.github.io/multideconv/reference/replicate_deconvolution_subgroups.html)
    (each subgroup is the median of its member features);

2.  for each cell group, its features are scaled with the training means
    and standard deviations;

3.  they are projected with the training CCA weights of the cell group,
    giving one score per test sample.

Scoring the training samples themselves reproduces their training
scores. Each test sample is scored on its own, so a single sample can be
scored (except for cell groups built with `batch`, see Details).

## Usage

``` r
compute.test.set(dt, cell_groups, group_names = NULL, deconv_test)
```

## Arguments

- dt:

  Deconvolution subgroups of the training set, as returned by
  [`multideconv::compute.deconvolution.analysis()`](https://verapancaldilab.github.io/multideconv/reference/compute.deconvolution.analysis.html)
  (the elements `Deconvolution matrix` and
  `Deconvolution subgroups composition` are used).

- cell_groups:

  The cell groups of the training set, as returned by
  [`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md):
  its elements `Composition` (features of each cell group) and `Weights`
  (CCA projection parameters `xcoef`, `train_means`, `train_sds` of each
  cell group) are used. Older unnamed lists, with these as second and
  third elements, are also accepted.

- group_names:

  Optional character vector with the names of the cell groups to score.
  If `NULL` (default), all cell groups are scored. Names that are not
  cell groups are ignored with a message; an error is raised if none of
  them is a cell group.

- deconv_test:

  A data frame or matrix with the deconvolution features of the test set
  before subgrouping (the same kind of matrix as given to
  [`multideconv::compute.deconvolution.analysis()`](https://verapancaldilab.github.io/multideconv/reference/compute.deconvolution.analysis.html)
  for the training set, with the same column names), samples as rows.

## Value

A data frame with one row per test sample (sample names as row names)
and one column per scored cell group, on the scale of the training
scores. A cell group is scored only if all its features are in the test
set without missing values: a score computed from part of them can be
very different from the training score (even reversed), so the other
cell groups are skipped, with a warning listing them. An error is raised
if no cell group can be scored.

## Details

If the cell groups were built with batch correction (`batch` in
[`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md)),
the features of the test samples are z-scored on their own means and SDs
(as each training cohort was), so `deconv_test` should contain a single
cohort with several samples: the scores then describe each sample
relative to the other test samples. A warning is given when such cell
groups are scored with fewer than 10 test samples.

## Examples

``` r
data("network.tuto")
data("deconv_subgroups.tuto")
data("deconv.tuto")
cell_groups <- construct_cell_groups(network.tuto, deconv_subgroups.tuto)
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> ℹ The deprecated feature was likely used in the dendextend package.
#>   Please report the issue at <https://github.com/talgalili/dendextend/issues>.
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Dendrogram plots saved to: Results/Dendrogram_color_clusters.pdf
#> 44 of 897 features were not assigned to any cluster by the tree cut and are not part of any cell group (pink: 2, yellow: 8, red: 2, blue: 6, brown: 4, turquoise: 12, black: 4, green: 6)
test_scores <- compute.test.set(deconv_subgroups.tuto, cell_groups,
                                deconv_test = deconv.tuto[1:20, ])
test_scores[1:3, 1:3]
#>                 Dendrogram_pink_group_1 Dendrogram_pink_group_2
#> SAM7f0d9cc7f001             -0.14044970             -0.05427430
#> SAM4305ab968b90              0.05344901              0.06916341
#> SAMcf018fee2acd             -0.10579244             -0.06208404
#>                 Dendrogram_pink_group_3
#> SAM7f0d9cc7f001             -0.09753325
#> SAM4305ab968b90              0.05193819
#> SAMcf018fee2acd             -0.17637402
```
