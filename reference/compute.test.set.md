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
compute.test.set(deconv_res, cell_groups, features = NULL, deconvolution_test)
```

## Arguments

- deconv_res:

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

- features:

  Optional character vector with the names of the cell groups to score.
  If `NULL` (default), all cell groups are scored. Names that are not
  cell groups are ignored with a message; an error is raised if none of
  them is a cell group.

- deconvolution_test:

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
the test samples are centred on their own means (as each training cohort
was), so `deconvolution_test` should contain a single cohort with
several samples: the scores then describe each sample relative to the
other test samples. A warning is given when such cell groups are scored
with fewer than 10 test samples.

## Examples

``` r
data("network.tuto")
data("deconv_subgroups.tuto")
data("deconv.tuto")
cell_groups <- construct_cell_groups(network.tuto, deconv_subgroups.tuto)
#> 20 of 764 features were not assigned to any cluster by the tree cut and are not part of any cell group (blue: 6, pink: 2, red: 2, yellow: 8, black: 2)
test_scores <- compute.test.set(deconv_subgroups.tuto, cell_groups, deconvolution_test = deconv.tuto[1:20, ])
test_scores[1:3, 1:3]
#>                 Dendrogram_brown_group_1 Dendrogram_brown_group_2
#> SAM7f0d9cc7f001               0.03998411               0.06757056
#> SAM4305ab968b90              -0.13368092              -0.14208692
#> SAMcf018fee2acd              -0.02670350               0.03518956
#>                 Dendrogram_brown_group_3
#> SAM7f0d9cc7f001               0.06172535
#> SAM4305ab968b90              -0.15458310
#> SAMcf018fee2acd              -0.03169306
```
