# Compute composite scores on test set based on previous cell groups

This function simulates cell subgroups compositions from deconvolution
results in the test deconvolution data, and calculates composite scores
using provided cell groups and features in order to replicate cell
groups coming from a training set into an independent set. The composite
scores represent summarized information from cell subgroup profiles.

## Usage

``` r
compute.test.set(deconv_res, cell_groups, features, deconvolution_test)
```

## Arguments

- deconv_res:

  Deconvolution subgroups of the training set, as returned by
  [`multideconv::compute.deconvolution.analysis()`](https://verapancaldilab.github.io/multideconv/reference/compute.deconvolution.analysis.html)
  (the elements `Deconvolution matrix` and
  `Deconvolution subgroups composition` are used).

- cell_groups:

  A list with three elements, as returned by
  [`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md)
  on the training set:

  - cell groups scores

  - composition: A named list of character vectors where each element
    represents cells belonging to a specific group.

  - loadings: A corresponding list of CCA projection parameters
    (`xcoef`, `train_means`, `train_sds`) for each cell group.

- features:

  A character vector of feature names to select relevant cell groups. An
  error is raised if none of them match a cell group name.

- deconvolution_test:

  A data frame or matrix of deconvolution features for the test set,
  with cells as columns and samples as rows.

## Value

A data frame where each column corresponds to a composite score
calculated for each feature group in the test set. If no composite
scores can be computed due to zero variance, returns an empty data frame
with a printed message.

## Details

The function first replicates the training cell subgroups in the test
deconvolution data with
[`multideconv::replicate_deconvolution_subgroups()`](https://verapancaldilab.github.io/multideconv/reference/replicate_deconvolution_subgroups.html)
(each subgroup is the median of its member features; subgroups or
features missing in the test set are `NA` and are skipped). Then it
extracts the relevant cells for each feature and calculates composite
scores. If the cell groups were built with batch correction, the test
samples are centred on their own means (as each training cohort was), so
`deconvolution_test` should contain a single cohort.
