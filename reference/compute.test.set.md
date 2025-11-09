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

  A list containing deconvolution results, including subgroup
  compositions (list of data frames or matrices).

- cell_groups:

  A list with three elements:

  - cell groups scores

  - composition: A named list of character vectors where each element
    represents cells belonging to a specific group.

  - loadings: A corresponding list of numeric vectors (loadings) for
    each cell group.

- features:

  A character vector of feature names to select relevant cell groups.

- deconvolution_test:

  A data frame or matrix of deconvolution features for the test set,
  with cells as columns and samples as rows.

## Value

A data frame where each column corresponds to a composite score
calculated for each feature group in the test set. If no composite
scores can be computed due to zero variance, returns an empty data frame
with a printed message.

## Details

The function first simulates cell subgroups by computing median values
across specified iterations and joins them with the original test
deconvolution data. Then it extracts the relevant cells for each feature
and calculates composite scores.
