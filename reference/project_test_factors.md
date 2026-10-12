# Project test-set samples onto training NMF factors

Convenience wrapper that combines
[`compute.test.set()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.test.set.md)
and
[`project_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/project_factors.md)
into a single call. Given a trained
[`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)
result and a test deconvolution matrix, it (1) reconstructs cell group
composite scores for the test samples using the training CCA loadings,
and (2) projects those scores onto the trained NMF latent space.

## Usage

``` r
project_test_factors(train_processed, test_deconv)
```

## Arguments

- train_processed:

  A list returned by
  [`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md),
  containing at minimum `Processed_deconvolution`, `Cell_groups`, and
  `Latent_spaces`.

- test_deconv:

  A numeric matrix or data frame of the deconvolution features of the
  test samples before subgrouping (samples x features; the same kind of
  matrix as given to
  [`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)
  as `deconv`), with the column names used during training. If the model
  was trained with batch correction, the test samples are centred on
  their own means (as each training cohort was), so they should come
  from a single cohort.

## Value

A numeric matrix (test samples x NMF factors) of non-negative projected
factor scores. Cell groups with features missing in `test_deconv` are
skipped by
[`compute.test.set()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.test.set.md)
(with a warning), and the projection uses the other cell groups.

## Examples

``` r
if (FALSE) { # \dontrun{
res <- CellTFusion(raw.counts_train, deconv = deconv_train)
test_factors <- project_test_factors(res, deconv_test)
} # }
```
