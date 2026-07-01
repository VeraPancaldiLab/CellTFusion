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

  A numeric matrix or data frame of deconvolution features for the test
  samples (samples x cell types). Column names must match those used
  during training.

## Value

A numeric matrix (test samples x NMF factors) of non-negative projected
factor scores.
