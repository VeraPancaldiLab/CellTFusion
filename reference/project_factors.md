# Project cell group scores onto trained NMF latent factors

Applies the positive-negative split to new samples and projects them
onto the latent space learned during training using the Moore-Penrose
pseudoinverse of the scaled basis matrix W.

## Usage

``` r
project_factors(latent_spaces, scores_test)
```

## Arguments

- latent_spaces:

  A list returned by
  [`compute.latent_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.latent_factors.md),
  containing at minimum `W` (features x rank), `rank`, and the RcppML
  model object with `nmf_model$d` scaling vector.

- scores_test:

  A samples x cell groups matrix of signed CCA composite scores for the
  test cohort. Column names must match those used in training (before
  the \_pos/\_neg suffix was added).

## Value

A non-negative matrix of size samples x rank containing the projected
NMF factor scores for the test cohort.

## Details

The RcppML NMF decomposition is A = W %*% diag(d) %*% H, not A = W %\*%
H. The diagonal scaling vector d must be absorbed into W before
computing the pseudoinverse, otherwise the projected scores are on the
wrong scale.
