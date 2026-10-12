# Project test-set cell group scores using training CCA parameters

Scales a test-set cell group matrix using the mean and standard
deviation stored from training, then projects it onto the first
canonical component learned during training. All the features of the
cell group are needed: a score computed from part of them (the others
set to their training mean) can be very different from the training
score, even reversed, because the CCA weights of correlated features can
have opposite signs.

## Usage

``` r
compute.test.score(cell_group, projection_params)
```

## Arguments

- cell_group:

  A numeric matrix of cell deconvolution features for test samples
  (samples x features). It must contain every feature of the cell group
  (the row names of `xcoef`); other columns are ignored.

- projection_params:

  A list containing:

  xcoef

  :   Named numeric matrix of CCA canonical weights (features x 1).

  train_means

  :   Named numeric vector of training column means, or `NULL` if the
      training used batch correction, in which case the test data are
      centred on their own column means.

  train_sds

  :   Named numeric vector of training column standard deviations.

## Value

A numeric matrix (samples x 1) of projected canonical scores for the
test cohort, or `NULL` (with a warning) if some features of the cell
group are missing.
