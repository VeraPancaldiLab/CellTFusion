# Project test-set cell group scores using training CCA parameters

Scales a test-set cell group matrix using the mean and standard
deviation stored from training, then projects it onto the first
canonical component learned during training.

## Usage

``` r
# S3 method for class 'test.score'
compute(cell_group, projection_params)
```

## Arguments

- cell_group:

  A numeric matrix of cell deconvolution features for test samples
  (samples x features). Column names must overlap with training
  features.

- projection_params:

  A list containing:

  xcoef

  :   Named numeric matrix of CCA canonical weights (features x 1).

  train_means

  :   Named numeric vector of training column means.

  train_sds

  :   Named numeric vector of training column standard deviations.

## Value

A numeric matrix (samples x 1) of projected canonical scores for the
test cohort, or `NULL` if no features overlap.
