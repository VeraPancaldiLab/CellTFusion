# Compute latent factors from cell group scores using NMF

Decomposes signed CCA composite scores into non-negative latent factors
by splitting each score into its positive and negative direction
components, then applying Non-negative Matrix Factorization via RcppML
(fast C++ backend).

## Usage

``` r
# S3 method for class 'latent_factors'
compute(X, rank = NULL, seed = 123, file_name = NULL, return = TRUE)
```

## Arguments

- X:

  Numeric matrix of size samples x cell groups (signed CCA composite
  scores).

- rank:

  Integer; number of NMF factors. If NULL, estimated automatically via
  elbow on reconstruction MSE across ranks 2:8.

- seed:

  Random seed. Default 123.

- file_name:

  Optional character. If provided, saves results to this file path.

- return:

  Logical. If TRUE (default), returns the result as an R object.

## Value

A named list with:

- Z:

  Sample-level NMF factor scores (samples x rank). Non-negative.

- W:

  Feature weights per factor (2 x n_CGs x rank). Non-negative.

- nmf_input:

  The positive-negative split matrix fed to NMF (samples x 2 x n_CGs).

- rank:

  The rank used.

## Details

Signed CCA scores are decomposed as: score_pos = max(score, 0) – patient
aligned with TF program score_neg = max(-score, 0) – patient
anti-aligned with TF program Both are concatenated column-wise before
NMF. Column names are suffixed with "\_pos" and "\_neg" to track
direction.
