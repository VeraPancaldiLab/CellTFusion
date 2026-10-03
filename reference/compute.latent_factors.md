# Compute latent factors from cell group scores using NMF

Decomposes signed CCA composite scores into non-negative latent factors
by splitting each score into its positive and negative direction
components, then applying Non-negative Matrix Factorization via RcppML
(fast C++ backend).

## Usage

``` r
compute.latent_factors(
  X,
  rank = NULL,
  seed = 123,
  file_name = NULL,
  return = TRUE
)
```

## Arguments

- X:

  Numeric matrix of size samples x cell groups (signed CCA composite
  scores).

- rank:

  Integer; number of NMF factors. If NULL, estimated automatically at
  the elbow of the reconstruction MSE across ranks 2:8 (the rank with
  the largest second difference of the MSE curve).

- seed:

  Random seed used for the NMF fits. Default 123. The caller's random
  number generator state is restored when the function returns.

- file_name:

  Optional character suffix for the saved patient-mixture plot.

- return:

  Logical. If TRUE (default), saves the patient-mixture barplot to
  `Results/NMF_patient_mixture_<file_name>.pdf`
  (`Results/NMF_patient_mixture.pdf` if `file_name` is `NULL`).

## Value

A named list with:

- Z:

  Sample-level NMF factor scores (samples x rank). Non-negative.

- W:

  Feature weights per factor ((2 x n_CGs) x rank). Non-negative.

- nmf_input:

  The positive-negative split matrix fed to NMF (samples x (2 x n_CGs)).

- nmf_model:

  The [`RcppML::nmf()`](https://rdrr.io/pkg/RcppML/man/nmf.html) model
  object (includes the scaling vector `d`).

- patient_mixture:

  Long-format data frame of per-sample factor proportions used for the
  mixture plot.

## Details

Signed CCA scores are decomposed as:

- `score_pos = max(score, 0)`: patient aligned with TF program

- `score_neg = max(-score, 0)`: patient anti-aligned with TF program

Both are concatenated column-wise before NMF. Column names are suffixed
with "\_pos" and "\_neg" to track direction.
