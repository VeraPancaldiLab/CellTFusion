# Compute latent factors from cell group scores using NMF

Decomposes signed CCA composite scores into non-negative latent factors
by splitting each score into its positive and negative direction
components, then applying Non-negative Matrix Factorization via RcppML
(fast C++ backend).

## Usage

``` r
compute.latent_factors(
  scores,
  rank = NULL,
  seed = 123,
  file_name = NULL,
  return = TRUE
)
```

## Arguments

- scores:

  Numeric matrix of size samples x cell groups (signed CCA composite
  scores).

- rank:

  Integer; number of NMF factors. If NULL, estimated automatically at
  the elbow of the reconstruction MSE across ranks 2:8 (the rank with
  the largest second difference of the MSE curve). Ranks are limited to
  the number of cell groups (no more factors than cell groups); with
  fewer than 4 cell groups there are too few candidate ranks for an
  elbow, and rank = 2 is used (rank = 1 with a single cell group).

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

  Sample-level NMF factor scores (samples x rank). Non-negative. Each
  factor (column) sums to 1 over the samples: the size of each factor is
  in `nmf_model$d` (the NMF is `W %*% diag(d) %*% t(Z)`). Tests on one
  factor at a time do not depend on this scale; to compare factors
  within a sample, multiply each column by its `d`, as done for
  `patient_mixture`.

- W:

  Feature weights per factor ((2 x n_CGs) x rank). Non-negative.

- nmf_input:

  The positive-negative split matrix fed to NMF (samples x (2 x n_CGs)).

- nmf_model:

  The [`RcppML::nmf()`](https://rdrr.io/pkg/RcppML/man/nmf.html) model
  object (includes the scaling vector `d`).

- patient_mixture:

  Long-format data frame of per-sample factor proportions used for the
  mixture plot: the factor scores multiplied by the size of each factor
  (`nmf_model$d`), divided by their sum in each sample.

## Details

Signed CCA scores are decomposed as:

- `score_pos = max(score, 0)`: how far the sample is on the active side
  of the TF module of that cell group (aligned with the TF program)

- `score_neg = max(-score, 0)`: how far the sample is on the inactive
  side (anti-aligned)

Both are concatenated column-wise before NMF. Column names are suffixed
with "\_pos" and "\_neg" to track direction. The scores must not contain
missing values.

## Examples

``` r
data("network.tuto")
data("deconv_subgroups.tuto")
cell_groups <- construct_cell_groups(network.tuto, deconv_subgroups.tuto)
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Dendrogram plots saved to: Results/Dendrogram_color_clusters.pdf
#> 44 of 897 features were not assigned to any cluster by the tree cut and are not part of any cell group (pink: 2, yellow: 8, red: 2, blue: 6, brown: 4, turquoise: 12, black: 4, green: 6)
latent <- compute.latent_factors(cell_groups$Cell_groups, return = FALSE)
#> Estimating optimal rank via elbow on reconstruction MSE (ranks 2:8)...
#> MSE by rank: k2=0.00129  k3=0.00114  k4=0.00105  k5=0.000978  k6=0.000933  k7=0.000893  k8=0.000857
#> Selected rank (elbow): 3
#> Running RcppML NMF with rank = 3...
head(latent$Z)
#>                      Factor1     Factor2     Factor3
#> SAM7f0d9cc7f001 0.0001788103 0.008323831 0.005546469
#> SAM4305ab968b90 0.0194986831 0.000000000 0.000000000
#> SAMcf018fee2acd 0.0018946456 0.006329754 0.005301130
#> SAMcc4675f394a1 0.0047871944 0.008222151 0.000000000
#> SAM49f9b2e57aa5 0.0086963389 0.000615321 0.001548553
#> SAM2e7aa8fa0ab3 0.0046124714 0.001796520 0.001894464
latent$nmf_model$d # size of each factor
#> [1] 811.4351 642.3944 383.9934
```
