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
  ranks = 2:15,
  n_boot = 100,
  recovery_thresh = 0.7,
  similarity_thresh = 0.8,
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

  Integer; number of NMF factors. If NULL (default), chosen among
  `ranks` with bootstrap NMF (see Details). If given, the stability of
  its factors is still computed (`stability`).

- ranks:

  Candidate ranks when `rank = NULL`. Default 2:15, limited to the
  number of cell groups (rank 1 with a single cell group).

- n_boot:

  Integer; number of bootstrap NMF runs per rank. Default 100.

- recovery_thresh:

  Numeric; minimum recovery frequency (proportion of bootstrap runs in
  which the factor is found) for a factor to be robust. Default 0.7.

- similarity_thresh:

  Numeric; minimum cosine similarity for a bootstrap factor to count as
  the same factor, and minimum mean similarity for a factor to be
  robust. Two factors of the same rank more similar than this are
  considered redundant. Default 0.8. Both thresholds are practical
  starting points, not universal cutoffs.

- seed:

  Random seed used for the NMF fits and the bootstrap samples.
  Default 123. The caller's random number generator state is restored
  when the function returns.

- file_name:

  Optional character suffix for the saved patient-mixture plot.

- return:

  Logical. If TRUE (default), saves the patient-mixture barplot to
  `Results/NMF_patient_mixture_<file_name>.pdf`
  (`Results/NMF_patient_mixture.pdf` if `file_name` is `NULL`), and the
  stability plot to `Results/NMF_rank_stability_<file_name>.pdf`: the
  lowest recovery, lowest mean similarity and highest similarity between
  factors of each candidate rank (with the thresholds and the selected
  rank; only when the rank is selected), and the recovery and mean
  similarity of each factor of the selected rank.

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

  Weights of the cell groups per factor (cell groups x rank).
  Non-negative.

- nmf_input:

  The matrix fed to NMF: the positive part of the cell group scores
  (samples x cell groups).

- nmf_model:

  The [`RcppML::nmf()`](https://rdrr.io/pkg/RcppML/man/nmf.html) model
  object (includes the scaling vector `d`).

- patient_mixture:

  Long-format data frame of per-sample factor proportions used for the
  mixture plot: the factor scores multiplied by the size of each factor
  (`nmf_model$d`), divided by their sum in each sample.

- stability:

  One row per factor of the selected rank: `recovery` (proportion of
  bootstrap runs in which a matching factor is found, cosine similarity
  \>= `similarity_thresh`), `mean_similarity` (mean cosine similarity of
  the matched bootstrap factors), `max_similarity_other` (highest cosine
  similarity with another factor of the same rank: a high value suggests
  redundancy) and `robust` (recovery \>= `recovery_thresh` and mean
  similarity \>= `similarity_thresh`).

- rank_selection:

  One row per candidate rank (`NULL` if `rank` is given): the number of
  robust factors, the lowest recovery and mean similarity of its
  factors, the highest similarity between two of its factors, and
  `stable_distinct` (all factors robust and none redundant).

## Details

Only the positive part of the signed cell group scores is used:
`max(score, 0)`, how far the sample is on the active side of the TF
module of the cell group (negative scores, on the inactive side, become
0). The factors therefore describe the high-activity direction of the
cell groups. The scores must not contain missing values.

Rank selection by bootstrap NMF (consensus of the factors across
resamples of the patients). The profile of a factor is its column of `W`
(the weights of the cell groups). For each candidate rank:

1.  NMF is fitted on all the samples (reference factors) and on `n_boot`
    bootstrap samples of the patients (sampled with replacement);

2.  the factors of each bootstrap run are matched one to one to the
    reference factors by the cosine similarity of their profiles (most
    similar pairs first), as NMF factors can come in any order;

3.  each reference factor gets its recovery frequency and mean
    similarity across runs, and the factors of the rank are compared
    with each other to detect redundancy (two factors with almost the
    same profile).

The selected rank is the largest one whose factors are all robust and
none redundant: as many factors as the data support reliably (low ranks
capture the dominant patterns and are almost always stable, so the
smallest stable rank would favour 2). If no rank qualifies, the rank
with the most robust factors is used (the smallest in case of ties),
with a message; check `rank_selection` and whether the main factors
persist at the neighbouring ranks.

## Examples

``` r
data("network.tuto")
data("deconv_subgroups.tuto")
cell_groups <- construct_cell_groups(network.tuto, deconv_subgroups.tuto)
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> ℹ The deprecated feature was likely used in the dendextend package.
#>   Please report the issue at <https://github.com/talgalili/dendextend/issues>.
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Dendrogram plots saved to: Results/Dendrogram_color_clusters.pdf
#> 44 of 897 features were not assigned to any cluster by the tree cut and are not part of any cell group (pink: 2, yellow: 8, red: 2, blue: 6, brown: 4, turquoise: 12, black: 4, green: 6)
latent <- compute.latent_factors(cell_groups$Cell_groups, n_boot = 20, return = FALSE)
#> Selecting the rank with bootstrap NMF (20 runs per rank, ranks 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)...
#> Selected rank: 5 (largest rank with all factors robust and none redundant)
#> Running RcppML NMF with rank = 5...
head(latent$Z)
#>                      Factor1     Factor2      Factor3      Factor4      Factor5
#> SAM7f0d9cc7f001 0.0000000000 0.000000000 0.0160659311 0.0067204288 1.237777e-02
#> SAM4305ab968b90 0.0000000000 0.002928540 0.0008473742 0.0000000000 0.000000e+00
#> SAMcf018fee2acd 0.0000000000 0.003929063 0.0142868585 0.0052293929 6.173109e-03
#> SAMcc4675f394a1 0.0013695473 0.000000000 0.0161753036 0.0000000000 3.847780e-03
#> SAM49f9b2e57aa5 0.0006725578 0.003517198 0.0017914287 0.0000000000 3.909520e-03
#> SAM2e7aa8fa0ab3 0.0006385870 0.005290507 0.0026158628 0.0002566458 4.774013e-05
latent$nmf_model$d     # size of each factor
#> [1] 220.7830 214.0921 191.4122 181.4714 118.0781
latent$rank_selection  # stability and redundancy of each candidate rank
#>    rank robust_factors min_recovery min_mean_similarity
#> 1     2              2         1.00           0.9621334
#> 2     3              3         0.90           0.9300359
#> 3     4              4         0.75           0.8717509
#> 4     5              5         0.85           0.8314428
#> 5     6              5         0.35           0.5643514
#> 6     7              4         0.35           0.6243858
#> 7     8              4         0.20           0.6317879
#> 8     9              4         0.30           0.6276593
#> 9    10              6         0.30           0.6813313
#> 10   11              4         0.35           0.6818662
#> 11   12              5         0.30           0.6606198
#> 12   13              6         0.05           0.4379847
#> 13   14              5         0.35           0.5666457
#> 14   15              7         0.05           0.5939563
#>    max_similarity_between_factors stable_distinct
#> 1                       0.4651282            TRUE
#> 2                       0.4458516            TRUE
#> 3                       0.5357415            TRUE
#> 4                       0.5297968            TRUE
#> 5                       0.5842894           FALSE
#> 6                       0.5043616           FALSE
#> 7                       0.5312403           FALSE
#> 8                       0.4934716           FALSE
#> 9                       0.5777689           FALSE
#> 10                      0.4718015           FALSE
#> 11                      0.4859773           FALSE
#> 12                      0.4562006           FALSE
#> 13                      0.5400378           FALSE
#> 14                      0.4882961           FALSE
latent$stability       # stability of the factors of the selected rank
#>    factor recovery mean_similarity max_similarity_other robust
#> 1 Factor1     0.95       0.8796546            0.5065433   TRUE
#> 2 Factor2     1.00       0.9458671            0.3888573   TRUE
#> 3 Factor3     0.85       0.8314428            0.5065433   TRUE
#> 4 Factor4     0.95       0.9215742            0.5297968   TRUE
#> 5 Factor5     0.90       0.9253520            0.5297968   TRUE
```
