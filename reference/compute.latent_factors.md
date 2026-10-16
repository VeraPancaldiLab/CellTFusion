# Compute latent factors from cell group scores using NMF

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
  (samples x cell groups). nmf_modelThe
  [`RcppML::nmf()`](https://rdrr.io/pkg/RcppML/man/nmf.html) model
  object (includes the scaling vector `d`). patient_mixtureLong-format
  data frame of per-sample factor proportions used for the mixture plot:
  the factor scores multiplied by the size of each factor
  (`nmf_model$d`), divided by their sum in each sample. stabilityOne row
  per factor of the selected rank: `recovery` (proportion of bootstrap
  runs in which a matching factor is found, cosine similarity \>=
  `similarity_thresh`), `mean_similarity` (mean cosine similarity of the
  matched bootstrap factors), `max_similarity_other` (highest cosine
  similarity with another factor of the same rank: a high value suggests
  redundancy) and `robust` (recovery \>= `recovery_thresh` and mean
  similarity \>= `similarity_thresh`). rank_selectionOne row per
  candidate rank (`NULL` if `rank` is given): the number of robust
  factors, the lowest recovery and mean similarity of its factors, the
  highest similarity between two of its factors, and `stable_distinct`
  (all factors robust and none redundant).

Decomposes signed CCA composite scores into non-negative latent factors
by splitting each score into its positive and negative direction
components, then applying Non-negative Matrix Factorization via RcppML
(fast C++ backend). Only the positive part of the signed cell group
scores is used: `max(score, 0)`, how far the sample is on the active
side of the TF module of the cell group (negative scores, on the
inactive side, become 0). The factors therefore describe the
high-activity direction of the cell groups. The scores must not contain
missing values.Rank selection by bootstrap NMF (consensus of the factors
across resamples of the patients). The profile of a factor is its column
of `W` (the weights of the cell groups). For each candidate rank:

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
persist at the neighbouring ranks. data("network.tuto")
data("deconv_subgroups.tuto") cell_groups \<-
construct_cell_groups(network.tuto, deconv_subgroups.tuto) latent \<-
compute.latent_factors(cell_groups\$Cell_groups, n_boot = 20, return =
FALSE) head(latent\$Z) latent\$nmf_model\$d \# size of each factor
latent\$rank_selection \# stability and redundancy of each candidate
rank latent\$stability \# stability of the factors of the selected rank
