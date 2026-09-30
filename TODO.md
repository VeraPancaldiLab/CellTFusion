# TODO

## CCA composite scores are not robust to sample order

Found during the September 2026 audit; fix deferred. Both problems are
in
[`compute_composite_score()`](https://verapancaldilab.github.io/CellTFusion/reference/compute_composite_score.md)
(`R/CellTFusion.R`) and exist in the original code. They cause no
errors, and a run with the same samples in the same order always gives
the same results, but some cell group scores change when the sample
order or the sample subset changes.

### 1. The sign of composite scores is arbitrary

Orient composite scores consistently

**What happens.** CCA weights are only defined up to a sign, and the
sign returned by
[`stats::cancor()`](https://rdrr.io/r/stats/cancor.html) depends on the
order of the samples. The composite score
(`scaled_obj %*% cca_result$xcoef[, 1]`) is used as returned, so the
same cell group can get the opposite sign when the samples are given in
another order or subset (e.g. cross-validation folds).

**Evidence.** Original code, tutorial data (192 samples, two cohorts),
same samples grouped by cohort in two orders (the second reversed within
each cohort): of 161 cell groups, 86 had identical scores, 61 were
sign-flipped (correlation −1) and 14 were unrelated (see problem 2).

**Impact.** A flipped score moves its values between the `_pos` and
`_neg` NMF inputs, so latent factors, and which factor captures a cell
group, can change between runs that differ only in sample order or
subset. The direction (positive/negative) of a cell group’s association
with its TF module is not stable.

**Possible fix (about 2 lines).** Flip `xcoef` when needed so that each
score correlates with its TF module eigengene in the direction of the
association (`signs[1]`).

### 2. Scores are arbitrary when there are as many variables as samples

Handle cell groups whose CCA has at least as many variables as samples

**What happens.** When a cell group’s deconvolution features plus the
TFs of its module reach the number of samples, the CCA fits perfectly
(first canonical correlation = 1) and its weights are not unique, so the
score is decided by numerical noise.

**Evidence.** `Dendrogram_turquoise_group_2` (10 features; turquoise
module with 184 TFs; 192 samples): first canonical correlation 1.000000.
Adding noise of 1e-12 to the TF matrix changed its score completely
(correlation 0.08 with the original score).

**Impact.** These groups look strongly associated with their module but
are overfitted. They are more frequent with large modules and with small
cohorts or cross-validation training folds (e.g. LODO).

**Options (a choice of method is needed).** - Summarize the TFs of each
module with a few principal components before the CCA. - Use only the
module eigengene. - Use a regularized CCA (e.g. `CCA::rcc()`). - At
least warn about, or skip, groups where the number of variables reaches
the number of samples.

### How to reproduce

``` r

devtools::load_all()
data(tfs.tuto); data(deconv_subgroups.tuto)
n <- nrow(tfs.tuto)
b <- rep(c("A", "B"), length.out = n)

# Same samples, grouped by cohort, in two different orders
run <- function(o) {
  net <- compute.WTCNA(split(tfs.tuto[o, ], b[o]), batch = TRUE, minMod = 3, return = FALSE)
  dt <- deconv_subgroups.tuto; dt[[1]] <- dt[[1]][o, ]
  set.seed(5)
  construct_cell_groups(net, dt, batch = b[o], n_perm = 99)
}
r1 <- run(order(b))
r2 <- run(order(b, -seq_len(n)))

# Match cell groups by module and composition, then compare their scores:
# correlation 1 = identical, -1 = sign flip (problem 1), ~0 = arbitrary (problem 2)
key <- function(cg) paste(sub("^Dendrogram_(.*)_group_.*$", "\\1", names(cg[[2]])),
                          sapply(cg[[2]], function(x) paste(sort(x), collapse = "|")))
k1 <- setNames(key(r1), names(r1[[2]])); k2 <- setNames(key(r2), names(r2[[2]]))
cors <- sapply(intersect(k1, k2), function(k) {
  s1 <- r1[[1]][, names(k1)[k1 == k][1]]; names(s1) <- rownames(r1[[1]])
  s2 <- r2[[1]][, names(k2)[k2 == k][1]]; names(s2) <- rownames(r2[[1]])
  cor(s1, s2[names(s1)])
})
table(round(cors))
```
