# TODO

## CCA composite scores are not robust to sample order

Found during the September 2026 audit; fix deferred. Both problems are
in
[`compute_composite_score()`](https://verapancaldilab.github.io/CellTFusion/reference/compute_composite_score.md)
(`R/CellTFusion.R`) and exist in the original code. They cause no
errors, and a run with the same samples in the same order always gives
the same results, but some cell group scores change when the sample
order or the sample subset changes.

**Status (5 October 2026).** Postponed again: to be dealt with when the
review reaches
[`compute_composite_score()`](https://verapancaldilab.github.io/CellTFusion/reference/compute_composite_score.md).
Measured again on the code of that day (tutorial data, single cohort,
123 cell groups, 192 samples):

- *Sign.* Same samples in original and reversed order: 72 groups with
  the same score, 50 with the opposite sign, 1 different. Of 66 groups
  whose features are positively related to their module, 30 have a score
  negatively related to it; of 57 negatively related groups, 32 have a
  score positively related.

- *Many TFs.* First canonical correlation with the real features against
  random numbers in their place: purple (26 TFs) 0.85 / 0.62; pink (29)
  0.87 / 0.47; red (74) 0.92 / 0.68; brown (92) 0.96 / 0.79;
  turquoise (181) 1.00 / 1.00. In turquoise, noise of 0.01 on the TF
  activities changes the score (correlation 0.11 with the original); in
  purple it does not (1.00). 6 of 123 groups have a canonical
  correlation above 0.99, all in turquoise.

- *What the sign changes downstream (measured 5 October 2026, correcting
  the “Impact” paragraph of problem 1 below).* Latent factors and niches
  did **not** change when 50 of 123 scores flipped: the three factors of
  the two runs correlate 1.00 and the niches are identical, because a
  flip only swaps the `_pos` and `_neg` columns of that group in the NMF
  input. What changes is the cell group score table itself and which of
  `_pos` / `_neg` a group appears under in `W`.

- *How a score could be oriented.* Its link with the module score is too
  weak in some groups to decide the sign reliably (absolute correlation:
  minimum 0.01, 10% of groups below 0.18). Its link with the mean of the
  group’s standardized features is stronger (minimum 0.14, 10% below
  0.47).

### 1. The sign of composite scores is arbitrary

Orient composite scores consistently — **fixed 5–6 October 2026.** After
the CCA,
[`compute_composite_score()`](https://verapancaldilab.github.io/CellTFusion/reference/compute_composite_score.md)
flips the weights and the score when needed so that the score is aligned
with the TF module: it increases with the features of the cell group
when they are positively related to the module and decreases with them
when negatively related (a positive score = the sample is on the active
side of the module). The `* -1` step was removed. The relation of each
group to its module was returned in a `Module_relation` table of
[`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md)
(removed 14 October 2026 at the user’s request). Tutorial data: the
score follows this rule in all 123 groups (66 positive, 57 negative);
122 of 123 scores are identical when the samples are reordered (72
before); 121 of 123 scores are positively related to their module score
(the two others have almost no link with it, \|r\| 0.11 and 0.05). The
text below describes the behaviour before the fix.

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
  construct_cell_groups(net, dt, batch = b[o])
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

## Idea: cluster features on their correlation with each TF of the module

Not a bug; kept as an idea for the future (October 2026). The aim of a
cell group is to gather features with a similar biological activity.
[`identify.cell.groups()`](https://verapancaldilab.github.io/CellTFusion/reference/identify.cell.groups.md)
measures that with a single number per feature: its correlation with the
module score. Two features with the same value can still be related to
different TFs of that module.

Consider clustering the features of each module on their correlations
with each TF of the module (one value per TF instead of one value per
module), a finer measure of the same idea. It would change the cell
groups and everything computed from them.

## Multiple testing in `compute.modules.relationship()` (noted 11 October 2026)

To deal with later. The `padj` option was removed, so the p-values are
compared with `pval` without any correction.
[`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md)
uses these p-values to choose the deconvolution features of each TF
module.

Measured on the tutorial data (10 TF modules x 184 deconvolution
features = 1,840 tests, `pval = 0.05`):

| `matB` | p \<= 0.05 | Benjamini-Hochberg \<= 0.05 |
|----|----|----|
| real deconvolution features | 764 | 614 |
| 184 random columns (no real relation) | 101 | 0 |

With the random columns, all 10 modules had at least two “significant”
features, which is enough to build cell groups.

Decide where to correct the p-values (in
[`compute.modules.relationship()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.modules.relationship.md)
for every use, or only in
[`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md)),
with which method, and across which tests (all pairs, or per module). It
would change the cell groups and everything computed from them. **Done
14 October 2026:** `padj` (TRUE / FALSE, default FALSE) in
[`compute.modules.relationship()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.modules.relationship.md),
[`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md)
and
[`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md);
Benjamini-Hochberg across the features for each TF module; the adjusted
p-values are returned as `p_adjusted`. With the default, results are
unchanged.

Decide whether `padj = TRUE` should become the default (changes the cell
groups for every caller).

## Scaling of TF and pathway activity scores (noted 12 October 2026)

To check later: whether the scores should be scaled, and why it is a
problem. It concerns both
[`compute.TFs.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.TFs.activity.md)
and
[`compute.pathway.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.pathway.activity.md),
which z-score every TF / pathway across the samples of the call
(`scale = TRUE`, the default).

What was measured on the tutorial data:

The score of a sample depends on the other samples given in the same
call. The same samples scored alone or within a larger set keep their
ranking, but the values differ (by up to 2.12 for TFs with 20 vs 60
samples, by up to 3.03 for pathways with 30 vs 192 samples). With
`scale = FALSE` the values are identical.

With one sample the scaled values were all `NaN`, and with two samples
all +0.707 / -0.707. Since 12 October 2026 the scores are not scaled
below three samples (a message is given); the threshold of three was
chosen without a specific reason.

Scores computed in separate calls (training and test sets, or several
cohorts) are therefore not on the same scale. `scripts/src/utils.R` of
the paper does this for training and test sets.

[`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md)
and
[`compute_composite_score()`](https://verapancaldilab.github.io/CellTFusion/reference/compute_composite_score.md)
assume scaled TF activities (comments in the code: module eigengenes are
computed with `scale = F` because the TF scores are already scaled; the
CCA uses the TF scores as they are).

Decide whether `scale = TRUE` stays the default, what to do for test
sets and multi-cohort runs (TF activity is computed and scaled per
cohort in batch mode), and whether the downstream steps should scale by
themselves instead.

## Consistency of parameters shared between functions (noted 14 October 2026)

To check later: the same setting has different names, values and
defaults in different functions, and some are hard-coded in one place
and settable in another. Inventory of `R/CellTFusion.R` on 14 October
2026 (`*` = internal function):

**Correlation type.**

- Names: `cor_type`
  ([`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md),
  [`compute.modules.relationship()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.modules.relationship.md),
  [`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md)),
  `corr_type`
  ([`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md),
  passed to
  [`multideconv::compute.deconvolution.analysis()`](https://verapancaldilab.github.io/multideconv/reference/compute.deconvolution.analysis.html)),
  `corr_method`
  ([`compute.metadata.association()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.metadata.association.md)),
  `method` (`mergeModules()*`, where `method` means the scoring
  statistic in
  [`compute.TFs.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.TFs.activity.md)
  /
  [`compute.pathway.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.pathway.activity.md)).
- Values: `"p"` / `"s"` in most functions, but `"spearman"` /
  `"pearson"` in `corr_type` and `mergeModules()*`.
- Defaults: Pearson everywhere except the deconvolution subgroups
  (Spearman).
- Hard-coded:
  - [`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)
    passes `cor_type = "p"` to
    [`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md).
    The `cor_type` of
    [`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md)
    (added 14 October 2026) is exposed in
    [`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)
    as `cor_type_groups`, next to `corr_type` (deconvolution subgroups).
  - The ARACNe signs without a `mor` column use Spearman.
- Threshold names: `corr`
  ([`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md):
  deconvolution subgroups; `mergeModules()*`: module merging),
  `corr_mod`
  ([`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md)).

**Tree cut (deep split).**

- [`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md)
  hard-codes `deepSplit = 2`, the hybrid method’s 0-4 scale, in both
  single- and multi-cohort mode.
- The cell groups always use `deep_split = TRUE`, a switch in
  `method = "tree"`. Since 14 October 2026 it is no longer a parameter
  of
  [`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md),
  [`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md)
  or
  [`cell.groups.computation()`](https://verapancaldilab.github.io/CellTFusion/reference/cell.groups.computation.md);
  `calculate_dendrogram_cuts()*` keeps it.
- Decide whether this stays hard-coded in both places, and document that
  the two values have different meanings.

**Minimum sizes.**

- `minMod` (TF modules): default 15 in
  [`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md)
  but 5 in
  [`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md).
- `min_cluster_size` (cell groups), `min_targets_size` (regulons),
  `min_genes` (`module_enrich()*`).

**Clustering method.** `clustering.method` is `"average"` in
[`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md)
and `"ward.D2"` in
[`identify.cell.groups()`](https://verapancaldilab.github.io/CellTFusion/reference/identify.cell.groups.md)
/
[`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md);
[`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)
hard-codes both.

**Other common parameters.**

- **Output file:** `file_name` vs `file.name` vs `dendrogram_file`.
- **Saving and plots:**
  - `return` defaults to `TRUE` in most functions but `FALSE` in
    [`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md)
    and
    [`prepare_celltfusion_folds()`](https://verapancaldilab.github.io/CellTFusion/reference/prepare_celltfusion_folds.md).
  - Plots are switched by `plot`, `plot_dot`, `plot_grid`, `plot.grid`
    or `return_dendrogram`.
- **`verbose`:** default `TRUE` in
  [`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md),
  `FALSE` in
  [`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md)
  and
  [`prepare_celltfusion_folds()`](https://verapancaldilab.github.io/CellTFusion/reference/prepare_celltfusion_folds.md).
- **`batch`:** a cohort vector everywhere except
  [`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md)
  (logical).
- **Multiple-testing correction:** `p_adjust` (any
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) method,
  default `"none"`) in
  [`compute.metadata.association()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.metadata.association.md),
  [`compute.survival.analysis()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.survival.analysis.md)
  and the `scores.*` tests, but `padj` (TRUE / FALSE, always
  Benjamini-Hochberg) in
  [`compute.modules.relationship()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.modules.relationship.md),
  [`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md)
  and
  [`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)
  since 14 October 2026.
- **`seed`:** only in
  [`compute.latent_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.latent_factors.md)
  and
  [`compute_factor_gsea()`](https://verapancaldilab.github.io/CellTFusion/reference/compute_factor_gsea.md).

Renaming or changing defaults breaks callers (MCP server, paper scripts,
Shiny app): list them before deciding.

## Deferred: to do after the review of `CellTFusion()` (noted October 2026)

**TCGA references.** Regenerate
`inst/extdata/TCGA_meta_programs_*.RData` and
`CellTFusion_TCGA_*.RData`: they were built while cell groups were
mislabelled (2026-04-29 to 2026-10-05), so `TME_states` maps against
references from wrong cell groups. They must also be rebuilt with the
current defaults (12 October 2026): ulm TF activities, TF tree on the
TOM dissimilarity with average linkage, merging repeated until no pair
is above `corr_mod`, cell groups cut into single ranges of correlation,
rank limited to the number of cell groups, GSEA with the cohort
covariate and a fixed seed. The current blca reference labels 5 of its 7
meta-programs as IE.

**Shiny app.** Launch and test it after the changes to `batch`,
`normalized` and the removed arguments; remove the unused
“High-correlation cutoff” slider; update the “Minimum TF targets”
default (3, package default is 5) and the “Normalize counts (TPM +
log2)” label (unticked now means the input is TPM).

**Vignettes.** Regenerate the figures of cell groups, latent factors,
niches and TME states (made before the cell group fix).

**Data.** Regenerate `deconv.tuto` with multideconv \>= 0.2.0 naming;
add survival columns to the tutorial data, or a small survival example,
so the survival examples can run. `tfs.tuto` and `network.tuto` were
computed before the default scoring method of
[`compute.TFs.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.TFs.activity.md)
changed from the decoupleR consensus to `"ulm"` (11 October 2026) and
before
[`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md)
built its TF tree on the TOM dissimilarity with average linkage (12
October 2026): regenerate them too.

**Docs.** Review them again once the items above are done.

**Tests.** Add unit tests (`tests/testthat/` does not exist yet).

**Plot of cell group scores with their relation to the TF module.** To
decide later (noted 9 October 2026). Four alternatives were drafted and
discarded for now: (A) heatmap of all cell groups x samples annotated
with module and relation; (B) one page per module, samples ordered by
module activity, one row per cell group labelled “more”/“fewer” with its
cell types; (C) cell groups as points left/right of zero by their
correlation with the module, labelled with cell types; (D) one sample,
its strongest cell groups per module coloured by more/fewer cells.

**Scoring method in
[`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md).**
[`compute.TFs.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.TFs.activity.md)
and
[`compute.pathway.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.pathway.activity.md)
have a `method` argument (default `"ulm"`), but
[`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)
does not pass one, so the pipeline always uses `"ulm"`. Decide whether
to add the argument to
[`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)
(noted 11 October 2026).

Check behavior across all pipeline when batch is set to TRUE (noted 12
October 2026).

**Revise
[`calculate_dendrogram_cuts()`](https://verapancaldilab.github.io/CellTFusion/reference/calculate_dendrogram_cuts.md)
to really understand how it works** (noted 12 October 2026): how
`dynamicTreeCut::cutreeDynamic(method = "tree")` gives out cluster
numbers (on the tutorial data, 55 of 123 cell groups were made of
separated ranges of correlation), and whether the step added that splits
clusters into single ranges and joins small pieces to their closest
neighbour is the right answer.

**Cell niches: more or fewer cells.**
[`compute_cells_niches()`](https://verapancaldilab.github.io/CellTFusion/reference/compute_cells_niches.md)
adds together the selected cell group parts that mean more cells and
those that mean fewer (it uses neither `_pos` / `_neg` nor the relation
of each group to its TF module). Measured 12 October 2026 on the
tutorial data (3 factors): the selected parts meaning more / fewer cells
are 41 / 25, 29 / 40 and 38 / 17; the top cell type of Factor1, CAF
(weight 0.14), gets 0.11 from parts meaning fewer CAFs and 0.03 from
parts meaning more. Proposed: give each part its direction (from the
relation of its group to the TF module; `Module_relation`, which gave
it, was removed on 14 October 2026, so it would have to be computed
here), compute weights and enrichment separately for more and fewer,
return per factor a table (cell type, direction, weight, enrichment) and
colour the plot edges by direction; this changes the structure of
`Cells_niches`.
