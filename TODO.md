# TODO

## CCA composite scores are not robust to sample order

Found during the September 2026 audit; fix deferred. Both problems are in `compute_composite_score()` (`R/CellTFusion.R`) and exist in the original code. They cause no errors, and a run with the same samples in the same order always gives the same results, but some cell group scores change when the sample order or the sample subset changes.

**Status (5 October 2026).** Postponed again: to be dealt with when the review reaches `compute_composite_score()`. Measured again on the code of that day (tutorial data, single cohort, 123 cell groups, 192 samples):

- *Sign.* Same samples in original and reversed order: 72 groups with the same score, 50 with the opposite sign, 1 different. Of 66 groups whose features are positively related to their module, 30 have a score negatively related to it; of 57 negatively related groups, 32 have a score positively related.
- *Many TFs.* First canonical correlation with the real features against random numbers in their place: purple (26 TFs) 0.85 / 0.62; pink (29) 0.87 / 0.47; red (74) 0.92 / 0.68; brown (92) 0.96 / 0.79; turquoise (181) 1.00 / 1.00. In turquoise, noise of 0.01 on the TF activities changes the score (correlation 0.11 with the original); in purple it does not (1.00). 6 of 123 groups have a canonical correlation above 0.99, all in turquoise.

- *What the sign changes downstream (measured 5 October 2026, correcting the "Impact" paragraph of problem 1 below).* Latent factors and niches did **not** change when 50 of 123 scores flipped: the three factors of the two runs correlate 1.00 and the niches are identical, because a flip only swaps the `_pos` and `_neg` columns of that group in the NMF input. What changes is the cell group score table itself and which of `_pos` / `_neg` a group appears under in `W`.
- *How a score could be oriented.* Its link with the module score is too weak in some groups to decide the sign reliably (absolute correlation: minimum 0.01, 10% of groups below 0.18). Its link with the mean of the group's standardized features is stronger (minimum 0.14, 10% below 0.47).

### 1. The sign of composite scores is arbitrary

- [x] Orient composite scores consistently — **fixed 5–6 October 2026.** After the CCA, `compute_composite_score()` flips the weights and the score when needed so that the score is aligned with the TF module: it increases with the features of the cell group when they are positively related to the module and decreases with them when negatively related (a positive score = the sample is on the active side of the module). The `* -1` step was removed. The relation of each group to its module was returned in a `Module_relation` table of `construct_cell_groups()` (removed 14 October 2026 at the user's request). Tutorial data: the score follows this rule in all 123 groups (66 positive, 57 negative); 122 of 123 scores are identical when the samples are reordered (72 before); 121 of 123 scores are positively related to their module score (the two others have almost no link with it, |r| 0.11 and 0.05). The text below describes the behaviour before the fix.

**What happens.** CCA weights are only defined up to a sign, and the sign returned by `stats::cancor()` depends on the order of the samples. The composite score (`scaled_obj %*% cca_result$xcoef[, 1]`) is used as returned, so the same cell group can get the opposite sign when the samples are given in another order or subset (e.g. cross-validation folds).

**Evidence.** Original code, tutorial data (192 samples, two cohorts), same samples grouped by cohort in two orders (the second reversed within each cohort): of 161 cell groups, 86 had identical scores, 61 were sign-flipped (correlation −1) and 14 were unrelated (see problem 2).

**Impact.** A flipped score moves its values between the `_pos` and `_neg` NMF inputs, so latent factors, and which factor captures a cell group, can change between runs that differ only in sample order or subset. The direction (positive/negative) of a cell group's association with its TF module is not stable.

**Possible fix (about 2 lines).** Flip `xcoef` when needed so that each score correlates with its TF module eigengene in the direction of the association (`signs[1]`).

### 2. Scores are arbitrary when there are as many variables as samples

- [x] Handle cell groups whose CCA has at least as many variables as samples — **largely resolved 15 October 2026:** the CCA is now between the features of the cell group and the TF module score only (not the TFs of the module), so it has (number of features + 1) variables instead of (features + TFs of the module); it is the linear regression of the module score on the features. A cell group with as many features as samples would still be degenerate (rare).

**What happens.** When a cell group's deconvolution features plus the TFs of its module reach the number of samples, the CCA fits perfectly (first canonical correlation = 1) and its weights are not unique, so the score is decided by numerical noise.

**Evidence.** `Dendrogram_turquoise_group_2` (10 features; turquoise module with 184 TFs; 192 samples): first canonical correlation 1.000000. Adding noise of 1e-12 to the TF matrix changed its score completely (correlation 0.08 with the original score).

**Impact.** These groups look strongly associated with their module but are overfitted. They are more frequent with large modules and with small cohorts or cross-validation training folds (e.g. LODO).

**Options (a choice of method is needed).**
- Summarize the TFs of each module with a few principal components before the CCA.
- Use only the module eigengene.
- Use a regularized CCA (e.g. `CCA::rcc()`).
- At least warn about, or skip, groups where the number of variables reaches the number of samples.

### How to reproduce

```r
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

Not a bug; kept as an idea for the future (October 2026). The aim of a cell group is to gather features with a similar biological activity. `identify.cell.groups()` measures that with a single number per feature: its correlation with the module score. Two features with the same value can still be related to different TFs of that module.

- [x] (Closed 15 October 2026, not pursued.) Consider clustering the features of each module on their correlations with each TF of the module (one value per TF instead of one value per module), a finer measure of the same idea. It would change the cell groups and everything computed from them.

## Multiple testing in `compute.modules.relationship()` (noted 11 October 2026)

To deal with later. The `padj` option was removed, so the p-values are compared with `pval` without any correction. `construct_cell_groups()` uses these p-values to choose the deconvolution features of each TF module.

Measured on the tutorial data (10 TF modules x 184 deconvolution features = 1,840 tests, `pval = 0.05`):

| `matB` | p <= 0.05 | Benjamini-Hochberg <= 0.05 |
| --- | --- | --- |
| real deconvolution features | 764 | 614 |
| 184 random columns (no real relation) | 101 | 0 |

With the random columns, all 10 modules had at least two "significant" features, which is enough to build cell groups.

- [x] Decide where to correct the p-values (in `compute.modules.relationship()` for every use, or only in `construct_cell_groups()`), with which method, and across which tests (all pairs, or per module). It would change the cell groups and everything computed from them. **Done 14 October 2026:** `padj` (TRUE / FALSE, default FALSE) in `compute.modules.relationship()`, `construct_cell_groups()` and `CellTFusion()`; Benjamini-Hochberg across the features for each TF module; the adjusted p-values are returned as `p_adjusted`. With the default, results are unchanged.
- [x] Decide whether `padj = TRUE` should become the default: **no** (user, 15 October 2026); it stays FALSE.

## Scaling of TF and pathway activity scores (noted 12 October 2026)

To check later: whether the scores should be scaled, and why it is a problem. It concerns both `compute.TFs.activity()` and `compute.pathway.activity()`, which z-score every TF / pathway across the samples of the call (`scale = TRUE`, the default).

What was measured on the tutorial data:

- The score of a sample depends on the other samples given in the same call. The same samples scored alone or within a larger set keep their ranking, but the values differ (by up to 2.12 for TFs with 20 vs 60 samples, by up to 3.03 for pathways with 30 vs 192 samples). With `scale = FALSE` the values are identical.
- With one sample the scaled values were all `NaN`, and with two samples all +0.707 / -0.707. Since 12 October 2026 the scores are not scaled below three samples (a message is given); the threshold of three was chosen without a specific reason.
- Scores computed in separate calls (training and test sets, or several cohorts) are therefore not on the same scale. `scripts/src/utils.R` of the paper does this for training and test sets.
- `compute.WTCNA()` assumes scaled TF activities (comment in the code: module eigengenes are computed with `scale = F` because the TF scores are already scaled). Since 15 October 2026 `compute_composite_score()` no longer uses the TF activities (its CCA uses the module score).

- [x] Decide whether `scale = TRUE` stays the default, what to do for test sets and multi-cohort runs, and whether the downstream steps should scale by themselves instead. **Done 15 October 2026:** `CellTFusion()` has no `scale` argument (user's decision): it returns raw TF and pathway scores in single-cohort mode (they do not depend on the other samples, e.g. for ML features scaled with the training parameters) and scores z-scored within each cohort with `batch` (correction of cohort shifts). `compute.TFs.activity()` / `compute.pathway.activity()` keep their own `scale` argument (default TRUE). `compute.WTCNA()` standardizes the TFs itself for the module eigengenes and `identify_hub_TFs()` standardizes them within each cohort for the membership, so modules, hub TFs, cell groups and factors are identical with either setting (checked on the tutorial data, single- and multi-cohort).
- [x] (Closed 15 October 2026, not pursued.) Optional: a pipeML fold construction function, or a reference-scaling option, to scale TF / pathway scores used directly as ML features with the training samples of each fold.

## Consistency of parameters shared between functions (noted 14 October 2026)

**Done 15 October 2026** (user's decisions):

- [x] **Correlation type.** One argument everywhere, `corr_type`, with the values `"spearman"` (default everywhere) or `"pearson"`: `compute.WTCNA()`, `compute.modules.relationship()`, `construct_cell_groups()`, `compute.metadata.association()` (was `corr_method`), `mergeModules()*` (was `method`). In `CellTFusion()`, the single `corr_type` drives the deconvolution subgroups, the TF network and the cell groups (`cor_type_groups` / `cor_type_modules` removed). The TF modules and cell groups change with the new Spearman default (`corr_type = "pearson"` reproduces the previous cell groups). Spearman is available in multi-cohort mode too.
- [x] **Tree cut (deep split).** Never a user parameter: `deepSplit = 2` (hybrid cut of the TF modules) and `deepSplit = TRUE` (tree cut of the cell groups) are hard-coded where `cutreeDynamic()` is called; `calculate_dendrogram_cuts()*` no longer has the argument.
- [x] **Minimum sizes.** `minMod` renamed `min_module_size`, default 5 in `compute.WTCNA()` (was 15) and `CellTFusion()`; `min_cluster_size` and `min_targets_size` unchanged.
- [x] **Clustering method.** One argument, `clustering_method`, default `"ward.D2"` everywhere (`compute.WTCNA()` was `"average"`, so the TF modules change), exposed in `CellTFusion()` and used for both the TF modules and the cell groups.
- [x] **Multiple-testing correction.** `padj` (TRUE / FALSE, Benjamini-Hochberg, default FALSE) everywhere: also `compute.metadata.association()`, `compute.survival.analysis()`, the `scores.*` tests and `scores.stat.analysis()` (was `p_adjust`).
- [x] **Setting arguments in snake_case.** `file.name` and `dendrogram_file` -> `file_name`; `cancer.type` and `cancer_name` (`map_factors_to_TME()`) -> `cancer_type`; `softPower` -> `soft_power`; `network.type` -> `network_type`; `plot.grid` -> `plot_grid`; `ncol.grid` and `ncol` (`compute.metadata.association()`) -> `ncol_grid`; `plot_dot` -> `plot` (`compute_factor_gsea()`); `MM_thresh` -> `mm_thresh`; `cbsx.mail` / `cbsx.token` -> `cbsx_mail` / `cbsx_token`; `doParallel` -> `do_parallel` (`CellTFusion()`); `p.value` -> `pval` (`compute.survival.analysis()`). Data inputs (`RNA.tpm`, `TFs.matrix`, `raw.counts`, ...) and `TF.collection` keep their names. The vignettes, README and Shiny app were updated.
- [x] **Scoring method in `CellTFusion()`.** `method` (default `"ulm"`) passed to `compute.TFs.activity()` and `compute.pathway.activity()`.

Done 15 October 2026 (second pass):

- [x] **Saving:** one switch `return` (default `TRUE`) everywhere saves a function's tables and figures. `plot` was renamed `return` in `compute.modules.relationship()`, `compute_factor_gsea()`, `derive_meta_programs()`, `map_factors_to_metaprograms()` and `map_factors_to_TME()`; `return_dendrogram` was folded into `return` (`construct_cell_groups()`, `cell.groups.computation()`); `construct_cell_groups()` and `prepare_celltfusion_folds()` now default to `TRUE`. Optional extra figures keep `plot_grid`.
- [x] **`verbose`:** default `TRUE` everywhere (`CellTFusion()`, `compute.WTCNA()`, `prepare_celltfusion_folds()`).
- [x] **`batch`:** a cohort vector everywhere; `compute.WTCNA()` now takes one samples x TFs matrix plus the cohort vector (was a list of matrices with `batch = TRUE`) and returns the module scores in the input row order.
- [x] **`seed`:** only where something is random: `compute.latent_factors()` (NMF), `compute_factor_gsea()` (fgsea permutations) and `CellTFusion()` (new `seed`, default 123, passed to both and used for the deconvolution methods with random steps). TF and pathway scoring are deterministic (decoupleR fixes its own seed for `wsum`/`consensus`); the `geom_jitter()` of the `scores.*` and ANOVA boxplots only moves the plotted points.
- [x] **Threshold names:** `corr` now only means the deconvolution subgroup threshold; the module merging threshold is `corr_mod` everywhere (`mergeModules()*` too).
- [x] **Other:** `scores.stat.analysis(method)` renamed `test` (`method` always means the decoupleR statistic); `compute.survival.analysis(thres)` renamed `quantile_cutoff`; `compute.metadata.association()` `file_name` optional (default `NULL`, as everywhere) and `ncol_grid` default `NULL` (automatic, as in `compute.modules.relationship()`); `T`/`F` defaults written `TRUE`/`FALSE`.

Still open:

- [x] **Data inputs: one name per object** (done 15 October 2026): `raw.counts` -> `raw_counts`; log expression `RNA.tpm` / `RNA.counts` / `tpm.counts*` -> `counts_norm`; `TFs.matrix` -> `tfs`; `compute.WTCNA()` output `tfs.module.network` / `TF.network` -> `network`; deconvolution matrix `deconvolution` -> `deconv`; test deconvolution `deconvolution_test` / `test_deconv` -> `deconv_test`; deconvolution analysis output `deconv_res` / `deconvolution.subgroupped` -> `dt`; `cell.groups` -> `cell_groups`; dendrograms `cell.dendrograms` / `cell.group.dendrogram*` -> `cell_dendrograms`; `latent_spaces` -> `latent_factors` (`project_factors()`); score matrices `X` / `scores_test` / `features` / `features_df` / `Z` / `tfs.modules` -> `scores`; `survival.data` -> `coldata`; `PFS` / `PFS_event` -> `time_column` / `event_column`; `gsea_study` -> `gsea_results`; `features` -> `relationship` (`identify.cell.groups()`) and `group_names` (`compute.test.set()`); `groups` -> `features` (`extract_cells()`); `module_group` -> `module`; `hub_genes*` / `tfs_universe*` -> `hub_tfs` / `universe`. Output element names (e.g. `Latent_spaces`) and `matA` / `matB` unchanged.
- [x] **Parallel workers:** `prepare_celltfusion_folds()` has `do_parallel` (default FALSE) and `workers` (NULL = available cores minus one), with the same meaning as in `CellTFusion()`; `ncores` removed.
- [x] **Correlations not controlled by `corr_type`** (decided 15 October 2026): the CCA of `compute_composite_score()` and its orientation stay Pearson (user's decision: a rank-based CCA would make scores depend on the other samples, so single samples could not be projected); `identify_hub_TFs()` module membership now follows the network's `corr_type`.
- [x] (Closed 15 October 2026: kept as a known behaviour, user's decision.) **Batch-mode TF modules depend on which cohort comes first** (noted 15 October 2026): the first cohort is the calibration reference of the consensus TOM (as in WGCNA); with tutorial data, cohort B first gave 8 modules instead of 9. The user chose not to change it.
- [ ] **External callers to update** (not in this repository): MCP server (`mcp-tme-servers`) and paper scripts (`CellTFusion_paper/scripts`), see the lists given on 15 October 2026.

## Deferred: to do after the review of `CellTFusion()` (noted October 2026)

- [ ] **TCGA references.** On 16 October 2026 the Hallmark GSEA (`compute_factor_gsea()`), the TCGA meta-program mapping (`derive_meta_programs()`, `build_nes_matrix()`, `map_factors_to_metaprograms()`, `annotate_metaprograms_TME()`) and `inst/extdata/TCGA_meta_programs_*.RData` were removed (user: the mapping was too coarse). `inst/extdata/CellTFusion_TCGA_*.RData` are no longer used by any function: decide whether to delete them or rebuild them (they were built while cell groups were mislabelled, 2026-04-29 to 2026-10-05, and with older defaults).
- [x] **Shiny app.** Done 15 October 2026: new controls (scoring method, clustering method, BH adjustment, consensus quantile, seed, cancer type), defaults 5 / 5, the unused slider removed, the counts checkbox relabelled; the cell group table and the TME panel (latent factor scores labelled with meta-program and TME subtype) rewritten for the current results; TF activities of all cohorts shown in batch mode; the deconvolution PCA no longer fails on constant features. Tested with shiny::testServer() (outputs and a full run on the demo data). Original note: Launch and test it after the changes to `batch`, `normalized` and the removed arguments; remove the unused "High-correlation cutoff" slider; update the "Minimum TF targets" default (3, package default is 5) and the "Normalize counts (TPM + log2)" label (unticked now means the input is TPM).
- [x] **Vignettes.** Done 15 October 2026 (figures regenerated from CellTFusion() on the tutorial data with the current defaults; the meta-program mapping figure still uses the TCGA references to be rebuilt). Original note: Regenerate the figures of cell groups, latent factors, niches and TME states (made before the cell group fix).
- [x] **Data.** Done 15 October 2026: `traitdata.tuto` has the overall survival (`os`, `censOS`) of the same study, `tfs.tuto` and `network.tuto` were rebuilt with the current defaults (`data-raw/tutorial_data.R`); `deconv.tuto` did not need rebuilding (raw feature names are unchanged in multideconv 0.2.0; `deconv_subgroups.tuto` was already rebuilt). Original note: Regenerate `deconv.tuto` with multideconv >= 0.2.0 naming; add survival columns to the tutorial data, or a small survival example, so the survival examples can run. `tfs.tuto` and `network.tuto` were computed before the default scoring method of `compute.TFs.activity()` changed from the decoupleR consensus to `"ulm"` (11 October 2026) and before `compute.WTCNA()` built its TF tree on the TOM dissimilarity with average linkage (12 October 2026): regenerate them too.
- [ ] **Docs.** Review them again once the items above are done.
- [x] **Tests.** Done 15 October 2026: `tests/testthat/` with 21 test blocks (77 expectations, offline, about 5 minutes): TF modules (defaults, independence from TF scaling, batch cohort vector), cell groups and their projection (single cohort and batch), latent factors and niches, statistics (`padj`, `test`), survival on the tutorial overall survival, raw TF scores and batch pathway scores.
- [x] (Closed 15 October 2026, not pursued.) **`scores.*` tests always save their plots** (noted 15 October 2026): `scores.stat.analysis()` and the `scores.*` tests have no `return` argument and fail if `Results/` does not exist, unlike the other functions (`return = TRUE` by default).
- [x] **Scoring method in `CellTFusion()`.** Done 15 October 2026 (`method` argument). `compute.TFs.activity()` and `compute.pathway.activity()` have a `method` argument (default `"ulm"`), but `CellTFusion()` does not pass one, so the pipeline always uses `"ulm"`. Decide whether to add the argument to `CellTFusion()` (noted 11 October 2026).
- [x] Check behavior across all pipeline when batch is set to TRUE (noted 12 October 2026). Closed 15 October 2026: consensus modules, per-cohort TF and pathway scores, within-cohort z-scoring before the CCA, projection of new cohorts and the small-cohort warning were reviewed and fixed.
- [ ] **Revise `calculate_dendrogram_cuts()` to really understand how it works** (noted 12 October 2026): how `dynamicTreeCut::cutreeDynamic(method = "tree")` gives out cluster numbers (on the tutorial data, 55 of 123 cell groups were made of separated ranges of correlation), and whether the step added that splits clusters into single ranges and joins small pieces to their closest neighbour is the right answer.
- [ ] **Cell niches: more or fewer cells.** Since 15 October 2026 the factors use only the positive part of the cell group scores (high-activity side of the TF modules, no `_pos` / `_neg` split). A high-activity cell group can still mean more or fewer cells (depending on whether its features are positively or negatively related to its TF module), which `compute_cells_niches()` does not take into account. Earlier note: `compute_cells_niches()` adds together the selected cell group parts that mean more cells and those that mean fewer (it uses neither `_pos` / `_neg` nor the relation of each group to its TF module). Measured 12 October 2026 on the tutorial data (3 factors): the selected parts meaning more / fewer cells are 41 / 25, 29 / 40 and 38 / 17; the top cell type of Factor1, CAF (weight 0.14), gets 0.11 from parts meaning fewer CAFs and 0.03 from parts meaning more. Proposed: give each part its direction (from the relation of its group to the TF module; `Module_relation`, which gave it, was removed on 14 October 2026, so it would have to be computed here), compute weights and enrichment separately for more and fewer, return per factor a table (cell type, direction, weight, enrichment) and colour the plot edges by direction; this changes the structure of `Cells_niches`.

