# Statistical analysis

``` r

library(CellTFusion)
#> Warning: replacing previous import 'S4Arrays::makeNindexFromArrayViewport' by
#> 'DelayedArray::makeNindexFromArrayViewport' when loading 'HDF5Array'
#> Warning: replacing previous import 'S4Arrays::makeNindexFromArrayViewport' by
#> 'DelayedArray::makeNindexFromArrayViewport' when loading 'SummarizedExperiment'
#> 
```

This tutorial shows how to test associations between `CellTFusion`
features — latent factors, TF modules or cell groups — and clinical
variables. It uses the result of
[`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md):

``` r

raw.counts <- CellTFusion::raw.counts.tuto
traitdata  <- CellTFusion::traitdata.tuto

res <- CellTFusion(
  raw_counts  = raw.counts,
  normalized  = TRUE,
  return      = TRUE
)
```

## **Latent factors and clinical traits**

[`scores.stat.analysis()`](https://verapancaldilab.github.io/CellTFusion/reference/scores.stat.analysis.md)
tests each latent factor against a clinical variable and saves a plot
for every significant association in `Results/`. `scores` can be the
output of
[`compute.latent_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.latent_factors.md)
(`res$Latent_spaces`), a samples x features matrix, or a list whose
first element is that matrix (e.g. `res$Cell_groups`). The samples are
matched to the rows of `coldata` by row name (without row names in
`coldata`, its rows must be in the order of the scores), and samples
with a missing trait are left out. The trait must be categorical.

Available tests (`test`):

- `"wilcox"` — Wilcoxon rank-sum test (two groups).
- `"ttest"` — Welch’s t-test (two groups, unequal variances).
- `"kruskal"` — Kruskal-Wallis test with Dunn post-hoc tests (two or
  more groups).
- `"anova"` — one-way ANOVA with Tukey post-hoc tests (two or more
  groups).
- `"fisher"` — Fisher’s exact test, with scores split at the median into
  High/Low.

``` r

sig_factors <- scores.stat.analysis(
  scores  = res$Latent_spaces,
  coldata = traitdata,
  trait   = "Best.Confirmed.Overall.Response",
  test    = "wilcox",
  pval    = 0.05
)
```

The result keeps only the significant features (its first element has no
column if none is significant), and its element `Statistics` gives the
test of every feature (statistic, p-value, adjusted p-value). Use
`padj = TRUE` to correct the p-values for the number of features tested
(Benjamini-Hochberg); it is `FALSE` by default. The individual tests
([`scores.wilcox.test()`](https://verapancaldilab.github.io/CellTFusion/reference/scores.wilcox.test.md),
[`scores.ttest()`](https://verapancaldilab.github.io/CellTFusion/reference/scores.ttest.md),
[`scores.kruskal.test()`](https://verapancaldilab.github.io/CellTFusion/reference/scores.kruskal.test.md),
[`scores.anova.test()`](https://verapancaldilab.github.io/CellTFusion/reference/scores.anova.test.md),
[`scores.fisher.test()`](https://verapancaldilab.github.io/CellTFusion/reference/scores.fisher.test.md))
can also be called directly with a list whose first element is the score
matrix.

The factor scores themselves are in `res$Latent_spaces$Z`:

``` r

factor_scores <- data.frame(res$Latent_spaces$Z)
head(factor_scores)
```

## **TF modules and clinical traits**

[`compute.metadata.association()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.metadata.association.md)
relates TF module scores to all clinical variables at once. The rows of
`coldata` must be the same samples, in the same order, as the rows of
the module matrix (an error is raised otherwise):

- **Numeric traits** are correlated with each module (Spearman by
  default, `corr_type = "pearson"` for Pearson) and shown in a heatmap,
  annotated for p-values at or below `pval`. For traits with missing
  values, each p-value uses the number of samples with both values.
  Numeric columns are treated as continuous: convert number-coded
  categories (e.g. a stage coded 1 to 4) to factor to test them as
  groups.
- **Categorical traits** (text, factor and TRUE/FALSE columns) are
  tested with one-way ANOVA and Tukey post-hoc tests when
  `plot_grid = TRUE`, and saved as boxplot grids with one panel per
  significant module.

Traits that cannot be tested are skipped with a message: numeric traits
without variation, and categorical traits with a single level or with
one value per sample (e.g. a sample identifier). Use `padj = TRUE` to
correct the p-values for multiple testing across the modules, for each
trait (Benjamini-Hochberg).

``` r

assoc <- compute.metadata.association(
  scores = res$TF_network[[1]],
  coldata     = traitdata,
  pval        = 0.05,
  file_name   = "Tutorial",
  plot_grid   = TRUE
)

assoc$correlation  # module-trait correlations (numeric traits)
assoc$p_value      # their p-values
assoc$anova        # ANOVA of every module for every categorical trait
```

## **TF modules and pathways**

[`compute.modules.relationship()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.modules.relationship.md)
correlates two matrices of continuous features over the same samples,
for example TF modules and pathway activities (Spearman by default,
`corr_type = "pearson"` for Pearson), and saves a heatmap of the
significant correlations. It is for numeric variables only; use
[`compute.metadata.association()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.metadata.association.md)
for categorical ones. With `plot_grid = TRUE`, a scatter plot is also
drawn for the significant pairs (the 100 strongest at most).

``` r

rel <- compute.modules.relationship(
  matA      = res$TF_network[[1]],
  matB      = res$Pathways_scores,
  file_name = "Pathways_vs_TF_modules"
)

rel$correlation  # correlations (TF modules x pathways)
rel$p_value      # their p-values
rel$significant  # for each TF module, the pathways significantly correlated with it
```

## **Survival analysis**

[`compute.survival.analysis()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.survival.analysis.md)
compares Kaplan-Meier survival curves with a log-rank test (packages
`survival`, `survminer` and `gridExtra`). It works in two modes:

1.  **Predefined groups** (`group_column`) — e.g. a clinical risk group
    or a cluster assignment in the clinical data.
2.  **Feature screening** (`scores`) — e.g. `res$Latent_spaces$Z`. Each
    feature is split into High/Low at a quantile (`quantile_cutoff`),
    and the features with a significant log-rank test (`pval`,
    optionally corrected with `padj = TRUE`) are returned, with a
    `Statistics` table giving the test of every feature. Samples without
    a value for `group_column` are left out. The rows of `features` must
    be the same samples, in the same order, as the rows of the clinical
    data.

Both modes need a time-to-event column and an event column (1 = event, 0
= censored) in the clinical data. The tutorial `traitdata` includes the
overall survival of the patients: `os` (months) and `censOS` (1 = death,
0 = censored).

``` r

# Predefined groups
surv_group <- compute.survival.analysis(
  coldata       = traitdata,
  time_column   = "os",
  event_column  = "censOS",
  group_column  = "Best.Confirmed.Overall.Response",
  file_name     = "Tutorial"
)

surv_group$p_value     # log-rank p-value
surv_group$median_PFS  # median survival per group (column name kept from the function output)
```

``` r

# Feature screening on latent factors
surv_factors <- compute.survival.analysis(
  coldata       = traitdata,
  time_column   = "os",
  event_column  = "censOS",
  scores        = res$Latent_spaces$Z,
  pval          = 0.05,
  quantile_cutoff = 0.5,   # median split; e.g. 0.75 for top quartile vs rest
  file_name     = "Tutorial"
)
```

A Kaplan-Meier plot with a risk table is saved (always for predefined
groups, and for every significant feature in feature screening) to
`Results/SurvPlot_<group-or-feature>_<file_name>.svg`
(`Results/SurvPlot_<group-or-feature>.svg` without `file_name`).
