# Multi-cohort analysis

``` r

library(CellTFusion)
#> Warning: replacing previous import 'S4Arrays::makeNindexFromArrayViewport' by
#> 'DelayedArray::makeNindexFromArrayViewport' when loading 'HDF5Array'
#> Warning: replacing previous import 'S4Arrays::makeNindexFromArrayViewport' by
#> 'DelayedArray::makeNindexFromArrayViewport' when loading 'SummarizedExperiment'
#> 
```

When samples come from several cohorts, batches or sequencing runs,
cohort-level differences can dominate TF activity and cell type
correlations and mask the biological signal. `CellTFusion` corrects for
this with the `batch` argument. This tutorial explains what `batch` does
at each step and how to use it.

## **What `batch` does at each step**

| Function | `batch` | Effect |
|----|----|----|
| `CellTFusion(batch = batch_vec)` | Cohort vector | Computes **TF activity separately per cohort** and builds a **consensus TF network**; the cohort vector is passed to every step below. |
| `compute.WTCNA(batch = batch_vec)` | Cohort vector | Takes the TF activities of all the samples (each cohort scored and scaled separately), splits them by cohort, picks a soft-thresholding power per cohort, combines the per-cohort topological overlaps into a consensus network (`consensus_quantile`, default 0: the minimum across cohorts) and builds the modules from it as in single-cohort mode. |
| `compute.deconvolution.analysis(batch = batch_vec)` | Cohort vector | Correlations used to build deconvolution subgroups control for cohort. |
| `compute.modules.relationship(batch = batch_vec)` | Cohort vector | TF module-feature correlations are **partial correlations** controlling for cohort (one indicator variable per cohort beyond the first), Pearson or Spearman according to `corr_type`. |
| `construct_cell_groups(batch = batch_vec)` | Cohort vector | Passes the cohort vector to the correlations and to the composite scores. Required for a network built with `compute.WTCNA(batch = batch_vec)`; it can be named by sample. |
| `compute_composite_score(batch = batch_vec)` | Cohort vector | Z-scores the cell group features within each cohort before the CCA with the TF module score, so composite scores are **corrected for cohort**. Only the CCA weights (`xcoef`) are kept: new cohorts are z-scored on their own means and SDs. |

In every function, including
[`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md),
`batch` is the cohort vector itself (one value per sample, or named by
sample), used to split the samples and to remove cohort effects. Cohort
is treated as categorical, and the samples do not need to be ordered by
cohort.

## **When to use it**

Use `batch` when samples do not come from a single homogeneous cohort —
for example, RNA-seq from different studies, sequencing batches or
clinical sites — and cohort identity could create TF co-activity or cell
type correlation patterns unrelated to biology. Without it, TF modules
and cell groups may reflect which cohort a sample comes from rather than
shared TME biology.

## **Implications**

- **Consensus modules are more conservative.** Only TF correlation
  structure shared by all cohorts is kept, so modules may be fewer or
  larger than with pooled samples.
- **Common TFs.** Only TFs present in every cohort are used to build the
  consensus network.
- **Correlation type.** `corr_type` (`"spearman"` by default, or
  `"pearson"`) is used in both single- and multi-cohort mode.
- **Partial correlation vs. correction.** Partial correlations only
  adjust a given correlation test, while the composite scores of the
  cell groups are themselves cohort-corrected, and so are the latent
  factors built from them.
- **Projection of new cohorts.** For models trained with `batch`, the
  features of new samples are z-scored on their own means and SDs (as
  each training cohort was), so each new cohort should be projected
  separately, with enough samples: a warning is given with fewer than 10
  (see [Machine learning
  workflows](https://verapancaldilab.github.io/CellTFusion/articles/a5_machine_learning.md)).

## **Example: multi-cohort pipeline**

``` r

# batch: cohort of each sample, in the same order as the columns of raw.counts
# (or named by sample, in which case it is matched to the samples by name)
res <- CellTFusion(
  raw_counts     = raw.counts,
  normalized     = TRUE,
  batch          = coldata$Cohort,
  deconv_methods = c("Quantiseq", "Epidish"),
  TF.collection  = "CollecTRI",
  corr           = 0.7,
  corr_mod       = 0.9,
  pval           = 0.05,
  file_name      = "Tutorial_batch",
  return         = TRUE
)
```

Internally, this:

1.  Splits the samples by cohort and runs
    [`compute.TFs.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.TFs.activity.md)
    (and
    [`compute.pathway.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.pathway.activity.md))
    on each cohort, then joins the TF activities of the TFs shared by
    all cohorts. TF and pathway scores are z-scored within each cohort
    (without `batch`,
    [`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)
    returns them unscaled).
2.  Runs `compute.WTCNA(batch = batch_vec)` on them to get consensus TF
    modules.
3.  Passes the cohort vector to `compute.deconvolution.analysis()` and
    [`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md),
    so subgroups, module-feature correlations and cell group scores are
    all cohort-corrected.

## **Example: running the steps individually**

When running the steps yourself (see [Feature
computation](https://verapancaldilab.github.io/CellTFusion/articles/a1_feature_computation.md)),
build the per-cohort TF list and the cohort vector:

``` r

batch_vec <- factor(coldata$Cohort)  # categorical
cohorts   <- split(seq_len(ncol(counts.norm)), batch_vec)

tfs_list <- lapply(cohorts, function(idx) {
  compute.TFs.activity(counts.norm[, idx, drop = FALSE], TF.collection = "CollecTRI")
})
common_tfs <- Reduce(intersect, lapply(tfs_list, colnames))
tfs <- do.call(rbind, unname(lapply(tfs_list, function(x) x[, common_tfs])))[colnames(counts.norm), ]

network <- compute.WTCNA(tfs = tfs, batch = batch_vec, min_module_size = 15)

dt <- multideconv::compute.deconvolution.analysis(deconv, corr = 0.7, batch = batch_vec)

cell_groups <- construct_cell_groups(network, dt, batch = batch_vec, pval = 0.05)
```

**Note:** the module matrix returned by
[`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md)
(`network[[1]]`) is in the row order of `tfs`, and
[`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md)
matches the samples by name.
