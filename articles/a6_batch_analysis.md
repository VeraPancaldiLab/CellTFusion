# Multi-cohort analysis

``` r

library(CellTFusion)
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
| `CellTFusion(batch = TRUE, batch_id = "Cohort")` | Logical + column of `coldata` | Computes **TF activity separately per cohort** and builds a **consensus TF network**; the cohort vector `coldata[, batch_id]` is passed to every step below. |
| `compute.WTCNA(batch = TRUE)` | Logical | Takes a **list** of per-cohort TF matrices, picks a soft-thresholding power per cohort and builds consensus modules with [`WGCNA::blockwiseConsensusModules()`](https://rdrr.io/pkg/WGCNA/man/blockwiseConsensusModules.html). |
| `compute.deconvolution.analysis(batch = batch_vec)` | Cohort vector | Correlations used to build deconvolution subgroups control for cohort. |
| `compute.modules.relationship(batch = batch_vec)` | Cohort vector | TF module-feature correlations are **partial correlations** controlling for cohort (one indicator variable per cohort). |
| `construct_cell_groups(batch = batch_vec)` | Cohort vector | Passes the cohort vector to the correlations and to the composite scores. |
| `compute_composite_score(batch = batch_vec)` | Cohort vector | Removes each cohort’s mean from the cell group features and TF activities before the CCA, so composite scores are **corrected for cohort**. |

In short, in
[`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)
and
[`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md),
`batch` is a switch that changes how the TF network is built. In every
later step, `batch` is the cohort vector itself, used to remove cohort
effects. Cohort is treated as categorical, and the samples do not need
to be ordered by cohort.

## **When to use it**

Use `batch = TRUE` when samples do not come from a single homogeneous
cohort — for example, RNA-seq from different studies, sequencing batches
or clinical sites — and cohort identity could create TF co-activity or
cell type correlation patterns unrelated to biology. Without it, TF
modules and cell groups may reflect which cohort a sample comes from
rather than shared TME biology.

## **Implications**

- **Consensus modules are more conservative.** Only TF correlation
  structure shared by all cohorts is kept, so modules may be fewer or
  larger than with pooled samples.
- **Common TFs.** Only TFs present in every cohort are used to build the
  consensus network.
- **Pearson only.** Consensus modules support Pearson correlation;
  `cor_type = "s"` is only available for a single cohort.
- **Partial correlation vs. correction.** Partial correlations only
  adjust a given correlation test, while the composite scores of the
  cell groups are themselves cohort-corrected, and so are the latent
  factors built from them.
- **Projection of new cohorts.** For models trained with `batch = TRUE`,
  new samples are centered on their own means (as each training cohort
  was), so each new cohort should be projected separately (see [Machine
  learning
  workflows](https://verapancaldilab.github.io/CellTFusion/articles/a5_machine_learning.md)).

## **Example: multi-cohort pipeline**

``` r

# coldata must contain a column with the cohort of each sample, e.g. "Cohort"
res <- CellTFusion(
  raw.counts     = raw.counts,
  normalized     = TRUE,
  coldata        = coldata,
  batch          = TRUE,
  batch_id       = "Cohort",
  deconv_methods = c("Quantiseq", "Epidish"),
  TF.collection  = "CollecTRI",
  cancer_type    = "skcm",
  corr           = 0.7,
  corr_mod       = 0.9,
  pval           = 0.05,
  file_name      = "Tutorial_batch",
  return         = TRUE
)
```

Internally, this:

1.  Splits the samples by `coldata$Cohort` and runs
    [`compute.TFs.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.TFs.activity.md)
    on each cohort, giving a list of TF activity matrices.
2.  Runs `compute.WTCNA(batch = TRUE)` on that list to get consensus TF
    modules.
3.  Passes `coldata[, "Cohort"]` to `compute.deconvolution.analysis()`
    and
    [`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md),
    so subgroups, module-feature correlations and cell group scores are
    all cohort-corrected.

## **Example: running the steps individually**

When running the steps yourself (see [Feature
computation](https://verapancaldilab.github.io/CellTFusion/articles/a1_feature_computation.md)),
build the per-cohort TF list and the cohort vector:

``` r

batch_vec <- coldata[, "Cohort"]
cohorts   <- split(seq_len(ncol(counts.norm)), batch_vec)

tfs_list <- lapply(cohorts, function(idx) {
  compute.TFs.activity(counts.norm[, idx, drop = FALSE], TF.collection = "CollecTRI")
})

network <- compute.WTCNA(TFs.matrix = tfs_list, batch = TRUE, minMod = 15)

dt <- multideconv::compute.deconvolution.analysis(deconv, corr = 0.7, batch = batch_vec)

cell_groups <- construct_cell_groups(network, dt, batch = batch_vec, pval = 0.05)
```
