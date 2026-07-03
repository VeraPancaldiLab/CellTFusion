# Batch/Multi-cohort Analysis

``` r

library(CellTFusion)
#> 
#> 
```

When samples come from multiple cohorts, batches, or sequencing runs,
technical and cohort-level variation can dominate TF activity and
cell-type correlation structure, masking real biological signal.
CellTFusion handles this in two distinct — and complementary — ways
depending on the pipeline stage, both controlled by the `batch`
argument. This article explains what `batch` actually does at each
stage, why it matters, and how to use it.

## Two different meanings of `batch`

`batch` is not a single mechanism reused everywhere — its type and
effect differ by function:

| Function | `batch` type | What it does |
|----|----|----|
| `CellTFusion(batch = TRUE, batch_id = "Cohort")` | Logical + column name | Splits samples by cohort and computes **TF activity separately per cohort**; runs **consensus WGCNA** instead of a single network. |
| `compute.WTCNA(batch = TRUE)` | Logical | Same as above at the network-construction level: expects `TFs.matrix` to be a **list** of per-cohort matrices and calls [`WGCNA::blockwiseConsensusModules()`](https://rdrr.io/pkg/WGCNA/man/blockwiseConsensusModules.html). |
| `compute.deconvolution.analysis(batch = batch_vec)` | Vector | Cell-type correlations used for subgrouping are computed as **partial correlations controlling for cohort**. |
| `compute.modules.relationship(batch = batch_vec)` | Vector | TF module vs. deconvolution (or pathway) correlations are computed as **partial correlations** via `ppcor::pcor.test(x, y, batch)`, removing the linear effect of cohort before testing association. |
| `construct_cell_groups(batch = batch_vec)` | Vector | Pass-through: forwards the batch vector to [`compute.modules.relationship()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.modules.relationship.md) and [`cell.groups.computation()`](https://verapancaldilab.github.io/CellTFusion/reference/cell.groups.computation.md) (which forwards it to [`compute_composite_score()`](https://verapancaldilab.github.io/CellTFusion/reference/compute_composite_score.md)). |
| [`compute_composite_score()`](https://verapancaldilab.github.io/CellTFusion/reference/compute_composite_score.md) (internal, called by [`cell.groups.computation()`](https://verapancaldilab.github.io/CellTFusion/reference/cell.groups.computation.md)) | Vector | Cell-group and TF-module scores are **residualized against batch** via `residuals(lm(x ~ batch))` before being used, removing any linear cohort effect from the composite scores themselves. |

In short: at the **TF network stage** (`CellTFusion`/`compute.WTCNA`),
`batch` is a **logical switch** that changes the algorithm (per-cohort
computation + consensus WGCNA). At every downstream stage (deconvolution
analysis, module-relationship correlations, cell group construction,
composite scores), `batch` is the **actual cohort vector**, used as a
covariate to regress out or partial out cohort effects.

## When to use it

Use `batch = TRUE` whenever your samples are not a single homogeneous
cohort — for example, RNAseq from different studies, sequencing batches,
or clinical sites — and cohort identity could plausibly create spurious
TF co-activity or cell-type correlation patterns unrelated to biology.

If you skip `batch` on multi-cohort data, TF modules and cell groups
risk being driven primarily by which cohort a sample came from, rather
than by shared TME biology.

## Implications

- **Consensus WGCNA is more conservative.** It only keeps TF-TF
  correlation structure that is consistent across all cohorts, so
  `batch = TRUE` can produce fewer, larger, or differently colored
  modules than pooling all samples into one matrix.
- **Common TF set.** When cohorts don’t share every TF (e.g. different
  regulon coverage), `compute.WTCNA(batch = TRUE)` restricts analysis to
  the TFs present in all cohorts before building the consensus network.
- **Partial correlation vs. residualization.** Partial correlations
  (`compute.deconvolution.analysis()`,
  [`compute.modules.relationship()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.modules.relationship.md))
  remove the *linear association with batch* only for the purpose of
  that specific correlation test; residualization
  ([`compute_composite_score()`](https://verapancaldilab.github.io/CellTFusion/reference/compute_composite_score.md))
  actually replaces the composite score values with their batch-adjusted
  residuals, which propagates downstream (e.g. into cell group scores
  and, from there, into latent factors).
- **`batch_id` must be provided consistently.**
  [`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)
  requires both `coldata` and `batch_id` when `batch = TRUE`; the same
  cohort vector (`coldata[, batch_id]`) is threaded through every
  downstream function automatically — you do not need to pass it
  manually if you use the
  [`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)
  wrapper.

## Example: multi-cohort pipeline

``` r

raw.counts <- CellTFusion::raw.counts.tuto
traitdata  <- CellTFusion::traitdata.tuto
# traitdata must contain a column identifying the cohort of each sample, e.g. "Cohort"

res <- CellTFusion(
  raw.counts     = raw.counts,
  normalized     = FALSE,
  coldata        = traitdata,
  task           = "unsupervised",
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

Internally this:

1.  Splits samples by `traitdata$Cohort` and computes
    [`compute.TFs.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.TFs.activity.md)
    separately for each cohort, yielding a **list** of TF activity
    matrices.
2.  Runs `compute.WTCNA(batch = TRUE)` on that list, producing a
    consensus TF module network shared across cohorts.
3.  Passes `batch_vec <- traitdata[, "Cohort"]` into
    `compute.deconvolution.analysis()`,
    [`compute.modules.relationship()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.modules.relationship.md),
    and
    [`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md),
    so subgroup correlations, module-deconvolution associations, and
    composite cell-group scores are all cohort-adjusted.

## Example: using individual functions directly

If you are running the feature computation steps manually (see [Feature
Computation](https://verapancaldilab.github.io/CellTFusion/articles/01-feature-computation.md))
rather than through
[`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md),
you must build the per-cohort TF list and batch vector yourself:

``` r

batch_vec <- traitdata[, "Cohort"]
cohorts   <- split(seq_len(ncol(counts.norm)), batch_vec)

tfs_list <- lapply(cohorts, function(idx) {
  compute.TFs.activity(counts.norm[, idx, drop = FALSE], TF.collection = "CollecTRI")
})

network <- compute.WTCNA(TFs.matrix = tfs_list, batch = TRUE, minMod = 15, corr_mod = 0.9)

dt <- multideconv::compute.deconvolution.analysis(deconv, corr = 0.7, batch = batch_vec)

cell_groups <- construct_cell_groups(network, dt, batch = batch_vec, pval = 0.05)
```
