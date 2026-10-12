# Compute Weighted TF-coactivity Network Analysis (WTCNA)

Construct a weighted signed or unsigned network using TF activity to
cluster protein regulators into modules that share similar activity
patterns. Each TF module will have a sample-level score represented by
the eigengene of the module.

## Usage

``` r
compute.WTCNA(
  TFs.matrix,
  batch = FALSE,
  network.type = "signed",
  clustering.method = "average",
  minMod = 15,
  corr_mod = 0.9,
  cor_type = "p",
  verbose = F,
  file.name = NULL,
  softPower = NULL,
  return = T
)
```

## Arguments

- TFs.matrix:

  Matrix of TF activity (samples x TFs). When `batch = TRUE`, a list of
  such matrices, one per cohort.

- batch:

  Logical; if TRUE, performs consensus WGCNA
  ([`WGCNA::blockwiseConsensusModules()`](https://rdrr.io/pkg/WGCNA/man/blockwiseConsensusModules.html))
  across cohorts provided as a list of matrices. In this mode
  `clustering.method` is not used (WGCNA builds the consensus tree
  itself), modules are merged by WGCNA with
  `mergeCutHeight = 1 - corr_mod`, and module eigengenes are scaled
  within each cohort.

- network.type:

  Network type: "signed" (default) or "unsigned".

- clustering.method:

  Linkage method of the hierarchical clustering of the TFs
  ([`stats::hclust()`](https://rdrr.io/r/stats/hclust.html)), applied to
  the TOM dissimilarity (1 - topological overlap) between TFs
  (single-cohort mode only). Default is "average", as in WGCNA: two
  groups of TFs are joined according to the average dissimilarity
  between their TFs.

- minMod:

  Minimum number of TFs per module. Default is 15.

- corr_mod:

  Correlation threshold (0-1): modules whose eigengenes have a
  correlation (of type `cor_type`) above it are merged. Default is 0.9.

- cor_type:

  Correlation used to pick the soft-threshold power, build the adjacency
  matrix, merge modules and compute the proportion of variance
  explained: "p" (Pearson) or "s" (Spearman). The module eigengene is
  always the first principal component of the module's TF activities.
  Spearman is only available when `batch = FALSE`, because WGCNA
  consensus modules only support Pearson correlation. Default is "p".

- verbose:

  Logical; whether to print the function messages. Default is FALSE.

- file.name:

  Optional character suffix used when writing WTCNA outputs (CSV files
  and, in single-cohort mode, the diagnostic plots
  `Soft_Threshold_<file.name>.pdf`,
  `Gene_dendrogram_and_module_colors_<file.name>.pdf` and
  `Gene_dendrogram_and_module_colors_after_merging_<file.name>.pdf`).

- softPower:

  Optional numeric value specifying the soft-thresholding power used to
  build the adjacency matrix (one value per cohort when `batch = TRUE`).
  If `NULL`, the lowest power whose scale-free fit \\R^2\\ reaches 0.9
  is chosen automatically (the power with the best fit if none reaches
  0.9, with a message).

- return:

  Logical, whether to save output plots, the module list
  (`Results/TFs_modules_<file.name>.csv`) and the module matrix
  (`Results/TF_module_matrix_<file.name>.csv`, with the sample names as
  row names) to "Results/" (without `file.name`, `TFs_modules.csv` and
  `TF_module_matrix.csv`). Default is TRUE.

## Value

A named list with:

- `TFs module matrix`: Scaled module eigengenes (samples x modules). In
  batch mode, samples are concatenated in cohort order (use the row
  names to match them to other tables).

- `TFs colors`: Vector of module colors assigned to each TF (`"grey"`
  for TFs not assigned to any module).

- `TFs per module`: List of TF names in each module.

- `Proportion of variance`: Variance explained per module (single-cohort
  mode only).

- `TFs_matrix`: The TF activity matrix used (a list of per-cohort
  matrices restricted to shared TFs in batch mode).

- `Network parameters`: A list with the soft-thresholding power(s) used
  (`softPower`, one per cohort in batch mode), `network.type` and
  `cor_type`, used by
  [`identify_hub_TFs()`](https://verapancaldilab.github.io/CellTFusion/reference/identify_hub_TFs.md).

## Details

TFs that WGCNA does not assign to any module (color `"grey"`) are not a
module: they are left out of the module matrix and of `TFs per module`,
with a message.

In batch mode, the soft-thresholding powers are chosen on the TFs shared
by all cohorts, the ones used to build the network.

In batch mode,
[`WGCNA::blockwiseConsensusModules()`](https://rdrr.io/pkg/WGCNA/man/blockwiseConsensusModules.html)
can fail with `"dim(X) must have a positive length"` when a consensus
module has exactly one TF to reassign to another module (a known WGCNA
issue). The problem goes away when the modules change. If the powers
were chosen automatically (`softPower = NULL`), `compute.WTCNA()`
retries moving the power of one cohort one step down or up, and gives a
warning with the powers that were finally used. If `softPower` was
given, or no neighbouring power works, it stops with a message reporting
the powers: set `softPower` to slightly different values or change
`minMod` (a small change may not be enough).

## References

Langfelder, P., & Horvath, S. (2008). WGCNA: an R package for weighted
correlation network analysis. BMC Bioinformatics, 9, 559.
https://doi.org/10.1186/1471-2105-9-559

## Examples

``` r

data("tfs.tuto")
network <- compute.WTCNA(tfs.tuto, corr_mod = 0.9, return = FALSE)
#> Warning: executing %dopar% sequentially: no parallel backend registered
```
