# Compute Weighted TF-coactivity Network Analysis (WTCNA)

Construct a weighted signed or unsigned network using TF activity to
cluster protein regulators into modules that share similar activity
patterns. Each TF module will have a sample-level score represented by
the eigenvalue of the module.

## Usage

``` r
compute.WTCNA(
  TFs.matrix,
  batch = FALSE,
  network.type = "signed",
  clustering.method = "ward.D2",
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

  Matrix of TF activity (samples x TFs).

- batch:

  Logical; if TRUE, performs consensus WGCNA
  ([`WGCNA::blockwiseConsensusModules()`](https://rdrr.io/pkg/WGCNA/man/blockwiseConsensusModules.html))
  across cohorts provided as a list of matrices. In this mode
  `clustering.method` and `corr_mod` are not used (modules are merged
  with a fixed `mergeCutHeight = 0.25`), and module eigengenes are
  scaled within each cohort.

- network.type:

  Network type: "signed", "unsigned", "signed hybrid", or "distance".
  Default is "signed".

- clustering.method:

  Clustering method for hierarchical clustering (single-cohort mode
  only). Default is "ward.D2".

- minMod:

  Minimum number of TFs per module. Default is 15.

- corr_mod:

  Correlation threshold (0-1) above which module eigengenes are merged
  (single-cohort mode only). Default is 0.9.

- cor_type:

  Correlation used to pick the soft-threshold power and build the
  adjacency matrix: "p" (Pearson) or "s" (Spearman). Spearman is only
  available when `batch = FALSE`, because WGCNA consensus modules only
  support Pearson correlation. Default is "p".

- verbose:

  Boolen value to whether print or no the function messages

- file.name:

  Optional character suffix used when writing WTCNA outputs.

- softPower:

  Optional numeric value specifying the soft-thresholding power used to
  build the adjacency matrix (one value per cohort when `batch = TRUE`).
  If `NULL`, the power whose scale-free fit \\R^2\\ is closest to 0.9 is
  chosen automatically.

- return:

  Logical, whether to save output plots and module list to "Results/".
  Default is TRUE.

## Value

A named list with:

- `TFs module matrix`: Scaled module eigengenes (samples x modules). In
  batch mode, samples are concatenated in cohort order.

- `TFs colors`: Vector of module colors assigned to each TF.

- `TFs per module`: List of TF names in each module.

- `Proportion of variance`: Variance explained per module (single-cohort
  mode only).

- `TFs_matrix`: The TF activity matrix used (a list of per-cohort
  matrices restricted to shared TFs in batch mode).

## References

Langfelder, P., & Horvath, S. (2008). WGCNA: an R package for weighted
correlation network analysis. BMC Bioinformatics, 9, 559.
https://doi.org/10.1186/1471-2105-9-559

## Examples

``` r

data("tfs.tuto")
network <- compute.WTCNA(tfs.tuto, corr_mod = 0.9, clustering.method = "ward.D2", return = FALSE)
#> Warning: executing %dopar% sequentially: no parallel backend registered
```
