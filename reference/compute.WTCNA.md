# Compute Weighted TF-coactivity Network Analysis (WTCNA)

Construct a weighted signed or unsigned network using TF activity to
cluster protein regulators into modules that share similar activity
patterns. Each TF module will have a sample-level score represented by
the eigengene of the module.

## Usage

``` r
compute.WTCNA(
  tfs,
  batch = NULL,
  network_type = "signed",
  clustering_method = "ward.D2",
  min_module_size = 5,
  corr_mod = 0.9,
  corr_type = "spearman",
  verbose = TRUE,
  file_name = NULL,
  soft_power = NULL,
  consensus_quantile = 0,
  return = TRUE
)
```

## Arguments

- tfs:

  Matrix of TF activity (samples x TFs), scaled or not: the TFs are
  standardized (within each cohort with `batch`) to compute the module
  eigengenes, and the other steps use correlations, which do not depend
  on the scaling. With `batch`, the TF activities of all the cohorts,
  each computed within its cohort (as
  [`CellTFusion()`](https://verapancaldilab.github.io/CellTFusion/reference/CellTFusion.md)
  does).

- batch:

  (Optional) A vector with the cohort of each sample, in the row order
  of `tfs`; if it is named, it is matched to the samples by name
  instead. It is always treated as categorical. If given, consensus TF
  modules are built across the cohorts (see Details) and module
  eigengenes are computed and scaled within each cohort. Default is
  `NULL` (single cohort); a vector with only one cohort is treated the
  same way.

- network_type:

  Network type: "signed" (default) or "unsigned".

- clustering_method:

  Linkage method of the hierarchical clustering of the TFs
  ([`stats::hclust()`](https://rdrr.io/r/stats/hclust.html)), applied to
  the TOM dissimilarity (1 - topological overlap) between TFs, or to the
  consensus TOM dissimilarity in batch mode. Default is "ward.D2"
  (WGCNA's own default is "average": two groups of TFs are joined
  according to the average dissimilarity between their TFs).

- min_module_size:

  Minimum number of TFs per module. Default is 5.

- corr_mod:

  Correlation threshold (0-1): modules whose eigengenes have a
  correlation (of type `corr_type`) above it are merged. In batch mode,
  the correlation is computed on the eigengenes of all the samples,
  scaled within each cohort. Default is 0.9.

- corr_type:

  Correlation used to pick the soft-threshold power, build the adjacency
  matrix, merge modules and compute the proportion of variance
  explained: "spearman" (default) or "pearson", in single-cohort and
  batch mode. The module eigengene is always the first principal
  component of the module's TF activities.

- verbose:

  Logical; whether to print the function messages. Default is TRUE.

- file_name:

  Optional character suffix used when writing WTCNA outputs (CSV files
  and, in single-cohort mode, the diagnostic plots
  `Soft_Threshold_<file_name>.pdf`,
  `Gene_dendrogram_and_module_colors_<file_name>.pdf` and
  `Gene_dendrogram_and_module_colors_after_merging_<file_name>.pdf`).

- soft_power:

  Optional numeric value specifying the soft-thresholding power used to
  build the adjacency matrix (one value per cohort, in the order the
  cohorts first appear in `batch`, in batch mode). If `NULL`, the lowest
  power whose scale-free fit \\R^2\\ reaches 0.9 is chosen automatically
  (the power with the best fit if none reaches 0.9, with a message).

- consensus_quantile:

  Numeric between 0 and 1 (batch mode only). For each pair of TFs, the
  consensus topological overlap is this quantile of its (calibrated)
  topological overlaps across the cohorts. With 0 (default), the
  minimum: two TFs are close only if they are close in every cohort, the
  strictest choice, which removes more TF relationships as more cohorts
  are added. Higher values are less strict: with `K` cohorts,
  `consensus_quantile = k / (K - 1)` takes the `(k + 1)`-th lowest
  value, so two TFs are close if they are close in all cohorts but `k`
  (e.g. 0.25 with 5 cohorts: in at least 4 of them); 0.5 is the median.

- return:

  Logical, whether to save output plots, the module list
  (`Results/TFs_modules_<file_name>.csv`) and the module matrix
  (`Results/TF_module_matrix_<file_name>.csv`, with the sample names as
  row names) to "Results/" (without `file_name`, `TFs_modules.csv` and
  `TF_module_matrix.csv`). Default is TRUE.

## Value

A named list with:

- `TFs module matrix`: Scaled module eigengenes (samples x modules), in
  the row order of `tfs`.

- `TFs colors`: Vector of module colors assigned to each TF (`"grey"`
  for TFs not assigned to any module).

- `TFs per module`: List of TF names in each module.

- `Proportion of variance`: Variance explained per module (single-cohort
  mode only).

- `TFs_matrix`: The TF activity matrix used (a list of per-cohort
  matrices in batch mode, in the order the cohorts first appear in
  `batch`).

- `Network parameters`: A list with the soft-thresholding power(s) used
  (`soft_power`, one per cohort in batch mode), `network_type` and
  `corr_type`, used by
  [`identify_hub_TFs()`](https://verapancaldilab.github.io/CellTFusion/reference/identify_hub_TFs.md),
  and in batch mode `consensus_quantile`.

## Details

TFs that WGCNA does not assign to any module (color `"grey"`) are not a
module: they are left out of the module matrix and of `TFs per module`,
with a message.

In batch mode, the modules are built with the same steps as in
single-cohort mode, on a consensus network:

1.  the samples are split by cohort, and a soft-thresholding power is
    chosen for each cohort;

2.  the adjacency and the topological overlap matrix (TOM) are computed
    in each cohort, with its own power;

3.  the TOMs are calibrated to a common scale as in
    [`WGCNA::blockwiseConsensusModules()`](https://rdrr.io/pkg/WGCNA/man/blockwiseConsensusModules.html):
    each TOM is raised to the power that makes its 0.95 quantile equal
    to that of the first cohort;

4.  the consensus TOM is, for each pair of TFs, the `consensus_quantile`
    quantile of its calibrated TOMs across cohorts;

5.  the TFs are clustered on 1 - consensus TOM with `clustering_method`
    and the tree is cut with the same dynamic tree cut as in
    single-cohort mode;

6.  the module eigengenes are computed and scaled within each cohort,
    and the modules are merged as in single-cohort mode (see
    `corr_mod`); merged scores are scaled within each cohort again.

Unlike
[`WGCNA::blockwiseConsensusModules()`](https://rdrr.io/pkg/WGCNA/man/blockwiseConsensusModules.html),
the TFs of each module are not filtered or reassigned afterwards by
their module membership (as in single-cohort mode). The diagnostic plots
are only saved in single-cohort mode.

## References

Langfelder, P., & Horvath, S. (2008). WGCNA: an R package for weighted
correlation network analysis. BMC Bioinformatics, 9, 559.
https://doi.org/10.1186/1471-2105-9-559

## Examples

``` r

data("tfs.tuto")
network <- compute.WTCNA(tfs.tuto, corr_mod = 0.9, return = FALSE)
#> Creating weighted TF-coactivity network......................................................................
#> 
#> Warning: executing %dopar% sequentially: no parallel backend registered
#> Automatically choosing 15 as soft-threshold......................................................................
#> 
#> Calculating nodes adjacency and topological overlapping nodes.................................................
#> 
#> Calculating eigenvectors from modules.................................................
#> 
#> [1] "Merging modules significantly correlated with 0.9........"
```
