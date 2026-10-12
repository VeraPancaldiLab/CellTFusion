# Derive TME meta-programs by clustering Hallmarks across NMF factors

Hierarchically clusters (Ward.D2 on Euclidean distance) Hallmark gene
sets by their NES profile across NMF factors to identify recurrent
transcriptional programs in the TME. Use
[`annotate_metaprograms_TME()`](https://verapancaldilab.github.io/CellTFusion/reference/annotate_metaprograms_TME.md)
afterwards to label each meta-program with a Bagaev et al. (2021) MFP
subtype. The meta-program numbers (`MP1`, `MP2`, ...) follow the order
of the clusters and have no meaning of their own.

## Usage

``` r
derive_meta_programs(gsea_results, k = NULL, file_name = NULL, plot = TRUE)
```

## Arguments

- gsea_results:

  Output of
  [`compute_factor_gsea()`](https://verapancaldilab.github.io/CellTFusion/reference/compute_factor_gsea.md)
  (a list with a `GSEA_results` element, one table per NMF factor). To
  derive meta-programs from several cohorts, put the `GSEA_results`
  lists of the runs together first, with distinct factor names, e.g.
  `list(GSEA_results = c(setNames(gsea_A$GSEA_results, paste0("A_", names(gsea_A$GSEA_results))), setNames(gsea_B$GSEA_results, paste0("B_", names(gsea_B$GSEA_results)))))`.

- k:

  Integer. Number of meta-programs (clusters) to extract. If `NULL`
  (default), chosen at the elbow of the within-cluster sum of squares
  (the k with the largest second difference) for k = 2 to min(10,
  n_hallmarks - 1).

- file_name:

  Optional character suffix for the saved heatmap
  (`Results/TCGA_meta_programs_<file_name>.pdf`, or
  `Results/TCGA_meta_programs.pdf` without it).

- plot:

  Logical. If `TRUE` (default), saves a clustering heatmap of the NES
  (Hallmarks x factors); NES beyond -3 and 3 are shown with the extreme
  colours. With a single factor, the factor columns are not clustered.

## Value

A data frame with one row per meta-program and columns `meta_program`
(`"MP1"`, `"MP2"`, ...) and `hallmarks` (comma-separated Hallmark gene
set names).

## Examples

``` r
if (FALSE) { # \dontrun{
gsea <- compute_factor_gsea(RNA.tpm = counts.norm, features_df = res$Latent_spaces$Z, plot_dot = FALSE)
meta_programs <- derive_meta_programs(gsea, file_name = "my_cohort")
} # }
```
