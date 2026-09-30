# Derive TME meta-programs by clustering Hallmarks across NMF factors

Hierarchically clusters (Ward.D2 on Euclidean distance) Hallmark gene
sets by their NES profile across NMF factors to identify recurrent
transcriptional programs in the TME. Use
[`annotate_metaprograms_TME()`](https://verapancaldilab.github.io/CellTFusion/reference/annotate_metaprograms_TME.md)
afterwards to label each meta-program with a Bagaev et al. (2021) MFP
subtype.

## Usage

``` r
derive_meta_programs(gsea_results, k = NULL, file_name = NULL, plot = TRUE)
```

## Arguments

- gsea_results:

  Output of
  [`compute_factor_gsea()`](https://verapancaldilab.github.io/CellTFusion/reference/compute_factor_gsea.md)
  (a list with a `GSEA_results` element, one table per NMF factor).

- k:

  Integer. Number of meta-programs (clusters) to extract. If `NULL`
  (default), chosen at the elbow of the within-cluster sum of squares
  (the k with the largest second difference) for k = 2 to min(10,
  n_hallmarks - 1).

- file_name:

  Optional character suffix for the saved heatmap
  (`Results/TCGA_meta_programs_<file_name>.pdf`).

- plot:

  Logical. If `TRUE` (default), saves a clustering heatmap.

## Value

A data frame with one row per meta-program and columns `meta_program`
(`"MP1"`, `"MP2"`, ...) and `hallmarks` (comma-separated Hallmark gene
set names).
