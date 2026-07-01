# Derive TME meta-programs by clustering Hallmarks across NMF factors

Hierarchically clusters Hallmark gene sets by their NES profile across
NMF factors to identify recurrent transcriptional programs in the TME.
Optionally annotates each NMF factor with a Bagaev et al. (2021) MFP
subtype (IE, IE/F, F, D) when `Z`, `annot`, and `cancer_name` are all
supplied.

## Usage

``` r
derive_meta_programs(gsea_results, k = NULL, file_name = NULL, plot = TRUE)
```

## Arguments

- gsea_results:

  A list of GSEA result data frames (one per NMF factor), as returned by
  [`compute_factor_gsea()`](https://verapancaldilab.github.io/CellTFusion/reference/compute_factor_gsea.md).

- k:

  Integer. Number of meta-programs (clusters) to extract. If `NULL`
  (default), estimated automatically from the dendrogram.

- file_name:

  Optional character. File path prefix for saving output plots.

- plot:

  Logical. If `TRUE` (default), saves a clustering heatmap.

## Value

A list with:

- `meta_programs`: Named list mapping meta-program labels to character
  vectors of Hallmark names.

- `hallmark_clusters`: Data frame with columns `Hallmark`,
  `meta_program`, and `mean_NES`.

- `heatmap`: The `pheatmap` object.

- `k`: The number of meta-programs used.

- `factor_mfp`: (Only when `Z`, `annot`, and `cancer_name` are provided)
  Data frame with columns `factor`, `best_MFP`, `cor_IE`, `cor_IEF`,
  `cor_F`, `cor_D`, `n_samples`.
