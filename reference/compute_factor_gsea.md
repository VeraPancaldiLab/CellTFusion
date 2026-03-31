# Run multivariate feature-based GSEA using limma and Hallmark gene sets

This function fits a multivariate linear model for each gene using all
features in `features_df` as continuous covariates. For each feature, it
extracts the moderated t-statistics and p-values, ranks genes, and
performs GSEA using the Hallmark gene sets from MSigDB. Optional
dotplots for the top enriched pathways can be saved as PDFs.

## Usage

``` r
compute_factor_gsea(
  RNA.tpm,
  features_df,
  plot_dot = TRUE,
  top_n = 10,
  file_name = NULL,
  width = 8,
  height = 10
)
```

## Arguments

- RNA.tpm:

  A numeric matrix or data frame of gene expression values (genes in
  rows, samples in columns).

- features_df:

  A data frame of continuous features (samples in rows, features in
  columns) to be modeled as covariates.

- plot_dot:

  Logical; if TRUE, generates and saves dotplots of top enriched
  Hallmark pathways for each feature. Default is TRUE.

- top_n:

  Integer; number of top pathways to display in the dotplot. Default is
  10.

- file_name:

  Character; optional suffix for saved PDF files. Default is NULL.

- width:

  Numeric; width of the PDF plot in inches. Default is 8.

- height:

  Numeric; height of the PDF plot in inches. Default is 10.

## Value

A list containing:

- DE_results:

  A named list of `topTable` results for each feature, including logFC,
  moderated t-statistics, p-values, and adjusted p-values.

- GSEA_results:

  A named list of `GSEA` results from `clusterProfiler` for each
  feature.

## Details

The function works as follows:

1.  Hallmark gene sets are retrieved from MSigDB using `msigdbr`.

2.  A multivariate linear model is fitted for each gene using
    [`limma::lmFit`](https://rdrr.io/pkg/limma/man/lmFit.html).

3.  Empirical Bayes moderation is applied via
    [`limma::eBayes`](https://rdrr.io/pkg/limma/man/ebayes.html).

4.  For each feature:

    1.  Differential expression results are extracted using `topTable`
        for the coefficient of that feature.

    2.  Genes are ranked by moderated t-statistics.

    3.  Hallmark GSEA is performed using the ranked gene list.

    4.  Optionally, a dotplot of the top enriched pathways is generated.

## Examples

``` r
if (FALSE) { # \dontrun{
  results <- compute_factor_gsea(RNA.tpm = expression_matrix,
                              features_df = feature_table,
                              plot_dot = TRUE,
                              top_n = 10,
                              file_name = "Feature1_vs_all")
} # }
```
