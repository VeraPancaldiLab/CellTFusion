# Run multivariate feature-based GSEA using limma and Hallmark gene sets

This function fits a multivariate linear model for each gene using all
features in `scores` as continuous covariates (each feature
standardized, and the cohort added when `batch` is given). For each
feature, it extracts the moderated t-statistics and p-values, ranks
genes, and performs GSEA (`fgsea`) using the human Hallmark gene sets
from MSigDB. Optional dotplots for the top enriched pathways can be
saved as PDFs.

## Usage

``` r
compute_factor_gsea(
  counts_norm,
  scores,
  return = TRUE,
  top_n = 10,
  file_name = NULL,
  width = 8,
  height = 10,
  batch = NULL,
  seed = 123
)
```

## Arguments

- counts_norm:

  A numeric matrix or data frame of log-normalized gene expression
  values (e.g. log2(TPM + 1); genes as rows, gene symbols as row names,
  samples as columns).

- scores:

  A data frame of continuous features (samples in rows, features in
  columns) to be modeled as covariates. Rows must be the same samples,
  in the same order, as the columns of `counts_norm` (checked by name;
  an error is raised otherwise). Features without variation are left out
  with a message.

- return:

  Logical; if TRUE, generates and saves dotplots of top enriched
  Hallmark pathways for each feature. Default is TRUE.

- top_n:

  Integer; number of top pathways to display in the dotplot. Default is
  10.

- file_name:

  Character; optional suffix for saved PDF files
  (`Results/GSEA_<feature>_<file_name>.pdf`, or
  `Results/GSEA_<feature>.pdf` without it). Default is NULL.

- width:

  Numeric; width of the PDF plot in inches. Default is 8.

- height:

  Numeric; height of the PDF plot in inches. Default is 10.

- batch:

  Optional vector with the cohort of each sample (in the row order of
  `scores`, or named by sample). With two or more cohorts, the cohort is
  added to the model as a categorical covariate, so that expression
  differences between cohorts are not attributed to the features.
  Default `NULL`.

- seed:

  Integer. Random seed used by
  [`fgsea::fgsea()`](https://rdrr.io/pkg/fgsea/man/fgsea.html) (its
  p-values come from random permutations), so that the results are
  reproducible. Default 123. The caller's random number generator state
  is not changed.

## Value

A list containing:

- DE_results:

  A named list of `topTable` results for each feature, including logFC
  (change in expression per standard deviation of the feature),
  moderated t-statistics, p-values, and adjusted p-values.

- GSEA_results:

  A named list of
  [`fgsea::fgsea()`](https://rdrr.io/pkg/fgsea/man/fgsea.html) result
  tables (columns include `pathway`, `pval`, `padj`, `NES`) for each
  feature, ordered by p-value. fgsea may not be able to estimate the
  p-value of some pathways: their `pval`, `padj` and `NES` are then `NA`
  (their `ES` is kept).

## Details

The function works as follows:

1.  Hallmark gene sets are retrieved from MSigDB using `msigdbr`.

2.  A multivariate linear model is fitted for each gene using
    [`limma::lmFit`](https://rdrr.io/pkg/limma/man/lmFit.html): all the
    features are in the same model, so the effect of each feature is
    adjusted for the others (and for the cohort with `batch`).

3.  Empirical Bayes moderation is applied via
    [`limma::eBayes`](https://rdrr.io/pkg/limma/man/ebayes.html).

4.  For each feature:

    1.  Differential expression results are extracted using `topTable`
        for the coefficient of that feature.

    2.  Genes are ranked by moderated t-statistics.

    3.  Hallmark GSEA is performed on the ranked gene list with
        [`fgsea::fgsea()`](https://rdrr.io/pkg/fgsea/man/fgsea.html).

    4.  Optionally, a dotplot of the top enriched pathways is generated.

## Examples

``` r
if (FALSE) { # \dontrun{
  results <- compute_factor_gsea(counts_norm = expression_matrix,
                              scores = feature_table,
                              return = TRUE,
                              top_n = 10,
                              file_name = "Feature1_vs_all")
} # }
```
