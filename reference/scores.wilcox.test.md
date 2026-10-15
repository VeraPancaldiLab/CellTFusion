# Wilcoxon rank-sum test for binary traits

Performs a Wilcoxon rank-sum (Mann-Whitney U) test comparing scores
between two levels of a binary clinical trait. Significant features are
plotted as boxplots and saved to the "Results/" folder
(`Results/Wilcoxon_<trait>_<feature>.pdf`).

## Usage

``` r
scores.wilcox.test(scores, coldata, trait, pval = 0.05, padj = FALSE)
```

## Arguments

- scores:

  A samples x features score matrix, a list whose first element is that
  matrix (e.g. the output of
  [`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md)),
  or the output of
  [`compute.latent_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.latent_factors.md)
  (its factor scores `Z` are tested).

- coldata:

  A data frame (or tibble) containing sample-level annotations including
  the trait to test. Its rows are matched to the samples of the scores
  by row name; without sample names, they must be in the order of the
  scores. Samples with a missing trait are left out (with a message).

- trait:

  Character. Name of the column in `coldata` used as the grouping
  variable (categorical, two groups).

- pval:

  Numeric. Threshold for significance: a feature is significant if its
  (adjusted) p-value is `<= pval` (default = 0.05).

- padj:

  Logical; if TRUE, the p-values are corrected for multiple testing
  across the features with the Benjamini-Hochberg method, and the
  adjusted p-values are compared with the threshold. Default FALSE.

## Value

See
[`scores.ttest()`](https://verapancaldilab.github.io/CellTFusion/reference/scores.ttest.md).

## Examples

``` r
data("network.tuto")
data("deconv_subgroups.tuto")
data("traitdata.tuto")
cell_groups <- construct_cell_groups(network.tuto, deconv_subgroups.tuto)
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Dendrogram plots saved to: Results/Dendrogram_color_clusters.pdf
#> 44 of 897 features were not assigned to any cluster by the tree cut and are not part of any cell group (pink: 2, yellow: 8, red: 2, blue: 6, brown: 4, turquoise: 12, black: 4, green: 6)
res <- scores.wilcox.test(cell_groups$Cell_groups[, 1:5], traitdata.tuto, "binaryResponse")
#> 2 significant feature(s) after Wilcoxon rank-sum test (p-value <= 0.05): Dendrogram_pink_group_1, Dendrogram_pink_group_3
res$Statistics
#>                   feature                   test statistic_name statistic
#> 1 Dendrogram_pink_group_1 Wilcoxon rank-sum test              W      2730
#> 2 Dendrogram_pink_group_2 Wilcoxon rank-sum test              W      2556
#> 3 Dendrogram_pink_group_3 Wilcoxon rank-sum test              W      2800
#> 4 Dendrogram_pink_group_4 Wilcoxon rank-sum test              W      2453
#> 5 Dendrogram_pink_group_5 Wilcoxon rank-sum test              W      1978
#>       p_value  p_adjusted significant
#> 1 0.013229564 0.013229564        TRUE
#> 2 0.070910652 0.070910652       FALSE
#> 3 0.006002137 0.006002137        TRUE
#> 4 0.158964122 0.158964122       FALSE
#> 5 0.674018222 0.674018222       FALSE
```
