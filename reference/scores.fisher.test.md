# Fisher's exact test for score-trait association

Performs a Fisher's exact test between each score, binarised at its
median into High/Low groups (a score equal to the median is "Low"), and
a categorical trait. Significant features are plotted as barplots and
saved to the "Results/" folder (`Results/Fisher_<trait>_<feature>.pdf`).
The statistic reported is the odds ratio for a two-group trait (`NA`
otherwise).

## Usage

``` r
scores.fisher.test(scores, coldata, trait, pval = 0.05, padj = FALSE)
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

  Character. Name of the column in `coldata` to test (categorical, two
  or more groups).

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
res <- scores.fisher.test(cell_groups$Cell_groups[, 1:5], traitdata.tuto, "binaryResponse")
#> 1 significant feature(s) after Fisher's exact test (p-value <= 0.05): Dendrogram_pink_group_3
res$Statistics
#>                   feature                test statistic_name statistic
#> 1 Dendrogram_pink_group_1 Fisher's exact test     odds ratio 2.3567205
#> 2 Dendrogram_pink_group_2 Fisher's exact test     odds ratio 2.3567205
#> 3 Dendrogram_pink_group_3 Fisher's exact test     odds ratio 2.9182035
#> 4 Dendrogram_pink_group_4 Fisher's exact test     odds ratio 1.9267619
#> 5 Dendrogram_pink_group_5 Fisher's exact test     odds ratio 0.7590797
#>      p_value p_adjusted significant
#> 1 0.08469494 0.08469494       FALSE
#> 2 0.08469494 0.08469494       FALSE
#> 3 0.03039497 0.03039497        TRUE
#> 4 0.19751821 0.19751821       FALSE
#> 5 0.66869828 0.66869828       FALSE
```
