# T-test for cell group comparisons

Performs a two-sample t-test (Welch's, unequal variances) comparing cell
group scores between two groups of a binary trait. Significant features
are plotted as boxplots and saved as PDF files in the "Results/"
directory (`Results/Ttest_<trait>_<feature>.pdf`).

## Usage

``` r
scores.ttest(scores, coldata, trait, pval = 0.05, padj = FALSE)
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

A list with the significant features: its first element is the score
matrix restricted to them (no column if none is significant), followed
by the matching parts of the other elements of `scores` (e.g.
`Composition` and `Weights` of
[`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md)),
with their names. The last element, `Statistics`, is a data frame with
the test of every feature (`feature`, `test`, `statistic_name`,
`statistic`, `p_value`, `p_adjusted`, `significant`).

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
res <- scores.ttest(cell_groups$Cell_groups[, 1:5], traitdata.tuto, "binaryResponse")
#> 3 significant feature(s) after Welch's t-test (p-value <= 0.05): Dendrogram_pink_group_1, Dendrogram_pink_group_2, Dendrogram_pink_group_3
res$Statistics
#>                   feature           test statistic_name   statistic     p_value
#> 1 Dendrogram_pink_group_1 Welch's t-test              t  2.78212854 0.008297032
#> 2 Dendrogram_pink_group_2 Welch's t-test              t  2.17742792 0.036433660
#> 3 Dendrogram_pink_group_3 Welch's t-test              t  2.72545759 0.010177443
#> 4 Dendrogram_pink_group_4 Welch's t-test              t  1.53700466 0.134682645
#> 5 Dendrogram_pink_group_5 Welch's t-test              t -0.06210619 0.950920316
#>    p_adjusted significant
#> 1 0.008297032        TRUE
#> 2 0.036433660        TRUE
#> 3 0.010177443        TRUE
#> 4 0.134682645       FALSE
#> 5 0.950920316       FALSE
```
