# One-way ANOVA test for multi-group comparisons

Performs one-way ANOVA to test for differences in scores across multiple
levels of a trait. Tukey post-hoc tests are used to identify pairwise
differences and significance is visualized as annotated boxplots
(`Results/ANOVA_<trait>_<feature>.pdf`).

## Usage

``` r
scores.anova.test(scores, coldata, trait, pval = 0.05, padj = FALSE)
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
  variable (categorical, two or more groups).

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
res <- scores.anova.test(cell_groups$Cell_groups[, 1:5], traitdata.tuto, "Enrollment.IC")
#> 5 significant feature(s) after one-way ANOVA (p-value <= 0.05): Dendrogram_pink_group_1, Dendrogram_pink_group_2, Dendrogram_pink_group_3, Dendrogram_pink_group_4, Dendrogram_pink_group_5
res$Statistics
#>                   feature          test statistic_name statistic      p_value
#> 1 Dendrogram_pink_group_1 one-way ANOVA              F  6.077982 2.765541e-03
#> 2 Dendrogram_pink_group_2 one-way ANOVA              F  4.343054 1.431839e-02
#> 3 Dendrogram_pink_group_3 one-way ANOVA              F 14.930184 9.546563e-07
#> 4 Dendrogram_pink_group_4 one-way ANOVA              F  8.533203 2.831134e-04
#> 5 Dendrogram_pink_group_5 one-way ANOVA              F  5.808542 3.563466e-03
#>     p_adjusted significant
#> 1 2.765541e-03        TRUE
#> 2 1.431839e-02        TRUE
#> 3 9.546563e-07        TRUE
#> 4 2.831134e-04        TRUE
#> 5 3.563466e-03        TRUE
```
