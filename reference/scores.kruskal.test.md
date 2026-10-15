# Kruskal-Wallis test for multi-group comparisons

Performs a Kruskal-Wallis test to compare scores across multiple trait
levels. Significant results are visualized as annotated boxplots with
Dunn post-hoc tests (BH-adjusted), saved as
`Results/Kruskal_<trait>_<feature>.pdf`.

## Usage

``` r
scores.kruskal.test(scores, coldata, trait, pval = 0.05, padj = FALSE)
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
res <- scores.kruskal.test(cell_groups$Cell_groups[, 1:5], traitdata.tuto, "Enrollment.IC")
#> 5 significant feature(s) after Kruskal-Wallis test (p-value <= 0.05): Dendrogram_pink_group_1, Dendrogram_pink_group_2, Dendrogram_pink_group_3, Dendrogram_pink_group_4, Dendrogram_pink_group_5
res$Statistics
#>                   feature                test statistic_name statistic
#> 1 Dendrogram_pink_group_1 Kruskal-Wallis test              H 10.091884
#> 2 Dendrogram_pink_group_2 Kruskal-Wallis test              H  8.075064
#> 3 Dendrogram_pink_group_3 Kruskal-Wallis test              H 29.615321
#> 4 Dendrogram_pink_group_4 Kruskal-Wallis test              H 15.049547
#> 5 Dendrogram_pink_group_5 Kruskal-Wallis test              H  8.563492
#>        p_value   p_adjusted significant
#> 1 6.435395e-03 6.435395e-03        TRUE
#> 2 1.764095e-02 1.764095e-02        TRUE
#> 3 3.707787e-07 3.707787e-07        TRUE
#> 4 5.395509e-04 5.395509e-04        TRUE
#> 5 1.381852e-02 1.381852e-02        TRUE
```
