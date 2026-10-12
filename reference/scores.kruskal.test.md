# Kruskal-Wallis test for multi-group comparisons

Performs a Kruskal-Wallis test to compare scores across multiple trait
levels. Significant results are visualized as annotated boxplots with
Dunn post-hoc tests (BH-adjusted), saved as
`Results/Kruskal_<trait>_<feature>.pdf`.

## Usage

``` r
scores.kruskal.test(scores, coldata, trait, pval = 0.05, p_adjust = "none")
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

- p_adjust:

  Character. Method to correct the p-values for multiple testing across
  the features (see
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html), e.g.
  `"BH"`). Default `"none"`.

## Value

See
[`scores.ttest()`](https://verapancaldilab.github.io/CellTFusion/reference/scores.ttest.md).

## Examples

``` r
data("network.tuto")
data("deconv_subgroups.tuto")
data("traitdata.tuto")
cell_groups <- construct_cell_groups(network.tuto, deconv_subgroups.tuto)
#> 20 of 764 features were not assigned to any cluster by the tree cut and are not part of any cell group (blue: 6, pink: 2, red: 2, yellow: 8, black: 2)
res <- scores.kruskal.test(cell_groups$Cell_groups[, 1:5], traitdata.tuto, "Enrollment.IC")
#> 3 significant feature(s) after Kruskal-Wallis test (p-value <= 0.05): Dendrogram_brown_group_3, Dendrogram_brown_group_4, Dendrogram_brown_group_5
res$Statistics
#>                    feature                test statistic_name statistic
#> 1 Dendrogram_brown_group_1 Kruskal-Wallis test              H  2.695207
#> 2 Dendrogram_brown_group_2 Kruskal-Wallis test              H  2.612820
#> 3 Dendrogram_brown_group_3 Kruskal-Wallis test              H 10.356225
#> 4 Dendrogram_brown_group_4 Kruskal-Wallis test              H 12.480027
#> 5 Dendrogram_brown_group_5 Kruskal-Wallis test              H 26.813980
#>        p_value   p_adjusted significant
#> 1 2.598623e-01 2.598623e-01       FALSE
#> 2 2.707904e-01 2.707904e-01       FALSE
#> 3 5.638640e-03 5.638640e-03        TRUE
#> 4 1.949829e-03 1.949829e-03        TRUE
#> 5 1.504590e-06 1.504590e-06        TRUE
```
