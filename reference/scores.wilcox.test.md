# Wilcoxon rank-sum test for binary traits

Performs a Wilcoxon rank-sum (Mann-Whitney U) test comparing scores
between two levels of a binary clinical trait. Significant features are
plotted as boxplots and saved to the "Results/" folder
(`Results/Wilcoxon_<trait>_<feature>.pdf`).

## Usage

``` r
scores.wilcox.test(scores, coldata, trait, pval = 0.05, p_adjust = "none")
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
res <- scores.wilcox.test(cell_groups$Cell_groups[, 1:5], traitdata.tuto, "binaryResponse")
#> 1 significant feature(s) after Wilcoxon rank-sum test (p-value <= 0.05): Dendrogram_brown_group_4
res$Statistics
#>                    feature                   test statistic_name statistic
#> 1 Dendrogram_brown_group_1 Wilcoxon rank-sum test              W      2137
#> 2 Dendrogram_brown_group_2 Wilcoxon rank-sum test              W      2034
#> 3 Dendrogram_brown_group_3 Wilcoxon rank-sum test              W      2577
#> 4 Dendrogram_brown_group_4 Wilcoxon rank-sum test              W      2855
#> 5 Dendrogram_brown_group_5 Wilcoxon rank-sum test              W      2551
#>       p_value  p_adjusted significant
#> 1 0.850018441 0.850018441       FALSE
#> 2 0.837938401 0.837938401       FALSE
#> 3 0.059148236 0.059148236       FALSE
#> 4 0.003077249 0.003077249        TRUE
#> 5 0.073977286 0.073977286       FALSE
```
