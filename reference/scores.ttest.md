# T-test for cell group comparisons

Performs a two-sample t-test (Welch's, unequal variances) comparing cell
group scores between two groups of a binary trait. Significant features
are plotted as boxplots and saved as PDF files in the "Results/"
directory (`Results/Ttest_<trait>_<feature>.pdf`).

## Usage

``` r
scores.ttest(scores, coldata, trait, pval = 0.05, p_adjust = "none")
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

A list with the significant features: its first element is the score
matrix restricted to them (no column if none is significant), followed
by the matching parts of the other elements of `scores` (e.g.
`Composition`, `Weights` and `Module_relation` of
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
#> 20 of 764 features were not assigned to any cluster by the tree cut and are not part of any cell group (blue: 6, pink: 2, red: 2, yellow: 8, black: 2)
res <- scores.ttest(cell_groups$Cell_groups[, 1:5], traitdata.tuto, "binaryResponse")
#> 2 significant feature(s) after Welch's t-test (p-value <= 0.05): Dendrogram_brown_group_3, Dendrogram_brown_group_4
res$Statistics
#>                    feature           test statistic_name  statistic
#> 1 Dendrogram_brown_group_1 Welch's t-test              t  0.2040300
#> 2 Dendrogram_brown_group_2 Welch's t-test              t -0.3177821
#> 3 Dendrogram_brown_group_3 Welch's t-test              t  2.4972292
#> 4 Dendrogram_brown_group_4 Welch's t-test              t  4.6047254
#> 5 Dendrogram_brown_group_5 Welch's t-test              t  1.7359607
#>        p_value   p_adjusted significant
#> 1 8.394441e-01 8.394441e-01       FALSE
#> 2 7.527480e-01 7.527480e-01       FALSE
#> 3 1.676793e-02 1.676793e-02        TRUE
#> 4 2.331147e-05 2.331147e-05        TRUE
#> 5 9.191050e-02 9.191050e-02       FALSE
```
