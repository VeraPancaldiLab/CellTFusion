# Wilcoxon rank-sum test for binary traits

Performs a Wilcoxon rank-sum (Mann-Whitney U) test comparing scores
between two levels of a binary clinical trait. Significant features are
plotted as boxplots and saved to the "Results/" folder.

## Usage

``` r
scores.wilcox.test(scores, coldata, trait, pval = 0.05)
```

## Arguments

- scores:

  A list whose first element is a samples x features score matrix, e.g.
  the output of
  [`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md)
  or
  [`compute.latent_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.latent_factors.md).
  To test a plain matrix, use
  [`scores.stat.analysis()`](https://verapancaldilab.github.io/CellTFusion/reference/scores.stat.analysis.md).

- coldata:

  A data frame containing sample annotations and clinical traits.

- trait:

  Character. Name of the column in `coldata` used as the binary grouping
  variable.

- pval:

  Numeric. P-value threshold for significance (default = 0.05).

## Value

A list containing significant features or `NULL` if none are found.
