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

  A list, NMF output from
  [`compute.latent_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.latent_factors.md),
  or a score matrix. When a list, the first element must be a samples x
  features score matrix.

- coldata:

  A data frame containing sample annotations and clinical traits.

- trait:

  Character. Name of the column in `coldata` used as the binary grouping
  variable.

- pval:

  Numeric. P-value threshold for significance (default = 0.05).

## Value

A list containing significant features or `NULL` if none are found.
