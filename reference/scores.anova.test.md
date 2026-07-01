# One-way ANOVA test for multi-group comparisons

Performs one-way ANOVA to test for differences in scores across multiple
levels of a trait. Tukey post-hoc tests are used to identify pairwise
differences and significance is visualized as annotated boxplots.

## Usage

``` r
scores.anova.test(scores, coldata, trait, pval = 0.05)
```

## Arguments

- scores:

  A list, NMF output from
  [`compute.latent_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.latent_factors.md),
  or a score matrix. When a list, the first element must be a samples x
  features score matrix.

- coldata:

  A data frame containing sample annotations including the grouping
  variable.

- trait:

  Character. Name of the column in `coldata` used for the grouping
  variable.

- pval:

  Numeric. P-value threshold for significance (default = 0.05).

## Value

A list of significant features or `NULL` if none are significant.
