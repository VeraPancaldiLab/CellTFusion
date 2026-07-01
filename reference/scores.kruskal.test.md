# Kruskal-Wallis test for multi-group comparisons

Performs a Kruskal-Wallis test to compare scores across multiple trait
levels. Significant results are visualized as annotated boxplots with
Dunn post-hoc tests.

## Usage

``` r
scores.kruskal.test(scores, coldata, trait, pval = 0.05)
```

## Arguments

- scores:

  A list, NMF output from
  [`compute.latent_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.latent_factors.md),
  or a score matrix. When a list, the first element must be a samples x
  features score matrix.

- coldata:

  A data frame containing sample annotations, including the grouping
  trait.

- trait:

  Character. Name of the column in `coldata` used as the grouping
  variable.

- pval:

  Numeric. P-value threshold for significance (default = 0.05).

## Value

A list containing only significant features after Kruskal-Wallis test.
Returns `NULL` if no significant features are found.
