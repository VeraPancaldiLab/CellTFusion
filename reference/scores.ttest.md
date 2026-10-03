# T-test for cell group comparisons

Performs a two-sample t-test (Welch's, unequal variances) comparing cell
group scores between two groups of a binary trait. Significant features
are plotted as boxplots and saved as PDF files in the "Results/"
directory.

## Usage

``` r
scores.ttest(scores, coldata, trait, pval = 0.05)
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

  A data frame containing sample-level annotations including the trait
  to test.

- trait:

  Character. Name of the column in `coldata` used as the grouping
  variable.

- pval:

  Numeric. P-value threshold for significance (default = 0.05).

## Value

A list containing only significant cell groups after the t-test. Returns
`NULL` if no significant groups are found.
