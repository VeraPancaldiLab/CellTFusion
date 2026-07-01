# Student's t-test for cell group comparisons

Performs a Student's t-test comparing cell group scores between two
groups of a binary trait. Significant features are plotted as boxplots
and saved as PDF files in the "Results/" directory.

## Usage

``` r
scores.ttest(scores, coldata, trait, pval = 0.05)
```

## Arguments

- scores:

  A list or matrix of cell group scores. When a list, the first element
  must be a data frame or matrix of scores (samples x features).

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
