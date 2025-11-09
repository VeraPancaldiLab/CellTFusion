# One-way ANOVA test for multi-group comparisons

Performs one-way ANOVA to test for differences in cell group scores
across multiple levels of a trait. Tukey post-hoc tests are used to
identify pairwise differences and significance is visualized as
annotated boxplots.

## Usage

``` r
cell.groups.anova.test(cell.groups, coldata, trait, pval = 0.05)
```

## Arguments

- cell.groups:

  A list containing cell group score data, metadata, and identifiers.

- coldata:

  A data frame containing sample annotations including the grouping
  variable.

- trait:

  Character. Name of the column in `coldata` used for the grouping
  variable.

- pval:

  Numeric. P-value threshold for significance (default = 0.05).

## Value

A list of significant cell groups or `NULL` if none are significant.
