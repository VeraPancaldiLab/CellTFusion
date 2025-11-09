# Kruskal–Wallis test for multi-group comparisons

Performs a Kruskal–Wallis test to compare cell group scores across
multiple trait levels. Significant results are visualized as annotated
boxplots with Dunn post-hoc tests.

## Usage

``` r
cell.groups.kruskal.test(cell.groups, coldata, trait, pval = 0.05)
```

## Arguments

- cell.groups:

  A list containing cell group score data, metadata, and identifiers.

- coldata:

  A data frame containing sample annotations, including the grouping
  trait.

- trait:

  Character. Name of the column in `coldata` used as the grouping
  variable.

- pval:

  Numeric. P-value threshold for significance (default = 0.05).

## Value

A list containing only significant cell groups after Kruskal–Wallis
test. Returns `NULL` if no significant groups are found.
