# Wilcoxon rank-sum test for binary traits

Performs a Wilcoxon rank-sum (Mann–Whitney U) test comparing cell group
scores between two levels of a binary clinical trait. Significant
features are plotted as boxplots and saved to the "Results/" folder.

## Usage

``` r
cell.groups.wilcox.test(cell.groups, coldata, trait, pval = 0.05)
```

## Arguments

- cell.groups:

  A list containing cell group scores and associated metadata.

- coldata:

  A data frame containing sample annotations and clinical traits.

- trait:

  Character. Name of the column in `coldata` used as the binary grouping
  variable.

- pval:

  Numeric. P-value threshold for significance (default = 0.05).

## Value

A list containing significant cell groups or `NULL` if none are found.
