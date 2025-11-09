# Extract significant features using Wilcoxon test

Performs Wilcoxon rank-sum test for each feature comparing groups
defined by the trait.

## Usage

``` r
extract_wilcox_significant(features, trait)
```

## Arguments

- features:

  A numeric data frame or matrix where columns are features and rows are
  samples.

- trait:

  A vector or factor defining group labels for each sample.

## Value

A character vector of feature names with significant difference between
trait groups (p \< 0.05).
