# Merge highly correlated TF modules

Identifies pairs of TF modules whose eigengene correlation exceeds
`corr` and merges them by averaging their columns.

## Usage

``` r
mergeModules(data, colors, corr)
```

## Arguments

- data:

  A data frame of TF module eigengenes (samples x modules), with columns
  named `ME<color>`.

- colors:

  A character vector with the module color label of each TF.

- corr:

  Numeric. Spearman correlation threshold above which two modules are
  merged.

## Value

A list of two elements:

- `[[1]]`: Data frame of merged module eigengenes (samples x modules).

- `[[2]]`: Updated character vector of module color labels after
  merging.
