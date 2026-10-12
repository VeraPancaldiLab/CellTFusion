# Merge highly correlated TF modules

Merges TF modules whose scores (eigengenes) are too similar, repeating
until no pair of modules has a correlation above `corr`: at each step
the most correlated pair above `corr` is merged, and the correlations
are computed again. The merged module keeps the name and color of the
first module of the pair; its score is the mean of the scores of all the
original modules it contains, scaled (mean 0, standard deviation 1).

## Usage

``` r
mergeModules(data, colors, corr, method = "spearman")
```

## Arguments

- data:

  A data frame of TF module eigengenes (samples x modules), with columns
  named `ME<color>`.

- colors:

  A character vector with the module color label of each TF.

- corr:

  Numeric. Correlation threshold above which two modules are merged.

- method:

  Character. Correlation used: "pearson" or "spearman" (default).

## Value

A list of two elements:

- `[[1]]`: Data frame of merged module eigengenes (samples x modules).

- `[[2]]`: Updated character vector of module color labels after
  merging.
