# Identify cell groups

Identifies cell dendrograms corresponding to each TF module group based
on correlation with deconvolution features.

## Usage

``` r
# S3 method for class 'cell.groups'
identify(
  features,
  clustering.method = "ward.D2",
  width = 12,
  height = 18,
  return = T
)
```

## Arguments

- features:

  A list of two elements (accessed by position), as returned by
  [`compute.modules.relationship()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.modules.relationship.md)
  with `return = TRUE`:

  - `[[1]]`: A matrix of correlations between TF modules (rows) and cell
    type features (columns).

  - `[[2]]`: A named list of significant features per TF module (e.g.,
    p-value \< 0.05).

- clustering.method:

  Clustering method used in hclust. Default: "ward.D2".

- width:

  Width (in inches) of the saved PDF plots. Default: 12.

- height:

  Height (in inches) of the saved PDF plots. Default: 18.

- return:

  Logical; whether to save dendrogram plots to the "Results/" folder.
  Default: TRUE.

## Value

A named list of dendrograms (`hclust` objects), one per TF module.
Modules with fewer than two significant features are discarded; `NULL`
is returned if no module remains.
