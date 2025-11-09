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

  A list with two named elements:

  correlations

  :   A matrix of correlations between TF modules and cell type
      features.

  significant

  :   A named list of significant features per TF module (e.g., p-value
      \< 0.05).

- clustering.method:

  Clustering method used in hclust. Default: "ward.D2".

- width:

  Width of the saved PDF plots.

- height:

  Height of the saved PDF plots.

- return:

  Logical; whether to save dendrogram plots to the "Results/" folder.

## Value

A named list of dendrograms for each TF module.
