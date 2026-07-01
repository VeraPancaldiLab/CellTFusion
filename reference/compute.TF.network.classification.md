# Compute TF Network Classification

Classifies transcription factor (TF) modules into clusters based on
their correlations with pathway activity values across samples. Uses
hierarchical clustering and silhouette width to determine the optimal
number of clusters.

## Usage

``` r
# S3 method for class 'TF.network.classification'
compute(tf.network, pathways.features, return = T)
```

## Arguments

- tf.network:

  A list containing TF module eigengenes, typically output from
  [`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md)

- pathways.features:

  A matrix with pathway activities, typically from
  [`compute.pathway.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.pathway.activity.md)

- return:

  Logical. If TRUE, intermediate plots (e.g. silhouette, dendrogram,
  PCA) are saved in the `Results/` directory. Default is TRUE.

## Value

A named list of TF module clusters.

## Examples

``` r
if (FALSE) { # \dontrun{
data("network.tuto")
data("counts.norm.tuto")
pathways <- compute.pathway.activity(counts.norm.tuto)
tfs.modules.clusters <- compute.TF.network.classification(tf.network = network.tuto,
                                                          pathways.features = pathways,
                                                          return = FALSE)
} # }
```
