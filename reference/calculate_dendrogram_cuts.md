# Calculate dendrogram cluster assignments

Cuts each cell-type dendrogram into clusters using dynamic tree cutting
([`dynamicTreeCut::cutreeDynamic()`](https://rdrr.io/pkg/dynamicTreeCut/man/cutreeDynamic.html),
`method = "tree"`).

## Usage

``` r
calculate_dendrogram_cuts(
  cell.group.dendrogram,
  deep_split = 4,
  min_cluster_size = 3
)
```

## Arguments

- cell.group.dendrogram:

  A list of `hclust`-convertible dendrogram objects, one per TF module,
  as returned by
  [`identify.cell.groups()`](https://verapancaldilab.github.io/CellTFusion/reference/identify.cell.groups.md).

- deep_split:

  Integer. Passed to
  [`dynamicTreeCut::cutreeDynamic()`](https://rdrr.io/pkg/dynamicTreeCut/man/cutreeDynamic.html)'s
  `deepSplit` argument; controls the sensitivity of cluster splitting.
  Default is 4.

- min_cluster_size:

  Integer. Passed to
  [`dynamicTreeCut::cutreeDynamic()`](https://rdrr.io/pkg/dynamicTreeCut/man/cutreeDynamic.html)'s
  `minClusterSize` argument; minimum number of elements per cluster.
  Default is 3.

## Value

A list of integer cluster label vectors, one per dendrogram, in the same
order as `cell.group.dendrogram` (label `0` marks unassigned elements).
