# Calculate dendrogram cluster assignments

Cuts each cell-type dendrogram into clusters (the cell groups of a TF
module) in two steps:

1.  the dendrogram is cut with dynamic tree cutting
    ([`dynamicTreeCut::cutreeDynamic()`](https://rdrr.io/pkg/dynamicTreeCut/man/cutreeDynamic.html),
    `method = "tree"`, `deepSplit = TRUE`, which splits clusters more
    finely than `FALSE`);

2.  each cluster is made to cover a single range of correlation values
    with the module. With `method = "tree"`, a cluster can be a branch
    from which smaller branches were split off as other clusters, so its
    features can lie in separated ranges of correlation (e.g. 0.14,
    0.17-0.19 and 0.26-0.39, with other clusters in between). The
    assigned features are sorted by their correlation, and each run of
    consecutive features with the same cluster becomes a piece. A piece
    with fewer than `min_cluster_size` features is joined to the
    neighbouring piece (just below or just above in correlation) with
    the smallest gap between their correlations, the smallest gap first,
    until every piece has at least `min_cluster_size` features. The
    pieces are the final clusters, numbered from the lowest to the
    highest correlation.

Features left unassigned by the dynamic tree cut stay unassigned. Step 2
needs the correlations stored in the dendrograms by
[`identify.cell.groups()`](https://verapancaldilab.github.io/CellTFusion/reference/identify.cell.groups.md);
for dendrograms without them, only step 1 is done (with a message).

## Usage

``` r
calculate_dendrogram_cuts(cell_dendrograms, min_cluster_size = 3)
```

## Arguments

- cell_dendrograms:

  A list of `hclust` objects (or objects convertible to `hclust`), one
  per TF module, as returned by
  [`identify.cell.groups()`](https://verapancaldilab.github.io/CellTFusion/reference/identify.cell.groups.md).

- min_cluster_size:

  Integer. Minimum number of features per cluster: passed to
  [`dynamicTreeCut::cutreeDynamic()`](https://rdrr.io/pkg/dynamicTreeCut/man/cutreeDynamic.html)'s
  `minClusterSize` argument, and used in step 2. Default is 3.

## Value

A list of integer cluster label vectors, one per dendrogram, in the same
order as `cell_dendrograms` and named like it, each in the order of the
dendrogram labels (label `0` marks unassigned elements).
