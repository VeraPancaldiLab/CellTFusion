# Plot cell-group dendrograms colored by cluster

Draws each cell-type dendrogram with its branches colored by the cluster
assignments from
[`calculate_dendrogram_cuts()`](https://verapancaldilab.github.io/CellTFusion/reference/calculate_dendrogram_cuts.md)
and saves all of them in a single PDF. The leaves are drawn in the order
of their correlation with the module as far as the tree allows (the two
sides of a branch are swapped, the tree itself is unchanged), so that
each cluster, a range of correlation values, is drawn in one place when
possible. A cluster that the tree splits into separate blocks of leaves
is drawn in each of them. The number of the cluster is written next to
each of its blocks: it is the number in the cell group name
(`Dendrogram_<module>_group_<number>`). Unassigned features are left
uncolored.

## Usage

``` r
plot_dendrogram_clusters(
  cell.group.dendrogram,
  cuts_per_dendrogram,
  file_name = NULL
)
```

## Arguments

- cell.group.dendrogram:

  A named list of `hclust` objects (or objects convertible to `hclust`),
  one per TF module, as returned by
  [`identify.cell.groups()`](https://verapancaldilab.github.io/CellTFusion/reference/identify.cell.groups.md)
  (the correlations they store are used to order the leaves; without
  them, the leaves keep the order of the tree).

- cuts_per_dendrogram:

  A list of integer cluster label vectors, one per dendrogram, as
  returned by
  [`calculate_dendrogram_cuts()`](https://verapancaldilab.github.io/CellTFusion/reference/calculate_dendrogram_cuts.md).

- file_name:

  Optional character. Suffix of the output file
  (`Results/Dendrogram_color_clusters_<file_name>.pdf`). If `NULL`
  (default), nothing is saved.

## Value

Called for its side effect (saves a PDF); returns `NULL` invisibly.
