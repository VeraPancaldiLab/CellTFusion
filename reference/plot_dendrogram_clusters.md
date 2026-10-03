# Plot cell-group dendrograms colored by cluster

Draws each cell-type dendrogram with its branches colored by the cluster
assignments from
[`calculate_dendrogram_cuts()`](https://verapancaldilab.github.io/CellTFusion/reference/calculate_dendrogram_cuts.md)
and saves all of them in a single PDF.

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

  A named list of `hclust`-convertible dendrogram objects, one per TF
  module, as returned by
  [`identify.cell.groups()`](https://verapancaldilab.github.io/CellTFusion/reference/identify.cell.groups.md).

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
