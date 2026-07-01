# Calculate dendrogram cut heights

Computes a sequence of candidate cut heights for each cell-type
dendrogram. Heights are distributed between a buffered minimum and
maximum derived from the dendrogram's own height distribution, avoiding
trivial cuts (single-element or all-in-one clusters).

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

  A list of `hclust` objects, one per TF module, as returned by
  [`identify.cell.groups()`](https://verapancaldilab.github.io/CellTFusion/reference/identify.cell.groups.md).

- n_cuts:

  Integer. Number of evenly spaced cut heights to generate per
  dendrogram. If `NULL` (default), the number is set proportional to the
  maximum dendrogram height.

## Value

A list of numeric vectors, one per dendrogram, containing the candidate
cut heights.
