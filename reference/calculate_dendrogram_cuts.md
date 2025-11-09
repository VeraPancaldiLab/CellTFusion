# Calculate dendrogram cuts

Calculate dendrogram cuts

## Usage

``` r
calculate_dendrogram_cuts(cell.group.dendrogram, n_cuts = NULL)
```

## Arguments

- cell.group.dendrogram:

  List with the cell dendrograms corresponding to each TF module
  obtained from identify.cell.groups()

- n_cuts:

  Optional parameter to limit the number of cuts the dendrogram needs to
  be cut (Default is NULL). If no parameter is set, number of cuts will
  be proportional to the height of the dendrogram.

## Value

A list with the sequence of numbers where each dendrogram will be cut
