# Identify cell-type niches from NMF latent factors

For each NMF factor, extracts the top-contributing cell groups (by basis
weight), maps them to their constituent cell types, and keeps only cell
types enriched relative to the background composition. Saves a
star-network PDF per factor and returns the weighted cell-type
associations.

## Usage

``` r
compute_cells_niches(
  latent_factors,
  dt,
  cell.groups,
  enrich_thresh = 1.5,
  quantile_cutoff = 0.7,
  cells_extra = NULL,
  return = TRUE,
  file_name = NULL
)
```

## Arguments

- latent_factors:

  A list returned by
  [`compute.latent_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.latent_factors.md),
  containing `W` (features x factors) and the NMF model object.

- dt:

  A named list of deconvolution subgroup results, used to build the
  cell-group composition matrix via
  [`compute.composition.matrix()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.composition.matrix.md).

- cell.groups:

  A list of cell group definitions as returned by
  [`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md).

- enrich_thresh:

  Numeric. Minimum enrichment ratio (foreground / background frequency)
  for a cell type to be retained per factor. Default 1.5.

- quantile_cutoff:

  Numeric between 0 and 1. Quantile threshold for selecting
  top-contributing cell groups per factor. Default 0.7.

- cells_extra:

  Optional character vector of additional cell-type columns to include
  in the composition matrix.

- return:

  Logical. If `TRUE`, saves network PDF plots to `Results/`. Default
  `TRUE`.

- file_name:

  Character. Suffix appended to output file names.

## Value

A named list (one element per factor) of named numeric vectors giving
the enriched cell types and their cumulative NMF edge weights, sorted
descending.
