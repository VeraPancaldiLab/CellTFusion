# Identify cell-type niches from NMF latent factors

For each NMF factor, extracts the top-contributing cell group parts
(`<group>_pos` / `<group>_neg`, by basis weight, see
[`extract_contributing_features()`](https://verapancaldilab.github.io/CellTFusion/reference/extract_contributing_features.md)),
maps them to their constituent cell types
([`compute.composition.matrix()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.composition.matrix.md)),
and keeps only cell types enriched relative to the background
composition. Saves a star-network PDF per factor and returns the
weighted cell-type associations.

## Usage

``` r
compute_cells_niches(
  latent_factors,
  dt,
  cell_groups,
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

- cell_groups:

  A list of cell group definitions as returned by
  [`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md).

- enrich_thresh:

  Numeric. Minimum enrichment ratio (see above) for a cell type to be
  retained per factor. Default 1.5.

- quantile_cutoff:

  Numeric between 0 and 1. Quantile threshold for selecting
  top-contributing cell group parts per factor (see
  [`extract_contributing_features()`](https://verapancaldilab.github.io/CellTFusion/reference/extract_contributing_features.md)).
  Default 0.7.

- cells_extra:

  Optional character vector of additional cell type names, added to the
  vocabulary used to recognise cell types in the feature names (see
  [`compute.composition.matrix()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.composition.matrix.md)).

- return:

  Logical. If `TRUE`, saves network PDF plots to
  `Results/Network_<factor>_<file_name>.pdf`
  (`Results/Network_<factor>.pdf` without `file_name`). Default `TRUE`.

- file_name:

  Character. Suffix appended to output file names.

## Value

A named list (one element per factor with at least one enriched cell
type) of named numeric vectors giving the enriched cell types and their
weights, sorted descending.

## Details

For each factor and cell type:

- the weight is the sum of the NMF weights of the selected parts whose
  cell group contains the cell type;

- the enrichment is the fraction of selected parts containing the cell
  type divided by the fraction of all cell groups containing it. It is a
  ratio, not a statistical test: with few selected parts, a cell type
  present in only one of them can pass the threshold.

The direction is not taken into account: a selected part can mean more
or fewer cells of its group (depending on `_pos` / `_neg` and on whether
the features of the group are positively or negatively related to its TF
module), and both are added together. A cell type of a niche is
therefore associated with the factor, without saying whether the factor
has more or fewer of it.

## Examples

``` r
data("network.tuto")
data("deconv_subgroups.tuto")
cell_groups <- construct_cell_groups(network.tuto, deconv_subgroups.tuto)
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Dendrogram plots saved to: Results/Dendrogram_color_clusters.pdf
#> 44 of 897 features were not assigned to any cluster by the tree cut and are not part of any cell group (pink: 2, yellow: 8, red: 2, blue: 6, brown: 4, turquoise: 12, black: 4, green: 6)
latent <- compute.latent_factors(cell_groups$Cell_groups, return = FALSE)
#> Estimating optimal rank via elbow on reconstruction MSE (ranks 2:8)...
#> MSE by rank: k2=0.00129  k3=0.00114  k4=0.00105  k5=0.000978  k6=0.000933  k7=0.000893  k8=0.000857
#> Selected rank (elbow): 3
#> Running RcppML NMF with rank = 3...
niches <- compute_cells_niches(latent, deconv_subgroups.tuto, cell_groups, return = FALSE)
#> Cell types not present in any cell group (left out of the composition matrix): T.cells.gamma.delta
niches
#> $Factor1
#>         Cancer      NKT.cells         Plasma     NK.resting T.cells.helper 
#>     0.23016502     0.08820450     0.07970784     0.05736662     0.05649637 
#> 
#> $Factor2
#>            Cancer Macrophages.cells    T.cells.helper       Eosinophils 
#>        0.20564591        0.09154775        0.06090095        0.03824595 
#> 
#> $Factor3
#> CD4.memory.activated       B.memory.cells             Myocytes 
#>           0.06358580           0.06054229           0.02153581 
#> 
```
