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
or fewer cells of its group (depending on `_pos` / `_neg` and on the
relation of the group to its TF module, see `Module_relation` in
[`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md)),
and both are added together. A cell type of a niche is therefore
associated with the factor, without saying whether the factor has more
or fewer of it.

## Examples

``` r
data("network.tuto")
data("deconv_subgroups.tuto")
cell_groups <- construct_cell_groups(network.tuto, deconv_subgroups.tuto)
#> 20 of 764 features were not assigned to any cluster by the tree cut and are not part of any cell group (blue: 6, pink: 2, red: 2, yellow: 8, black: 2)
latent <- compute.latent_factors(cell_groups$Cell_groups, return = FALSE)
#> Estimating optimal rank via elbow on reconstruction MSE (ranks 2:8)...
#> MSE by rank: k2=0.00127  k3=0.000964  k4=0.000837  k5=0.000754  k6=0.000708  k7=0.000666  k8=0.000635
#> Selected rank (elbow): 3
#> Running RcppML NMF with rank = 3...
niches <- compute_cells_niches(latent, deconv_subgroups.tuto, cell_groups, return = FALSE)
niches
#> $Factor1
#>                  CAF               Plasma           NK.resting 
#>           0.13981351           0.08590253           0.07909920 
#> Mast.activated.cells          Eosinophils 
#>           0.04980160           0.03468501 
#> 
#> $Factor2
#> Mast.activated.cells 
#>           0.05223202 
#> 
#> $Factor3
#>             B.cells T.cells.gamma.delta 
#>          0.35521646          0.02325468 
#> 
```
