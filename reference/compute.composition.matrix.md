# Compute a cell-type composition matrix from deconvolution subgroups

Builds a binary presence matrix indicating which original cell types are
present in each cell group. Each feature of a cell group is opened down
to its original deconvolution features (a subgroup is replaced by its
members), and the cell type of each feature is read from its name with
the cell type vocabulary of multideconv
([`multideconv::get_cell_type_nomenclature()`](https://verapancaldilab.github.io/multideconv/reference/get_cell_type_nomenclature.html)).
Cell types are therefore named as in that vocabulary: for example the
LM22 features `Plasma.cells` and `T.cells.CD4.memory.activated` are
counted as `Plasma` and `CD4.memory.activated`, as in the names of the
deconvolution subgroups.

## Usage

``` r
compute.composition.matrix(dt, cell_groups, cells_extra = NULL)
```

## Arguments

- dt:

  Deconvolution subgroups as returned by
  [`multideconv::compute.deconvolution.analysis()`](https://verapancaldilab.github.io/multideconv/reference/compute.deconvolution.analysis.html)
  (the elements `Deconvolution matrix` and
  `Deconvolution subgroups composition` are used).

- cell_groups:

  The cell groups returned by
  [`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md):
  its element `Composition` (the features of each cell group) is used.
  An unnamed list whose second element holds the composition is also
  accepted.

- cells_extra:

  Optional character vector of additional cell type names, added to the
  vocabulary used to recognise cell types in the feature names.

## Value

A binary matrix (data.frame) where rows are cell groups and columns are
cell types (1 = present, 0 = absent). Rows are named as the cell groups,
as the features of
[`compute.latent_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.latent_factors.md).
Cell types not present in any cell group are dropped, and cell groups
whose cell types are not recognised (their row is all zeros) are listed,
each with a message.

## Examples

``` r
data("network.tuto")
data("deconv_subgroups.tuto")
cell_groups <- construct_cell_groups(network.tuto, deconv_subgroups.tuto)
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
#> ℹ Please use `linewidth` instead.
#> ℹ The deprecated feature was likely used in the dendextend package.
#>   Please report the issue at <https://github.com/talgalili/dendextend/issues>.
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Warning: Length of color vector was shorter than the number of clusters - color vector was recycled
#> Dendrogram plots saved to: Results/Dendrogram_color_clusters.pdf
#> 44 of 897 features were not assigned to any cluster by the tree cut and are not part of any cell group (pink: 2, yellow: 8, red: 2, blue: 6, brown: 4, turquoise: 12, black: 4, green: 6)
composition <- compute.composition.matrix(deconv_subgroups.tuto, cell_groups)
#> Cell types not present in any cell group (left out of the composition matrix): T.cells.gamma.delta
composition[1:4, 1:5]
#>                         B.cells B.naive.cells B.memory.cells Macrophages.cells
#> Dendrogram_pink_group_1       0             0              0                 0
#> Dendrogram_pink_group_2       1             0              0                 0
#> Dendrogram_pink_group_3       1             1              0                 0
#> Dendrogram_pink_group_4       1             0              1                 0
#>                         Macrophages.M0
#> Dendrogram_pink_group_1              0
#> Dendrogram_pink_group_2              0
#> Dendrogram_pink_group_3              1
#> Dendrogram_pink_group_4              0
```
