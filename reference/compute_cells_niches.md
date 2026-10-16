# Identify cell niches from NMF latent factors

For each NMF factor, selects the cell groups that contribute most to it
(the `top_n` highest weights of the factor in `W`, see
[`extract_contributing_features()`](https://verapancaldilab.github.io/CellTFusion/reference/extract_contributing_features.md))
and describes them by their deconvolution features and CCA coefficients.
Deconvolution features are kept as they are (subgroups are not expanded
into their members).

## Usage

``` r
compute_cells_niches(
  latent_factors,
  cell_groups,
  top_n = 5,
  return = TRUE,
  file_name = NULL
)
```

## Arguments

- latent_factors:

  A list returned by
  [`compute.latent_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.latent_factors.md)
  (its `W` is used).

- cell_groups:

  The cell groups, as returned by
  [`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md)
  (its `Weights` are used).

- top_n:

  Integer. Number of cell groups selected per factor. Default 5.

- return:

  Logical. If `TRUE` (default), saves the networks in
  `Results/Cells_niches_<file_name>.pdf` (`Results/Cells_niches.pdf`
  without `file_name`), one page per factor, and a heatmap of the NMF
  weights (`W`) of all the cell groups (rows: cell groups, ordered by
  the factor where their weight is highest and by decreasing weight;
  columns: factors; outlined: the top cell groups of each factor) in
  `Results/Cells_niches_heatmap_<file_name>.pdf`.

- file_name:

  Character. Suffix appended to the output file name.

## Value

A data frame with one row per factor, selected cell group and feature:
`factor`, `cell_group`, `module` (TF module of the cell group), `weight`
(NMF weight of the cell group in the factor), `feature` and
`coefficient` (CCA coefficient of the feature in the cell group score).
Rows are ordered by factor and by decreasing weight.

## Details

With `return = TRUE`, a tree network of each factor is saved, one page
per factor: the factor on the left, its top cell groups in the middle
(filled with the colour of their TF module) and the deconvolution
features of each cell group listed from top to bottom on the right, next
to their cell group. The edge between a cell group and a feature is the
CCA coefficient of the feature in the cell group score (`xcoef` in
`Weights`, on the z-scored features): red if positive (the score
increases with the feature), blue if negative, and wider for a larger
absolute coefficient (coefficients close to zero are hairlines).

The factors describe the high-activity side of the TF modules of their
cell groups (the NMF input is the positive part of the cell group
scores). A red edge therefore means that the feature is higher in
samples where the cell group is active, a blue edge that it is lower.
The coefficients of features correlated with each other (e.g. several
estimates of the same cell type) are estimated jointly, so a single
coefficient can be small or have the opposite sign of the others of its
cell group even when the feature itself is related to the TF module like
them: read the sign of a single edge together with the rest of its cell
group.

## Examples

``` r
data("network.tuto")
data("deconv_subgroups.tuto")
cell_groups <- construct_cell_groups(network.tuto, deconv_subgroups.tuto, return = FALSE)
#> 44 of 897 features were not assigned to any cluster by the tree cut and are not part of any cell group (pink: 2, yellow: 8, red: 2, blue: 6, brown: 4, turquoise: 12, black: 4, green: 6)
latent <- compute.latent_factors(cell_groups$Cell_groups, return = FALSE)
#> Selecting the rank with bootstrap NMF (100 runs per rank, ranks 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15)...
#> Selected rank: 5 (largest rank with all factors robust and none redundant)
#> Running RcppML NMF with rank = 5...
niches <- compute_cells_niches(latent, cell_groups, return = FALSE)
head(niches)
#>    factor             cell_group module    weight
#> 1 Factor1 Dendrogram_red_group_3    red 0.0222516
#> 2 Factor1 Dendrogram_red_group_3    red 0.0222516
#> 3 Factor1 Dendrogram_red_group_3    red 0.0222516
#> 4 Factor1 Dendrogram_red_group_3    red 0.0222516
#> 5 Factor1 Dendrogram_red_group_3    red 0.0222516
#> 6 Factor1 Dendrogram_red_group_3    red 0.0222516
#>                                           feature   coefficient
#> 1                              B.cells_Subgroup.5 -0.0219763701
#> 2                          Neutrophils_Subgroup.3 -0.0091728791
#> 3                     DeconRNASeq_LM22_NK.resting -0.0160350224
#> 4 DeconRNASeq_CBSX.NSCLC.PBMCs.scRNAseq_NKT.cells  0.0003870828
#> 5 DeconRNASeq_CBSX.NSCLC.PBMCs.scRNAseq_CD4.cells  0.0014083163
#> 6                  DeconRNASeq_BPRNACan_CD8.cells -0.0396184688
```
