# Construct cell groups based on TF networks and deconvolution

A cell group is a set of deconvolution features with a similar
biological activity: features related to the same TF module to the same
degree (similar correlation with the module score). The features of a
cell group do not need to vary together across samples. The cell groups
are built in three steps:

1.  each TF module is correlated with every deconvolution feature
    ([`compute.modules.relationship()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.modules.relationship.md));

2.  for each module, the significantly correlated features are clustered
    on their correlation with the module
    ([`identify.cell.groups()`](https://verapancaldilab.github.io/CellTFusion/reference/identify.cell.groups.md));

3.  each dendrogram is cut into cell groups, and each cell group gets
    one score per sample from a CCA with the TFs of its module
    ([`cell.groups.computation()`](https://verapancaldilab.github.io/CellTFusion/reference/cell.groups.computation.md)).

Samples are matched by name: the module scores and TF activities of
`network` are put in the sample order of the deconvolution matrix, which
may contain only part of the samples of the network.

## Usage

``` r
construct_cell_groups(
  network,
  dt,
  batch = NULL,
  pval = 0.05,
  clustering.method = "ward.D2",
  dendrogram_file = NULL,
  return_dendrogram = FALSE,
  return = FALSE,
  cells_extra = NULL,
  min_cluster_size = 3,
  deep_split = TRUE
)
```

## Arguments

- network:

  A TF module network as returned by
  [`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md)
  (its module scores, TF activities and TFs per module are used).

- dt:

  Deconvolution subgroups as returned by
  [`multideconv::compute.deconvolution.analysis()`](https://verapancaldilab.github.io/multideconv/reference/compute.deconvolution.analysis.html);
  its element `Deconvolution matrix` (samples x features) is used. Every
  sample must be in `network`.

- batch:

  Optional vector with the cohort of each sample: in the row order of
  the deconvolution matrix, or named by sample. Required when `network`
  was built with `compute.WTCNA(batch = TRUE)`. With two or more
  cohorts, correlations and scores are corrected for cohort.

- pval:

  Numeric. P-value threshold used to select the deconvolution features
  significantly correlated with each TF module. Default: 0.05.

- clustering.method:

  Linkage method of the hierarchical clustering of the deconvolution
  features of each module on their correlation with it (see
  [`identify.cell.groups()`](https://verapancaldilab.github.io/CellTFusion/reference/identify.cell.groups.md)).
  Default: "ward.D2".

- dendrogram_file:

  Optional character. Suffix of the output file names (see
  `return_dendrogram` and `return`).

- return_dendrogram:

  Logical. If TRUE, saves a PDF of the colored cell-group dendrograms to
  `Results/Dendrogram_color_clusters_<dendrogram_file>.pdf` (only when
  `dendrogram_file` is set). Default FALSE.

- return:

  Logical. If TRUE, saves the composition and the scores of the cell
  groups to `Results/Cell.groups.composition_<dendrogram_file>.csv` and
  `Results/Cell.groups.scores_<dendrogram_file>.csv`. Default FALSE.

- cells_extra:

  Optional character vector with the cell names that are not in the
  nomenclature of multideconv (the same as in
  [`multideconv::compute.deconvolution.analysis()`](https://verapancaldilab.github.io/multideconv/reference/compute.deconvolution.analysis.html)),
  used to list the cell types of each cell group in `Module_relation`.

- min_cluster_size:

  Integer. Minimum number of features of a cell group (minimum cluster
  size of the dynamic tree cut of each dendrogram). Default is 3.

- deep_split:

  Logical. If TRUE (default), the dynamic tree cut splits the clusters
  more finely.

## Value

A named list of 4 elements:

- Cell_groups:

  A data frame with the projected cell group scores (samples x groups).
  A positive score means that the sample is on the active side of the TF
  module of that cell group, a negative score that it is on the inactive
  side.

- Composition:

  A named list where each element is a character vector of the
  deconvolution features in each group.

- Weights:

  A named list of CCA projection parameters (`xcoef`, `train_means`,
  `train_sds`) for each group. With `batch`, `train_means` is `NULL`:
  test data are centred on their own means.

- Module_relation:

  A data frame with one row per cell group: `cell_group`, `module` (its
  TF module), `relation` (`"positive"` or `"negative"`), `correlation`
  (mean correlation of its features with the module score), `n_features`
  and `cell_types` (the cell types of its features, the most frequent
  first). It translates a score into cell abundance: for a positively
  related group, a positive score means more of its features; for a
  negatively related group, a positive score means fewer of them.

Cell groups are named `Dendrogram_<module>_group_<n>`, where `n` is the
cluster number shown in the dendrogram PDF of that module. Features that
the tree cut does not assign to any cluster are not part of any cell
group (a message reports how many).

## Examples

``` r
data("network.tuto")
data("deconv_subgroups.tuto")
cell_groups <- construct_cell_groups(network.tuto, deconv_subgroups.tuto)
#> 20 of 764 features were not assigned to any cluster by the tree cut and are not part of any cell group (blue: 6, pink: 2, red: 2, yellow: 8, black: 2)
head(cell_groups$Module_relation)
#>                 cell_group module relation correlation n_features
#> 1 Dendrogram_brown_group_1  brown negative  -0.4082087          8
#> 2 Dendrogram_brown_group_2  brown negative  -0.3277192          5
#> 3 Dendrogram_brown_group_3  brown negative  -0.2895496          5
#> 4 Dendrogram_brown_group_4  brown negative  -0.2522286          4
#> 5 Dendrogram_brown_group_5  brown negative  -0.1882202         16
#> 6 Dendrogram_brown_group_6  brown positive   0.1479860          3
#>                                                                                                                                                                                              cell_types
#> 1                                                                                                                         Cancer, B.cells, CD8.cells, Eosinophils, Macrophages.M2, NK.cells, NK.resting
#> 2                                                                                                                                                CD4.cells, Cancer, NK.cells, NKT.cells, T.cells.helper
#> 3                                                                                                                                      Endothelial, Macrophages.M1, Macrophages.M2, Neutrophils, Plasma
#> 4                                                                                                                                                           CD4.cells, NK.cells, NKT.cells, Neutrophils
#> 5 Cancer, CD4.cells, CD4.memory.resting, CD4.naive, CD8.cells, Dendritic.activated.cells, Endothelial, Macrophages.M0, Macrophages.cells, Mast.resting.cells, NK.cells, NK.resting, Neutrophils, Plasma
#> 6                                                                                                                                                                                   Monocytes, NK.cells
```
