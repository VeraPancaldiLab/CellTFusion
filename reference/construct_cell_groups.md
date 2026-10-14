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
  padj = FALSE,
  cor_type = "p",
  clustering.method = "ward.D2",
  dendrogram_file = NULL,
  return_dendrogram = FALSE,
  return = FALSE,
  min_cluster_size = 3
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
  significantly correlated with each TF module (adjusted p-value with
  `padj = TRUE`). Default: 0.05.

- padj:

  Logical. If TRUE, the p-values of each TF module are corrected for
  multiple testing across the deconvolution features
  (Benjamini-Hochberg, see
  [`compute.modules.relationship()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.modules.relationship.md))
  before selecting the significant features. Default: FALSE.

- cor_type:

  Correlation used between the TF modules and the deconvolution features
  (passed to
  [`compute.modules.relationship()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.modules.relationship.md)):
  `"p"` (Pearson, default) or `"s"` (Spearman). It is used to select the
  significant features and to cluster them.

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

- min_cluster_size:

  Integer. Minimum number of features of a cell group (minimum cluster
  size of the dynamic tree cut of each dendrogram). Default is 3.

## Value

A named list of 3 elements:

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
cell_groups$Cell_groups[1:3, 1:3]
#>                 Dendrogram_brown_group_1 Dendrogram_brown_group_2
#> SAM7f0d9cc7f001               0.03998411               0.06757056
#> SAM4305ab968b90              -0.13368092              -0.14208692
#> SAMcf018fee2acd              -0.02670350               0.03518956
#>                 Dendrogram_brown_group_3
#> SAM7f0d9cc7f001               0.06172535
#> SAM4305ab968b90              -0.15458310
#> SAMcf018fee2acd              -0.03169306
```
