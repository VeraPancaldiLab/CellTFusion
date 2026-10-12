# Compute cell group scores from deconvolution and TF module network

This function identifies cell groups based on dendrogram cuts, computes
composite scores for each group using deconvolution features and TF
module networks, and optionally exports the results.

## Usage

``` r
cell.groups.computation(
  deconvolution,
  cell.dendrograms,
  tfs.module.network,
  batch = NULL,
  return = T,
  dendrogram_file = NULL,
  return_dendrogram = FALSE,
  cells_extra = NULL,
  min_cluster_size = 3,
  deep_split = TRUE
)
```

## Arguments

- deconvolution:

  A data frame with deconvolution features (samples as rows, cell-type
  features as columns). This is usually the first element returned by
  [`multideconv::compute.deconvolution.analysis()`](https://verapancaldilab.github.io/multideconv/reference/compute.deconvolution.analysis.html).
  It must contain all the features of `cell.dendrograms`.

- cell.dendrograms:

  A named list of dendrogram (`hclust`) objects, each corresponding to a
  TF module, typically returned by
  [`identify.cell.groups()`](https://verapancaldilab.github.io/CellTFusion/reference/identify.cell.groups.md).

- tfs.module.network:

  A list containing network information of transcription factor (TF)
  modules, as obtained from
  [`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md)
  (the elements `TFs module matrix`, `TFs per module` and `TFs_matrix`
  are used).

- batch:

  Optional vector indicating batch assignment for samples.

- return:

  Logical; if TRUE (default), writes CSV files with cell group
  compositions and scores to the "Results/" folder (named with
  `dendrogram_file`; the folder is created if needed).

- dendrogram_file:

  Optional character. Suffix of the output file names: the dendrogram
  PDF (see `return_dendrogram`) and the CSV files (see `return`).

- return_dendrogram:

  Logical. If TRUE, saves a PDF of the colored cell-group dendrograms to
  `Results/Dendrogram_color_clusters_<dendrogram_file>.pdf` (only when
  `dendrogram_file` is set). Default FALSE.

- cells_extra:

  Optional character vector with the cell names that are not in the
  nomenclature of multideconv, used to list the cell types of each cell
  group.

- min_cluster_size:

  Integer. Minimum number of features of a cell group (minimum cluster
  size of the dynamic tree cut of each dendrogram). Default is 3.

- deep_split:

  Logical. If TRUE (default), the dynamic tree cut splits the clusters
  more finely.

## Value

An unnamed list of four elements:

- `[[1]]`: A data frame with the composite scores of all identified cell
  groups across samples. A positive score means that the sample is on
  the active side of the TF module of that cell group.

- `[[2]]`: A named list of vectors indicating the composition (original
  features) of each cell group.

- `[[3]]`: A named list of CCA projection parameters (`xcoef`,
  `train_means`, `train_sds`) for each cell group.

- `[[4]]`: A data frame with the relation of each cell group to its TF
  module (see
  [`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md)).

Cell groups are named `Dendrogram_<module>_group_<n>`, where `n` is the
cluster number of the dendrogram of that module (the number shown in the
dendrogram PDF); the first three elements use the same names. An error
is raised if no cluster is found in any module. Features that the tree
cut does not assign to any cluster are not part of any cell group: a
message reports how many there are per module. If `return=TRUE`, two CSV
files will be created:

- `Results/Cell.groups.composition_<dendrogram_file>.csv`: A table
  showing the composition of each cell group, its TF module, its
  relation to it (positive or negative, and the correlation) and its
  cell types.

- `Results/Cell.groups.scores_<dendrogram_file>.csv`: A matrix of cell
  group scores across samples.

Without `dendrogram_file`, they are named `Cell.groups.composition.csv`
and `Cell.groups.scores.csv`.

## Examples

``` r
if (FALSE) { # \dontrun{
deconv_results <- multideconv::compute.deconvolution.analysis(...)
tf_network <- compute.WTCNA(...)
dendrograms <- identify.cell.groups(...)

cell.groups <- cell.groups.computation(
  deconvolution = deconv_results[[1]],
  cell.dendrograms = dendrograms,
  tfs.module.network = tf_network,
  return = TRUE
)
} # }
```
