# Compute cell group scores from deconvolution and TF module network

This function identifies cell groups based on dendrogram cuts, computes
composite scores for each group using deconvolution features and TF
module networks, and optionally exports the results.

## Usage

``` r
cell.groups.computation(
  deconv,
  cell_dendrograms,
  network,
  batch = NULL,
  return = TRUE,
  file_name = NULL,
  min_cluster_size = 3
)
```

## Arguments

- deconv:

  A data frame with deconvolution features (samples as rows, cell-type
  features as columns). This is usually the first element returned by
  [`multideconv::compute.deconvolution.analysis()`](https://verapancaldilab.github.io/multideconv/reference/compute.deconvolution.analysis.html).
  It must contain all the features of `cell_dendrograms`.

- cell_dendrograms:

  A named list of dendrogram (`hclust`) objects, each corresponding to a
  TF module, typically returned by
  [`identify.cell.groups()`](https://verapancaldilab.github.io/CellTFusion/reference/identify.cell.groups.md).

- network:

  A list containing network information of transcription factor (TF)
  modules, as obtained from
  [`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md)
  (the module scores, element `TFs module matrix`, are used).

- batch:

  Optional vector indicating batch assignment for samples.

- return:

  Logical; if TRUE (default), saves a PDF of the colored cell-group
  dendrograms (`Results/Dendrogram_color_clusters_<file_name>.pdf`) and
  CSV files with cell group compositions and scores to the "Results/"
  folder (the folder is created if needed).

- file_name:

  Optional character. Suffix of the output file names (see `return`).

- min_cluster_size:

  Integer. Minimum number of features of a cell group (minimum cluster
  size of the dynamic tree cut of each dendrogram). Default is 3.

## Value

An unnamed list of three elements:

- `[[1]]`: A data frame with the composite scores of all identified cell
  groups across samples. A positive score means that the sample is on
  the active side of the TF module of that cell group.

- `[[2]]`: A named list of vectors indicating the composition (original
  features) of each cell group.

- `[[3]]`: A named list of CCA projection parameters (`xcoef`,
  `train_means`, `train_sds`; only `xcoef` with `batch`) for each cell
  group.

Cell groups are named `Dendrogram_<module>_group_<n>`, where `n` is the
cluster number of the dendrogram of that module (the number shown in the
dendrogram PDF); the three elements use the same names. An error is
raised if no cluster is found in any module. Features that the tree cut
does not assign to any cluster are not part of any cell group: a message
reports how many there are per module. If `return=TRUE`, two CSV files
will be created:

- `Results/Cell.groups.composition_<file_name>.csv`: A table showing the
  composition of each cell group.

- `Results/Cell.groups.scores_<file_name>.csv`: A matrix of cell group
  scores across samples.

Without `file_name`, they are named `Cell.groups.composition.csv` and
`Cell.groups.scores.csv`.

## Examples

``` r
if (FALSE) { # \dontrun{
deconv_results <- multideconv::compute.deconvolution.analysis(...)
tf_network <- compute.WTCNA(...)
dendrograms <- identify.cell.groups(...)

cell.groups <- cell.groups.computation(
  deconv = deconv_results[[1]],
  cell_dendrograms = dendrograms,
  network = tf_network,
  return = TRUE
)
} # }
```
