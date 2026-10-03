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
  pval = 0.05,
  n_perm = 999,
  dendrogram_file = NULL,
  return_dendrogram = FALSE
)
```

## Arguments

- deconvolution:

  A data frame with deconvolution features (samples as rows, cell-type
  features as columns). This is usually the first element returned by
  [`multideconv::compute.deconvolution.analysis()`](https://verapancaldilab.github.io/multideconv/reference/compute.deconvolution.analysis.html).

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
  compositions and scores to the "Results/" folder.

- pval:

  Numeric. Significance threshold for the CCA permutation test of each
  cell group. Default is 0.05.

- n_perm:

  Integer. Number of permutations for significance testing. Default is
  999.

- dendrogram_file:

  Optional character. Suffix of the dendrogram PDF file name (see
  `return_dendrogram`).

- return_dendrogram:

  Logical. If TRUE, saves a PDF of the colored cell-group dendrograms to
  `Results/Dendrogram_color_clusters_<dendrogram_file>.pdf` (only when
  `dendrogram_file` is set). Default FALSE.

## Value

An unnamed list of three elements:

- `[[1]]`: A data frame with the composite scores of all identified cell
  groups across samples.

- `[[2]]`: A list of vectors indicating the composition (original
  features) of each cell group.

- `[[3]]`: A list of CCA projection parameters (`xcoef`, `train_means`,
  `train_sds`) for each cell group.

Cell groups made of a single feature, or that do not pass the CCA
permutation test, are discarded; an error is raised if no cell group
remains. If `return=TRUE`, two CSV files will be created:

- `Results/Cell.groups.composition.csv`: A table showing the composition
  of each cell group.

- `Results/Cell.groups.scores.csv`: A matrix of cell group scores across
  samples.

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
