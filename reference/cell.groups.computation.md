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

  A data frame with deconvolution features (typically a cell-type or
  cluster x sample matrix). This is usually the first element returned
  by
  [`multideconv::compute.deconvolution.analysis()`](https://rdrr.io/pkg/multideconv/man/compute.deconvolution.analysis.html).

- cell.dendrograms:

  A named list of dendrogram objects, each corresponding to a TF module,
  typically returned by
  [`identify.cell.groups()`](https://verapancaldilab.github.io/CellTFusion/reference/identify.cell.groups.md).

- tfs.module.network:

  A list containing network information of transcription factor (TF)
  modules, as obtained from
  [`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md).
  It should contain at least one element with TF module membership or
  connectivity.

- batch:

  Optional vector indicating batch assignment for samples.

- return:

  Logical; if TRUE (default), writes CSV files with cell group
  compositions and scores to the "Results/" folder.

- pval:

  Numeric. P-value threshold for statistical tests. Default is 0.05.

- n_perm:

  Integer. Number of permutations for significance testing. Default is
  999.

- dendrogram_file:

  Optional character. File path to save dendrogram plot.

- return_dendrogram:

  Logical. If TRUE, includes the dendrogram in the returned list.
  Default FALSE.

## Value

A list of three elements:

- scores:

  A data frame with the composite scores of all identified cell groups
  across samples.

- composition:

  A list of vectors indicating the composition (original features) of
  each cell group.

- loadings:

  A list of loadings (feature contributions) for each cell group.

If `return=TRUE`, two CSV files will be created:

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
