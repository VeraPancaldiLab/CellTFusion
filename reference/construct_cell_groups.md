# Construct cell groups based on TF networks and deconvolution

Identifies and projects cell groups using module relationships derived
from TF networks and deconvolution outputs.

## Usage

``` r
construct_cell_groups(
  network,
  dt,
  batch = NULL,
  pval = 0.05,
  clustering.method = "ward.D2",
  n_perm = 999,
  dendrogram_file = NULL,
  return_dendrogram = FALSE
)
```

## Arguments

- network:

  A TF module network as returned by
  [`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md).

- dt:

  Deconvolution subgroups as returned by
  [`multideconv::compute.deconvolution.analysis()`](https://verapancaldilab.github.io/multideconv/reference/compute.deconvolution.analysis.html).

- batch:

  Optional vector indicating batch assignment for samples.

- pval:

  Numeric. P-value threshold applied both to filter TF
  module-deconvolution feature correlations and as the significance
  cutoff for the CCA permutation test. Default: 0.05.

- clustering.method:

  Clustering method for hierarchical clustering. Default: "ward.D2".

- n_perm:

  Integer. Number of permutations for the CCA significance test per cell
  group. Higher values give more precise p-values but increase runtime.
  Default: 999.

- dendrogram_file:

  Optional character. Suffix of the dendrogram PDF file name (see
  `return_dendrogram`).

- return_dendrogram:

  Logical. If TRUE, saves a PDF of the colored cell-group dendrograms to
  `Results/Dendrogram_color_clusters_<dendrogram_file>.pdf` (only when
  `dendrogram_file` is set). Default FALSE.

## Value

A named list of 3 elements:

- Cell_groups:

  A data frame with the projected cell group scores (samples x groups).

- Composition:

  A named list where each element is a character vector of the
  deconvolution features in each group.

- Weights:

  A list of CCA projection parameters (`xcoef`, `train_means`,
  `train_sds`) for each group.
