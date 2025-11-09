# Construct cell groups based on TF networks and deconvolution

Identifies and projects cell groups using module relationships derived
from TF networks and deconvolution outputs. If a binary trait is
specified, the function splits the data and constructs cell groups for
both classes (supervised analysis).

## Usage

``` r
construct_cell_groups(
  counts,
  tfs,
  deconv,
  network,
  dt,
  clinical,
  pval = 0.05,
  high_corr_groups = 0.8,
  clustering.method = "ward.D2",
  trait = NULL,
  positive = NULL,
  TF.collection = "CollecTRI",
  min_targets_size = 10,
  tfs.pruned = FALSE,
  universe = NULL
)
```

## Arguments

- counts:

  A matrix of gene expression counts (genes x samples).

- tfs:

  A list or matrix of transcription factors used in the analysis.

- deconv:

  A matrix of cell type proportions from deconvolution (samples x cell
  types).

- network:

  A list containing TF networks for cell types.

- dt:

  A list containing deconvolution subgroup structures.

- clinical:

  A data frame with clinical metadata, including the trait of interest.

- pval:

  P-value threshold used to filter module relationships. Default: 0.05.

- high_corr_groups:

  Correlation threshold to merge or remove redundant cell groups.
  Default: 0.9.

- clustering.method:

  Clustering method for hierarchical clustering. Default: "ward.D2".

- trait:

  Optional character: column name in `clinical` for trait to split by
  and do a supervised cell group analysis (see paper for more info). If
  no provided, analysis will be unsupervised.

- positive:

  Optional value defining the positive class of the `trait`.

- TF.collection:

  Character. The source of the TF–target network. Options are
  `"CollecTRI"` (default), `"Dorothea"`, or `"ARACNE"`. Only needed when
  supervised analysis will be performed, if not, it will be ignored.

  - `"CollecTRI"` and `"Dorothea"` use prebuilt collections from
    OmnipathR.

  - `"ARACNE"` allows user input of a custom network file in a 3-column
    format: `regulator`, `target`, and `mutual information`.

- min_targets_size:

  Integer. Minimum number of target genes per regulon required for TF
  activity inference. Default is 5. Only needed when supervised analysis
  will be performed, if not, it will be ignored.

- tfs.pruned:

  Logical. Whether to prune TF regulons to limit the number of target
  genes, which helps reduce bias introduced by TFs with large regulons.
  If `TRUE`, the user will be prompted to input a maximum size for
  regulons. Default is `FALSE`. Only needed when supervised analysis
  will be performed, if not, it will be ignored.

- universe:

  Optional. A user-specified data frame of TF–target interactions. If
  not provided, the function will fetch the relevant network based on
  the `TF.collection` argument. Only needed when supervised analysis
  will be performed, if not, it will be ignored.

## Value

A list of 3 elements:

- scores:

  A data frame or matrix with the projected cell group scores (samples x
  groups).

- composition:

  A named list where each element is a character vector of original cell
  types per group.

- loadings:

  A list of numeric vectors indicating the loadings (feature
  contributions) for each group.
