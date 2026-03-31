# Compute one-step CellTFusion

Compute one-step CellTFusion

## Usage

``` r
CellTFusion(
  raw.counts,
  deconv = NULL,
  normalized = T,
  coldata = NULL,
  batch = F,
  batch_id = NULL,
  deconv_methods = c("Quantiseq", "CBSX", "Epidish", "DeconRNASeq", "DWLS"),
  cbsx.mail = NULL,
  cbsx.token = NULL,
  file_name = NULL,
  task = c("supervised", "unsupervised"),
  contrast = NULL,
  ref_level = NULL,
  TF.collection = "CollecTRI",
  min_targets_size = 3,
  universe = NULL,
  paths = NULL,
  minMod = 3,
  corr_mod = 0.9,
  corr = 0.7,
  corr_type = "spearman",
  cells_extra = NULL,
  pval = 0.05,
  high_corr_groups = 0.8,
  return = T,
  verbose = T
)
```

## Arguments

- raw.counts:

  A matrix of raw gene expression counts (genes as rows, samples as
  columns).

- deconv:

  A data frame with deconvolution features (cell-type proportions as
  columns x samples as rows).

- normalized:

  Logical; if TRUE, normalize raw counts to log-transformed TPM for TF
  computation. For deconvolution they are going to be normalize just as
  TPM. Default is TRUE.

- coldata:

  (Optional) A data frame containing clinical metadata for association
  analysis with TF modules.

- batch:

  Logical; whether batch correction should be applied where supported.
  Default is FALSE.

- batch_id:

  Optional character indicating the column name in coldata containing
  batch identifiers.

- deconv_methods:

  A character vector of deconvolution methods to apply. Default
  includes:
  `c("Quantiseq", "Epidish", "DeconRNASeq", "DWLS", "CibersortX")`.

- cbsx.mail:

  (Optional) Email credential for CIBERSORTx. Required if "CibersortX"
  is among deconv_methods.

- cbsx.token:

  (Optional) Token credential for CIBERSORTx. Required if "CibersortX"
  is among deconv_methods.

- file_name:

  (Optional) Prefix for output files saved in the "Results/" directory.

- task:

  Analysis mode. Choose between `"supervised"` and `"unsupervised"`.

- contrast:

  Optional character indicating the condition column used for supervised
  DEG analysis.

- ref_level:

  Optional character indicating the reference level for supervised DEG
  analysis.

- TF.collection:

  Character. The source of the TF–target network. Options are
  `"CollecTRI"` (default), `"Dorothea"`, or `"ARACNE"`.

  - `"CollecTRI"` and `"Dorothea"` use prebuilt collections from
    OmnipathR.

  - `"ARACNE"` allows user input of a custom network file in a 3-column
    format: `regulator`, `target`, and `mutual information`.

- min_targets_size:

  Integer. Minimum number of target genes per regulon required for TF
  activity inference. Default is 5.

- universe:

  Optional. A user-specified data frame of TF–target interactions. If
  not provided, the function will fetch the relevant network based on
  the `TF.collection` argument.

- paths:

  Optional. A user-specified data frame of pathways gene sets. If not
  provided, the function will fetch the relevant pathways based on
  `PROGENy`.

- minMod:

  Integer; minimum module size for WGCNA module detection.

- corr_mod:

  Numeric; correlation threshold for merging TF modules.

- corr:

  Numeric; correlation threshold used in the deconvolution analysis.

- corr_type:

  Correlation type used in deconvolution analysis. Default is
  `"spearman"`.

- cells_extra:

  A string specifying the cells names to consider and that are not
  including in the nomenclature of multideconv (see R package)

- pval:

  Numeric; p-value threshold for statistical tests (e.g., metadata and
  relationship associations).

- high_corr_groups:

  Numeric; correlation threshold to identify highly similar cell groups.

- return:

  Logical; if TRUE, returns intermediate results from internal
  functions. Default is TRUE.

- verbose:

  Boolen value to whether print or no the function messages

## Value

A list containing:

- Deconvolution:

  A matrix with cell-type proportions (samples as rows, cell types as
  columns).

- TFs_matrix:

  A matrix with TF activity scores (samples as rows, TFs as columns).

- TF_network:

  A list representing the TF module network and related WGCNA output.

- Pathways_scores:

  A matrix of pathway activity scores.

- Processed_deconvolution:

  An object with the processed deconvolution analysis results.

- Cell_groups:

  A matrix of scores representing the cell groups across samples.

## Examples

``` r
if (FALSE) { # \dontrun{
data("raw.counts.tuto")
data("traitdata.tuto")

res <- CellTFusion(
  raw.counts = raw.counts.tuto,
  normalized = TRUE,
  coldata = traitdata.tuto,
  deconv_methods = c("Quantiseq", "DeconRNASeq"),
  file_name = "TestRun",
  min_targets_size = 15,
  minMod = 20,
  corr_mod = 0.25,
  corr = 0.7,
  pval = 0.05,
  high_corr_groups = 0.85
)
} # }
```
