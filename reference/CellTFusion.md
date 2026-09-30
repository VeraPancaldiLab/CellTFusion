# Compute one-step CellTFusion

Compute one-step CellTFusion

## Usage

``` r
CellTFusion(
  raw.counts,
  deconv = NULL,
  dt = NULL,
  tfs = NULL,
  pathways = NULL,
  normalized = T,
  coldata = NULL,
  batch = F,
  batch_id = NULL,
  deconv_methods = c("Quantiseq", "Epidish", "DeconRNASeq", "DWLS"),
  cbsx.mail = NULL,
  cbsx.token = NULL,
  file_name = NULL,
  TF.collection = "CollecTRI",
  min_targets_size = 3,
  universe = NULL,
  paths = NULL,
  gene_sets = NULL,
  minMod = 3,
  corr_mod = 0.9,
  corr = 0.7,
  corr_type = "spearman",
  cells_extra = NULL,
  pval = 0.05,
  enrich_thresh = 1.5,
  quantile_cutoff = 0.7,
  cancer_type = NULL,
  return = T,
  verbose = T
)
```

## Arguments

- raw.counts:

  A matrix of raw gene expression counts (genes as rows, samples as
  columns). Always required: even when `dt`, `tfs` and `pathways` are
  all supplied precomputed, the (normalized) expression matrix is still
  needed for the downstream GSEA-based TME state characterization step.

- deconv:

  A data frame with deconvolution features (cell-type proportions as
  columns x samples as rows). Ignored if `dt` is supplied.

- dt:

  (Optional) A precomputed cell-subgroup object, typically the output of
  [`multideconv::compute.deconvolution.analysis()`](https://verapancaldilab.github.io/multideconv/reference/compute.deconvolution.analysis.html)
  (or the `Processed_deconvolution` element returned by a previous
  `CellTFusion()` run). If supplied, cell-type deconvolution and the
  deconvolution analysis step are both skipped and the pipeline proceeds
  straight to cell group construction using this object.

- tfs:

  (Optional) A precomputed TF activity matrix (samples as rows, TFs as
  columns), typically the `TFs_matrix` element returned by a previous
  `CellTFusion()` run or by
  [`compute.TFs.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.TFs.activity.md).
  If supplied, TF activity inference is skipped.

- pathways:

  (Optional) A precomputed pathway activity matrix, typically the
  `Pathways_scores` element returned by a previous `CellTFusion()` run
  or by
  [`compute.pathway.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.pathway.activity.md).
  If supplied, pathway activity inference is skipped.

- normalized:

  Logical; if TRUE, normalize raw counts to log-transformed TPM for TF
  computation. For deconvolution they are going to be normalize just as
  TPM. Default is TRUE.

- coldata:

  (Optional) A data frame of sample metadata (samples as rows). Only
  used when `batch = TRUE`, to read the batch column given by
  `batch_id`.

- batch:

  Logical; whether batch correction should be applied where supported.
  Default is FALSE.

- batch_id:

  Optional character indicating the column name in coldata containing
  batch identifiers.

- deconv_methods:

  A character vector of deconvolution methods to apply. Default is
  `c("Quantiseq", "Epidish", "DeconRNASeq", "DWLS")`. Add `"CBSX"` to
  also run CIBERSORTx (requires `cbsx.mail` and `cbsx.token`).

- cbsx.mail:

  (Optional) Email credential for CIBERSORTx. Required if "CBSX" is
  among deconv_methods.

- cbsx.token:

  (Optional) Token credential for CIBERSORTx. Required if "CBSX" is
  among deconv_methods.

- file_name:

  (Optional) Prefix for output files saved in the "Results/" directory.

- TF.collection:

  Character. The source of the TF-target network. Options are
  `"CollecTRI"` (default), `"Dorothea"`, or `"ARACNE"`.

  - `"CollecTRI"` and `"Dorothea"` use prebuilt collections from
    OmnipathR.

  - `"ARACNE"` reads a network file (tab-separated with `Regulator` and
    `Target` columns) from
    `input/ARACNE/<cancer_type>/network/network.txt`.

- min_targets_size:

  Integer. Minimum number of target genes per regulon, passed to
  [`compute.TFs.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.TFs.activity.md).
  Default is 3.

- universe:

  Optional. A user-specified data frame of TF-target interactions. If
  not provided, the function will fetch the relevant network based on
  the `TF.collection` argument.

- paths:

  Optional. A user-specified data frame of pathways gene sets. If not
  provided, the function will fetch the relevant pathways based on
  `PROGENy`.

- gene_sets:

  Optional. A named list of custom gene sets (character vectors of gene
  symbols) passed to
  [`compute.pathway.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.pathway.activity.md)'s
  `gene_sets` argument for GSVA-based scoring. If `NULL`, only PROGENy
  is used.

- minMod:

  Integer; minimum module size for WGCNA module detection. Default is 3.

- corr_mod:

  Numeric; correlation threshold for merging TF modules. Default is 0.9.

- corr:

  Numeric; correlation threshold used in the deconvolution analysis.

- corr_type:

  Correlation type used in deconvolution analysis. Default is
  `"spearman"`.

- cells_extra:

  A string specifying the cells names to consider and that are not
  including in the nomenclature of multideconv (see R package)

- pval:

  Numeric; p-value threshold used when building cell groups (TF
  module-deconvolution correlations and the CCA permutation test).

- enrich_thresh:

  Numeric. Minimum enrichment ratio (foreground/background cell-type
  frequency) required to include a cell type in a latent factor's niche.
  Default is 1.5.

- quantile_cutoff:

  Numeric between 0 and 1. Quantile threshold for selecting
  top-contributing cell groups per NMF factor. Default is 0.7.

- cancer_type:

  Character. TCGA cancer type abbreviation (one of `"blca"`, `"luad"`,
  `"skcm"`). Used for two purposes: (1) loading TCGA meta-programs for
  TME state mapping, and (2) when `TF.collection = "ARACNE"`, locating
  the ARACNe network at
  `input/ARACNE/<cancer_type>/network/network.txt`. If `NULL`, the
  meta-program mapping step is skipped (`TME_states` and
  `Metaprograms_reference` are `NULL`) and, for ARACNE, the network is
  auto-detected when only one exists under `input/ARACNE/`.

- return:

  Logical; if TRUE, intermediate matrices and plots are written to the
  "Results/" folder. Default is TRUE.

- verbose:

  Boolen value to whether print or no the function messages

## Value

A list containing:

- Deconvolution:

  A matrix with cell-type proportions (samples as rows, cell types as
  columns); `NULL` if `dt` was supplied.

- TFs_matrix:

  A matrix with TF activity scores (samples as rows, TFs as columns), or
  a list of matrices (one per cohort) when `batch = TRUE`.

- TF_network:

  The TF module network returned by
  [`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md).

- Pathways_scores:

  Pathway activity scores returned by
  [`compute.pathway.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.pathway.activity.md).

- Processed_deconvolution:

  The processed deconvolution subgroups (output of
  [`multideconv::compute.deconvolution.analysis()`](https://verapancaldilab.github.io/multideconv/reference/compute.deconvolution.analysis.html)).

- Cell_groups:

  Output of
  [`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md):
  cell group scores, compositions and projection weights.

- Latent_spaces:

  NMF latent factors returned by
  [`compute.latent_factors()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.latent_factors.md).

- Cells_niches:

  Enriched cell types per latent factor (see
  [`compute_cells_niches()`](https://verapancaldilab.github.io/CellTFusion/reference/compute_cells_niches.md)).

- TME_states:

  Factor-to-meta-program mapping (`factor_mapping` from
  [`map_factors_to_metaprograms()`](https://verapancaldilab.github.io/CellTFusion/reference/map_factors_to_metaprograms.md)),
  or `NULL` if `cancer_type` is `NULL`.

- Metaprograms_reference:

  The TCGA meta-program reference used for mapping, or `NULL` if
  `cancer_type` is `NULL`.

## Details

The pipeline runs with a fixed random seed (123) so results are
reproducible; the caller's random number generator state is restored
when the function returns. With `normalized = TRUE`, sample names are
converted with [`make.names()`](https://rdrr.io/r/base/make.names.html)
(e.g. `"TCGA-XX-1234"` becomes `"TCGA.XX.1234"`), as
[`multideconv::compute.deconvolution()`](https://verapancaldilab.github.io/multideconv/reference/compute.deconvolution.html)
also does; `deconv` rows are matched to the samples in either form.

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
  pval = 0.05
)

# Re-run with previously computed features, skipping straight to cell group construction
res2 <- CellTFusion(
  raw.counts = raw.counts.tuto,
  dt = res$Processed_deconvolution,
  tfs = res$TFs_matrix,
  pathways = res$Pathways_scores,
  normalized = TRUE,
  coldata = traitdata.tuto,
  file_name = "TestRun_rerun",
  minMod = 20,
  corr_mod = 0.25,
  pval = 0.05
)
} # }
```
