# Compute one-step CellTFusion

Runs the full CellTFusion pipeline in a single call: cell-type
deconvolution, TF activity inference, TF module construction, pathway
activity, cell group construction, NMF latent factors, cell niches and
TME state characterization (Hallmark GSEA per factor and mapping to TCGA
meta-programs).

## Usage

``` r
CellTFusion(
  raw_counts,
  deconv = NULL,
  dt = NULL,
  tfs = NULL,
  pathways = NULL,
  normalized = TRUE,
  batch = NULL,
  deconv_methods = c("Quantiseq", "CBSX", "Epidish", "DeconRNASeq", "DWLS"),
  cbsx_mail = NULL,
  cbsx_token = NULL,
  do_parallel = FALSE,
  workers = NULL,
  file_name = NULL,
  TF.collection = "CollecTRI",
  method = "ulm",
  min_targets_size = 5,
  universe = NULL,
  paths = NULL,
  gene_sets = NULL,
  min_module_size = 5,
  corr_mod = 0.9,
  soft_power = NULL,
  consensus_quantile = 0,
  clustering_method = "ward.D2",
  corr = 0.7,
  corr_type = "spearman",
  cells_extra = NULL,
  pval = 0.05,
  padj = FALSE,
  min_cluster_size = 3,
  enrich_thresh = 1.5,
  quantile_cutoff = 0.7,
  cancer_type = NULL,
  return = TRUE,
  verbose = TRUE,
  seed = 123
)
```

## Arguments

- raw_counts:

  A matrix of raw gene expression counts (genes as rows, samples as
  columns). If `normalized = FALSE`, a matrix of TPM-normalized values
  instead (not log-transformed). Always required: even when `dt`, `tfs`
  and `pathways` are all supplied precomputed, the (normalized)
  expression matrix is still needed for the downstream GSEA-based TME
  state characterization step.

- deconv:

  (Optional) A data frame with deconvolution features (cell-type
  proportions as columns x samples as rows). If `NULL` (default), it is
  computed with
  [`multideconv::compute.deconvolution()`](https://verapancaldilab.github.io/multideconv/reference/compute.deconvolution.html)
  using `deconv_methods`. Ignored if `dt` is supplied.

- dt:

  (Optional) A precomputed cell-subgroup object, typically the output of
  [`multideconv::compute.deconvolution.analysis()`](https://verapancaldilab.github.io/multideconv/reference/compute.deconvolution.analysis.html)
  (or the `Processed_deconvolution` element returned by a previous
  `CellTFusion()` run). If supplied, cell-type deconvolution and the
  deconvolution analysis step are both skipped and the pipeline proceeds
  straight to cell group construction using this object. It must have
  the same samples as `raw_counts`.

- tfs:

  (Optional) A precomputed TF activity matrix (samples as rows, TFs as
  columns), typically the `TFs_matrix` element returned by a previous
  `CellTFusion()` run or by
  [`compute.TFs.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.TFs.activity.md).
  When `batch` is given, a named list of such matrices (one per cohort).
  If supplied, TF activity inference is skipped.

- pathways:

  (Optional) A precomputed pathway activity matrix, typically the
  `Pathways_scores` element returned by a previous `CellTFusion()` run
  or by
  [`compute.pathway.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.pathway.activity.md).
  If supplied, pathway activity inference is skipped.

- normalized:

  Logical; if TRUE (default), `raw_counts` are raw counts and are
  normalized to log2(TPM + 1) for the TF activity, pathway activity and
  GSEA steps (for deconvolution they are normalized to TPM only). If
  FALSE, `raw_counts` must already be TPM-normalized and **not**
  log-transformed: the TPM values are used as they are for deconvolution
  and log2(TPM + 1) is computed internally for the other steps. Do not
  pass log-transformed values, as they would be log-transformed a second
  time (a warning is given when the values do not look like TPM, i.e. no
  value above 20 or negative values).

- batch:

  (Optional) A vector with the cohort of each sample (e.g.
  `coldata$Cohort`), in the same order as the columns of `raw_counts`;
  if it is named, it is matched to the samples by name instead. It is
  always treated as categorical. If given, the multi-cohort mode is run:
  TF and pathway activities are computed per cohort, TF modules are
  built from a consensus network across cohorts (see
  `consensus_quantile`), and the cohort is controlled for in the
  deconvolution analysis and cell group construction. Default is `NULL`
  (single cohort); a vector with only one cohort is treated the same
  way.

- deconv_methods:

  A character vector of deconvolution methods to apply. Default is
  `c("Quantiseq", "CBSX", "Epidish", "DeconRNASeq", "DWLS")`. `"CBSX"`
  (CIBERSORTx) is only run when `cbsx_mail` and `cbsx_token` are given;
  otherwise it is skipped with a message.

- cbsx_mail:

  (Optional) Email credential for CIBERSORTx. Needed to run "CBSX";
  without it, CIBERSORTx is skipped.

- cbsx_token:

  (Optional) Token credential for CIBERSORTx. Needed to run "CBSX";
  without it, CIBERSORTx is skipped.

- do_parallel:

  Logical; whether to run the deconvolution methods that support it
  (CIBERSORTx and DWLS) in parallel, one signature per worker (passed to
  [`multideconv::compute.deconvolution()`](https://verapancaldilab.github.io/multideconv/reference/compute.deconvolution.html)
  as `doParallel`). Default is FALSE.

- workers:

  (Optional) Number of parallel workers when `do_parallel = TRUE`. If
  `NULL` (default), the number of available cores minus one is used.

- file_name:

  (Optional) Name used in the output files saved in the "Results/"
  directory. If `NULL` and `return = TRUE`, `"CellTFusion"` is used.

- TF.collection:

  Character. The source of the TF-target network. Options are
  `"CollecTRI"` (default), `"Dorothea"`, or `"ARACNE"`.

  - `"CollecTRI"` uses the prebuilt collection from OmnipathR and
    `"Dorothea"` the one from the `dorothea` package.

  - `"ARACNE"` reads a network file (tab-separated with `Regulator` and
    `Target` columns, and a `mor` column with the sign of each
    interaction; see
    [`compute.TFs.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.TFs.activity.md))
    from `input/ARACNE/<cancer_type>/network/network.txt`.

- method:

  Character. The `decoupleR` statistic used to score the TFs and the
  PROGENy pathways: `"ulm"` (univariate linear model, default), `"mlm"`
  (multivariate linear model), `"wsum"` (weighted sum) or `"consensus"`
  (consensus of the three). Passed to
  [`compute.TFs.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.TFs.activity.md)
  and
  [`compute.pathway.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.pathway.activity.md).

- min_targets_size:

  Integer. Minimum number of target genes per regulon, passed to
  [`compute.TFs.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.TFs.activity.md).
  Default is 5.

- universe:

  Optional. A user-specified data frame of TF-target interactions
  (columns `source`, `target` and `mor`, the mode of regulation: 1
  activation, -1 repression). If not provided, the function will fetch
  the relevant network based on the `TF.collection` argument.

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

- min_module_size:

  Integer; minimum module size for WGCNA module detection. Default is 5.

- corr_mod:

  Numeric; correlation threshold for merging TF modules. Default is 0.9.

- soft_power:

  (Optional) Soft-thresholding power used to build the TF network,
  passed to
  [`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md):
  one value, or one value per cohort (in the order the cohorts first
  appear in `batch`) in multi-cohort mode. If `NULL` (default), it is
  chosen automatically.

- consensus_quantile:

  Numeric between 0 and 1, used in multi-cohort mode only: quantile of
  the topological overlaps of each pair of TFs across cohorts used as
  their consensus (0, the default, is the minimum: two TFs are close
  only if they are close in every cohort; higher values are less
  strict). Passed to
  [`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md).

- clustering_method:

  Linkage method of the hierarchical clusterings
  ([`stats::hclust()`](https://rdrr.io/r/stats/hclust.html)): of the
  TFs, to build the TF modules
  ([`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md)),
  and of the deconvolution features of each TF module, to build the cell
  groups
  ([`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md)).
  Default is `"ward.D2"`.

- corr:

  Numeric; correlation threshold used in the deconvolution analysis.
  Default is 0.7.

- corr_type:

  Correlation type used in every step that computes correlations:
  `"spearman"` (default) or `"pearson"`. Passed to
  [`multideconv::compute.deconvolution.analysis()`](https://verapancaldilab.github.io/multideconv/reference/compute.deconvolution.analysis.html)
  (deconvolution subgroups),
  [`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md)
  (TF network and module merging) and
  [`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md)
  (correlation between TF modules and deconvolution features).

- cells_extra:

  (Optional) A character vector with the cell names to consider that are
  not included in the nomenclature of multideconv (see R package). This
  includes the group names created with
  [`multideconv::aggregate_cell_groups()`](https://verapancaldilab.github.io/multideconv/reference/aggregate_cell_groups.html)
  on `deconv` (e.g. `"Lymphocytes"`): if they are not listed here they
  are discarded.

- pval:

  Numeric; p-value threshold used when building cell groups (selection
  of the deconvolution features significantly correlated with each TF
  module). Default is 0.05.

- padj:

  Logical; if TRUE, the p-values used when building cell groups are
  corrected for multiple testing (Benjamini-Hochberg, across the
  deconvolution features for each TF module) and compared with `pval`.
  Passed to
  [`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md).
  Default is FALSE.

- min_cluster_size:

  Integer; minimum number of features of a cell group (minimum cluster
  size of the dynamic tree cut), passed to
  [`construct_cell_groups()`](https://verapancaldilab.github.io/CellTFusion/reference/construct_cell_groups.md).
  Default is 3.

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
  `input/ARACNE/<cancer_type>/network/network.txt`. If `NULL`, or if no
  TCGA meta-programs exist for that cancer type (a warning is given),
  the meta-program mapping step is skipped (`TME_states` and
  `Metaprograms_reference` are `NULL`). For ARACNE, if `NULL`, the
  network is auto-detected when only one exists under `input/ARACNE/`.

- return:

  Logical; if TRUE, intermediate matrices and plots are written to the
  "Results/" folder. Default is TRUE.

- verbose:

  Logical; whether to print the progress and the console output and
  messages of every step (warnings are always shown). Default is TRUE.

- seed:

  Integer. Random seed of the pipeline (deconvolution methods with
  random steps, NMF and the GSEA permutations), so that results are
  reproducible. Default is 123.

## Value

A list containing:

- Deconvolution:

  A matrix with cell-type proportions (samples as rows, cell types as
  columns); `NULL` if `dt` was supplied without `deconv`.

- TFs_matrix:

  A matrix with TF activity scores (samples as rows, TFs as columns), or
  a list of matrices (one per cohort) when `batch` is given. Single
  cohort: raw scores (not scaled), which do not depend on the other
  samples (e.g. for machine learning features, to be scaled with the
  training samples). With `batch`: scores z-scored within each cohort
  (correction of cohort shifts).

- TF_network:

  The TF module network returned by
  [`compute.WTCNA()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.WTCNA.md).
  Its module matrix is always in the sample order of `raw_counts` (also
  in multi-cohort mode).

- Pathways_scores:

  Pathway activity scores returned by
  [`compute.pathway.activity()`](https://verapancaldilab.github.io/CellTFusion/reference/compute.pathway.activity.md):
  raw in single-cohort mode, z-scored within each cohort with `batch`
  (as `TFs_matrix`).

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
  or `NULL` if the mapping was skipped (see `cancer_type`).

- Metaprograms_reference:

  The TCGA meta-program reference used for mapping, or `NULL` if the
  mapping was skipped.

## Details

The pipeline runs with a fixed random seed (`seed`) so results are
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

res <- CellTFusion(
  raw_counts = raw.counts.tuto,
  normalized = TRUE,
  deconv_methods = c("Quantiseq", "DeconRNASeq"),
  file_name = "TestRun",
  min_targets_size = 15,
  min_module_size = 20,
  corr_mod = 0.25,
  corr = 0.7,
  pval = 0.05
)

# Re-run with previously computed features, skipping straight to cell group construction
res2 <- CellTFusion(
  raw_counts = raw.counts.tuto,
  dt = res$Processed_deconvolution,
  tfs = res$TFs_matrix,
  pathways = res$Pathways_scores,
  normalized = TRUE,
  file_name = "TestRun_rerun",
  min_module_size = 20,
  corr_mod = 0.25,
  pval = 0.05
)
} # }
```
