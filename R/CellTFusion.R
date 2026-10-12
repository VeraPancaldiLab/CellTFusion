
utils::globalVariables(c("Trait", "Value" ,"level", ".", "Cells_level", "PC1", "Features", "tf", "Module", "shap_value", "mean_shap", "direction", "Feature", "Impact", "values", "ind", "Freq", "new_column", ".data", "id", "value", "P", "p", "sig_p", "r", "i",
                         ":=", "MFP", "MP", "NES", "Regulator", "TME_subtype", "Target", "as.hclust",
                         "condition", "dist", "feature_label",
                         "group", "gs_name", "hclust", "kruskal.test",
                         "lm", "meta_programs", "module", "mor", "mtext", "order.dendrogram", "padj",
                         "pathway", "proportion", "quantile", "read.delim", "residuals", "score",
                         "score_level", "setNames", "svg", "target", "title",
                         "n", "n_distinct"))

#' @importFrom grDevices colorRampPalette svg
#' @importFrom graphics mtext title
#' @importFrom stats as.hclust dist hclust kruskal.test lm order.dendrogram quantile residuals setNames
#' @importFrom utils read.delim
#' @importFrom RColorBrewer brewer.pal
#' @importFrom colorspace rainbow_hcl
#' @importFrom fgsea fgsea
#' @importFrom ppcor pcor.test
#' @importFrom multideconv compute.deconvolution compute.deconvolution.analysis get_cell_type_nomenclature replicate_deconvolution_subgroups
#' @rawNamespace export(compute.TFs.activity)
#' @rawNamespace export(compute.WTCNA)
#' @rawNamespace export(compute.pathway.activity)
#' @rawNamespace export(compute.metadata.association)
#' @rawNamespace export(compute.latent_factors)
#' @rawNamespace export(compute.modules.relationship)
#' @rawNamespace export(compute.modules.enrichment)
#' @rawNamespace export(compute.composition.matrix)
#' @rawNamespace export(compute.test.set)
#' @rawNamespace export(compute.survival.analysis)
#' @rawNamespace export(identify.cell.groups)
NULL


#' Compute one-step CellTFusion
#'
#' Runs the full CellTFusion pipeline in a single call: cell-type deconvolution, TF activity inference,
#' TF module construction, pathway activity, cell group construction, NMF latent factors, cell niches and
#' TME state characterization (Hallmark GSEA per factor and mapping to TCGA meta-programs).
#'
#' @param raw.counts A matrix of raw gene expression counts (genes as rows, samples as columns). If
#'   \code{normalized = FALSE}, a matrix of TPM-normalized values instead (not log-transformed).
#'   Always required: even when \code{dt}, \code{tfs} and \code{pathways} are all supplied precomputed,
#'   the (normalized) expression matrix is still needed for the downstream GSEA-based TME state
#'   characterization step.
#' @param deconv (Optional) A data frame with deconvolution features (cell-type proportions as columns x samples as rows).
#'   If \code{NULL} (default), it is computed with \code{multideconv::compute.deconvolution()} using
#'   \code{deconv_methods}. Ignored if \code{dt} is supplied.
#' @param dt (Optional) A precomputed cell-subgroup object, typically the output of
#'   \code{multideconv::compute.deconvolution.analysis()} (or the \code{Processed_deconvolution} element
#'   returned by a previous \code{CellTFusion()} run). If supplied, cell-type deconvolution and the
#'   deconvolution analysis step are both skipped and the pipeline proceeds straight to cell group
#'   construction using this object.
#' @param tfs (Optional) A precomputed TF activity matrix (samples as rows, TFs as columns), typically
#'   the \code{TFs_matrix} element returned by a previous \code{CellTFusion()} run or by
#'   \code{compute.TFs.activity()}. When \code{batch} is given, a named list of such matrices (one per cohort).
#'   If supplied, TF activity inference is skipped.
#' @param pathways (Optional) A precomputed pathway activity matrix, typically the \code{Pathways_scores}
#'   element returned by a previous \code{CellTFusion()} run or by \code{compute.pathway.activity()}.
#'   If supplied, pathway activity inference is skipped.
#' @param normalized Logical; if TRUE (default), \code{raw.counts} are raw counts and are normalized to
#'   log2(TPM + 1) for the TF activity, pathway activity and GSEA steps (for deconvolution they are normalized to
#'   TPM only). If FALSE, \code{raw.counts} must already be TPM-normalized and \strong{not} log-transformed: the
#'   TPM values are used as they are for deconvolution and log2(TPM + 1) is computed internally for the other
#'   steps. Do not pass log-transformed values, as they would be log-transformed a second time (a warning is
#'   given when the values do not look like TPM, i.e. no value above 20 or negative values).
#' @param batch (Optional) A vector with the cohort of each sample (e.g. \code{coldata$Cohort}), in the same order
#'   as the columns of \code{raw.counts}; if it is named, it is matched to the samples by name instead. It is
#'   always treated as categorical. If given, the multi-cohort mode is run: TF activity is computed per cohort,
#'   TF modules are built with consensus WGCNA across cohorts, and the cohort is controlled for in the
#'   deconvolution analysis and cell group construction. Default is \code{NULL} (single cohort); a vector with
#'   only one cohort is treated the same way.
#' @param deconv_methods A character vector of deconvolution methods to apply. Default is
#'   \code{c("Quantiseq", "CBSX", "Epidish", "DeconRNASeq", "DWLS")}. \code{"CBSX"} (CIBERSORTx) is only run
#'   when \code{cbsx.mail} and \code{cbsx.token} are given; otherwise it is skipped with a message.
#' @param cbsx.mail (Optional) Email credential for CIBERSORTx. Needed to run "CBSX"; without it, CIBERSORTx is skipped.
#' @param cbsx.token (Optional) Token credential for CIBERSORTx. Needed to run "CBSX"; without it, CIBERSORTx is skipped.
#' @param doParallel Logical; whether to run the deconvolution methods that support it (CIBERSORTx and DWLS) in
#'   parallel, one signature per worker (passed to \code{multideconv::compute.deconvolution()}). Default is FALSE.
#' @param workers (Optional) Number of parallel workers when \code{doParallel = TRUE}. If \code{NULL} (default),
#'   the number of available cores minus one is used.
#' @param file_name (Optional) Name used in the output files saved in the "Results/" directory. If \code{NULL}
#'   and \code{return = TRUE}, \code{"CellTFusion"} is used.
#' @param TF.collection Character. The source of the TF-target network. Options are `"CollecTRI"` (default), `"Dorothea"`, or `"ARACNE"`.
#' - `"CollecTRI"` uses the prebuilt collection from OmnipathR and `"Dorothea"` the one from the `dorothea` package.
#' - `"ARACNE"` reads a network file (tab-separated with `Regulator` and `Target` columns, and a `mor` column with
#'   the sign of each interaction; see \code{compute.TFs.activity()}) from
#'   \code{input/ARACNE/<cancer_type>/network/network.txt}.
#' @param min_targets_size Integer. Minimum number of target genes per regulon, passed to
#'   \code{compute.TFs.activity()}. Default is 5.
#' @param universe Optional. A user-specified data frame of TF-target interactions (columns \code{source}, \code{target} and
#'   \code{mor}, the mode of regulation: 1 activation, -1 repression). If not provided, the function will fetch the relevant network based on the `TF.collection` argument.
#' @param paths Optional. A user-specified data frame of pathways gene sets. If not provided, the function will fetch the relevant pathways based on `PROGENy`.
#' @param gene_sets Optional. A named list of custom gene sets (character vectors of gene symbols) passed to
#'   \code{compute.pathway.activity()}'s \code{gene_sets} argument for GSVA-based scoring. If \code{NULL}, only PROGENy is used.
#' @param minMod Integer; minimum module size for WGCNA module detection. Default is 5.
#' @param corr_mod Numeric; correlation threshold for merging TF modules. Default is 0.9.
#' @param softPower (Optional) Soft-thresholding power used to build the TF network, passed to
#'   \code{compute.WTCNA()}: one value, or one value per cohort (in the order the cohorts first appear in
#'   \code{batch}) in multi-cohort mode. If \code{NULL} (default), it is chosen automatically.
#' @param corr Numeric; correlation threshold used in the deconvolution analysis. Default is 0.7.
#' @param corr_type Correlation type used in deconvolution analysis. Default is \code{"spearman"}.
#' @param cells_extra (Optional) A character vector with the cell names to consider that are not included in the nomenclature of multideconv (see R package).
#'   This includes the group names created with \code{multideconv::aggregate_cell_groups()} on \code{deconv} (e.g. \code{"Lymphocytes"}):
#'   if they are not listed here they are discarded.
#' @param pval Numeric; p-value threshold used when building cell groups (selection of the deconvolution
#'   features significantly correlated with each TF module). Default is 0.05.
#' @param min_cluster_size Integer; minimum number of features of a cell group (minimum cluster size of the
#'   dynamic tree cut), passed to \code{construct_cell_groups()}. Default is 3.
#' @param deep_split Logical; if TRUE (default), the dynamic tree cut that defines the cell groups splits the
#'   clusters more finely. Passed to \code{construct_cell_groups()}.
#' @param enrich_thresh Numeric. Minimum enrichment ratio (foreground/background cell-type frequency)
#'   required to include a cell type in a latent factor's niche. Default is 1.5.
#' @param quantile_cutoff Numeric between 0 and 1. Quantile threshold for selecting top-contributing
#'   cell groups per NMF factor. Default is 0.7.
#' @param return Logical; if TRUE, intermediate matrices and plots are written to the "Results/" folder. Default is TRUE.
#' @param cancer_type Character. TCGA cancer type abbreviation (one of \code{"blca"}, \code{"luad"}, \code{"skcm"}).
#'   Used for two purposes: (1) loading TCGA meta-programs for TME state mapping, and (2) when
#'   \code{TF.collection = "ARACNE"}, locating the ARACNe network at
#'   \code{input/ARACNE/<cancer_type>/network/network.txt}. If \code{NULL}, or if no TCGA meta-programs
#'   exist for that cancer type (a warning is given), the meta-program mapping step is skipped
#'   (\code{TME_states} and \code{Metaprograms_reference} are \code{NULL}). For ARACNE, if \code{NULL}, the
#'   network is auto-detected when only one exists under \code{input/ARACNE/}.
#' @param verbose Logical; whether to print the progress and the console output and messages of every step
#'   (warnings are always shown). Default is TRUE.
#'
#' @details
#' The pipeline runs with a fixed random seed (123) so results are reproducible; the caller's random
#' number generator state is restored when the function returns. With \code{normalized = TRUE}, sample
#' names are converted with \code{make.names()} (e.g. \code{"TCGA-XX-1234"} becomes \code{"TCGA.XX.1234"}),
#' as \code{multideconv::compute.deconvolution()} also does; \code{deconv} rows are matched to the samples
#' in either form.
#'
#' @return A list containing:
#' \describe{
#'   \item{Deconvolution}{A matrix with cell-type proportions (samples as rows, cell types as columns); \code{NULL} if \code{dt} was supplied without \code{deconv}.}
#'   \item{TFs_matrix}{A matrix with TF activity scores (samples as rows, TFs as columns), or a list of matrices (one per cohort) when \code{batch} is given.}
#'   \item{TF_network}{The TF module network returned by \code{compute.WTCNA()}. Its module matrix is always
#'     in the sample order of \code{raw.counts} (also in multi-cohort mode).}
#'   \item{Pathways_scores}{Pathway activity scores returned by \code{compute.pathway.activity()}.}
#'   \item{Processed_deconvolution}{The processed deconvolution subgroups (output of \code{multideconv::compute.deconvolution.analysis()}).}
#'   \item{Cell_groups}{Output of \code{construct_cell_groups()}: cell group scores, compositions, projection weights
#'     and the relation of each cell group to its TF module.}
#'   \item{Latent_spaces}{NMF latent factors returned by \code{compute.latent_factors()}.}
#'   \item{Cells_niches}{Enriched cell types per latent factor (see \code{compute_cells_niches()}).}
#'   \item{TME_states}{Factor-to-meta-program mapping (\code{factor_mapping} from \code{map_factors_to_metaprograms()}), or \code{NULL} if the mapping was skipped (see \code{cancer_type}).}
#'   \item{Metaprograms_reference}{The TCGA meta-program reference used for mapping, or \code{NULL} if the mapping was skipped.}
#' }
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#' data("raw.counts.tuto")
#'
#' res <- CellTFusion(
#'   raw.counts = raw.counts.tuto,
#'   normalized = TRUE,
#'   deconv_methods = c("Quantiseq", "DeconRNASeq"),
#'   file_name = "TestRun",
#'   min_targets_size = 15,
#'   minMod = 20,
#'   corr_mod = 0.25,
#'   corr = 0.7,
#'   pval = 0.05
#' )
#'
#' # Re-run with previously computed features, skipping straight to cell group construction
#' res2 <- CellTFusion(
#'   raw.counts = raw.counts.tuto,
#'   dt = res$Processed_deconvolution,
#'   tfs = res$TFs_matrix,
#'   pathways = res$Pathways_scores,
#'   normalized = TRUE,
#'   file_name = "TestRun_rerun",
#'   minMod = 20,
#'   corr_mod = 0.25,
#'   pval = 0.05
#' )
#'}
#'
CellTFusion = function(raw.counts, deconv = NULL, dt = NULL, tfs = NULL, pathways = NULL, normalized = T, batch = NULL, deconv_methods = c("Quantiseq", "CBSX", "Epidish", "DeconRNASeq", "DWLS"), cbsx.mail = NULL, cbsx.token = NULL, doParallel = F, workers = NULL, file_name = NULL,
                       TF.collection = "CollecTRI", min_targets_size = 5, universe = NULL, paths = NULL, gene_sets = NULL, minMod = 5, corr_mod = 0.9, softPower = NULL, corr = 0.7, corr_type = "spearman", cells_extra = NULL, pval = 0.05, min_cluster_size = 3, deep_split = T, enrich_thresh = 1.5, quantile_cutoff = 0.7, cancer_type = NULL, return = T, verbose = T){

  withr::local_preserve_seed()
  set.seed(123)

  # Run a step, silencing its console output and messages when verbose = FALSE
  quiet = function(expr) if(verbose) expr else { utils::capture.output(res <- suppressMessages(expr)); res }

  #Normalize counts
  if(normalized == T){
    counts.norm = quiet(data.frame(ADImpute::NormalizeTPM(raw.counts, log = T)))
  }else{
    if(max(raw.counts, na.rm = T) <= 20 || min(raw.counts, na.rm = T) < 0){ # log2(TPM + 1) cannot exceed 20
      warning("raw.counts does not look like TPM values (range ", signif(min(raw.counts, na.rm = T), 3), " to ", signif(max(raw.counts, na.rm = T), 3),
              "): with normalized = FALSE, TPM-normalized values that are not log-transformed are expected, and they are log-transformed here")
    }
    counts.norm = log2(raw.counts + 1) # raw.counts are TPM: log2(TPM + 1) as above (deconvolution uses the TPM as they are)
  }

  if(return && is.null(file_name)){
    file_name = "CellTFusion"
  }

  #Cohort of each sample (always categorical), in the sample order of the counts
  if(!is.null(batch)){
    if(is.data.frame(batch) && ncol(batch) == 1) batch = batch[[1]] # One-column data frame or tibble
    if(is.data.frame(batch) || length(batch) != ncol(counts.norm)){
      stop("batch must be a vector with the cohort of each sample (", ncol(counts.norm), " samples in raw.counts)")
    }
    if(!is.null(names(batch))){ # Named: matched to the samples by name
      idx = match(make.names(colnames(counts.norm)), make.names(names(batch)))
      if(anyNA(idx)) stop("The names of batch do not match the sample names of raw.counts")
      batch = batch[idx]
    }
    if(anyNA(batch)) stop("batch contains missing values")
    batch_vec = factor(as.character(batch), levels = unique(as.character(batch)))
    if(nlevels(batch_vec) < 2){ # A single cohort: same as no batch
      if(verbose){
        cat("batch contains a single cohort: running the single-cohort mode\n")
      }
      batch_vec = NULL
      if(is.list(tfs) && !is.data.frame(tfs) && length(tfs) == 1) tfs = tfs[[1]] # Precomputed per-cohort list with its only cohort
    }
  }else{
    batch_vec = NULL
  }

  #OmnipathR writes its own log lines to the console: only warnings and errors when verbose = FALSE
  if(!verbose && requireNamespace("OmnipathR", quietly = TRUE)){
    omnipath_level = getOption("omnipathr.console_loglevel", "success")
    OmnipathR::omnipath_set_console_loglevel("warn")
    on.exit(OmnipathR::omnipath_set_console_loglevel(omnipath_level), add = TRUE)
  }

  #TCGA meta-programs only exist for some cancer types: otherwise the mapping step is skipped
  map_metaprograms = !is.null(cancer_type) && system.file("extdata", paste0("TCGA_meta_programs_", tolower(cancer_type), ".RData"), package = "CellTFusion") != ""
  if(!is.null(cancer_type) && !map_metaprograms){
    warning("No TCGA meta-programs available for cancer type '", cancer_type, "' (available: blca, luad, skcm): the mapping to meta-programs is skipped")
  }

  #Deconvolution (skipped entirely if precomputed cell subgroups `dt` are supplied)
  if(is.null(dt)){
    if(is.null(deconv)){
      if(verbose){
        cat("Calculating cell type deconvolution............................................................\n")
      }
      #CIBERSORTx needs credentials: without them it is skipped
      if("CBSX" %in% deconv_methods && (is.null(cbsx.mail) || is.null(cbsx.token))){
        message("CIBERSORTx (CBSX) is skipped: it needs the credentials cbsx.mail and cbsx.token")
        deconv_methods = setdiff(deconv_methods, "CBSX")
        if(length(deconv_methods) == 0) stop("No deconvolution method to run: CBSX needs cbsx.mail and cbsx.token")
      }
      deconv = quiet(compute.deconvolution(raw.counts, normalized = normalized, methods = deconv_methods, credentials.mail = cbsx.mail, credentials.token = cbsx.token, doParallel = doParallel, workers = workers, file_name = file_name, return = return))
    }

    # Validate sample alignment between counts and deconv
    if (!isTRUE(all.equal(make.names(colnames(counts.norm), unique = TRUE), make.names(rownames(deconv), unique = TRUE)))) {
      stop("Sample mismatch between counts columns and deconvolution rows.")
    }
    rownames(deconv) = colnames(counts.norm) # multideconv converts sample names with make.names() (e.g. "TCGA-XX" -> "TCGA.XX")
  }

  #TF activity (skipped entirely if a precomputed `tfs` matrix is supplied)
  if(is.null(tfs)){
    if (verbose) {
      cat("\nCalculating TF activity............................................................\n")
    }

    if (is.null(batch_vec)) {  ## No batch -> single matrix
      tfs <- quiet(compute.TFs.activity(counts.norm, TF.collection = TF.collection, min_targets_size = min_targets_size, universe = universe, cancer.type = cancer_type,return = return, file.name = file_name))
    }else {
      cohorts <- split(seq_len(ncol(counts.norm)), batch_vec)

      tfs <- lapply(names(cohorts), function(cohort) {
        idx <- cohorts[[cohort]]
        quiet(compute.TFs.activity(counts.norm[, idx, drop = FALSE], TF.collection = TF.collection, min_targets_size = min_targets_size, universe = universe, cancer.type = cancer_type,file.name = paste0(file_name, "_", cohort), return = return))
      })

      names(tfs) <- names(cohorts)
    }
  }

  # 1. TFs network construction
  if(verbose){
    cat("\nConstructing TF network............................................................\n")
  }
  network = quiet(compute.WTCNA(TFs.matrix = tfs, batch = !is.null(batch_vec), network.type = "signed", clustering.method = "average", minMod = minMod, corr_mod = corr_mod, cor_type = "p", softPower = softPower, return = return, file.name = file_name))
  if(!is.null(batch_vec)){ # Consensus module scores are stacked cohort by cohort: back to the sample order of the counts
    idx = match(make.names(colnames(counts.norm)), make.names(rownames(network[["TFs module matrix"]])))
    if(!anyNA(idx)) network[["TFs module matrix"]] = network[["TFs module matrix"]][idx, , drop = F]
  }

  # 2. Pathways activity inference (skipped if a precomputed `pathways` matrix is supplied); returned for downstream analyses
  if(is.null(pathways)){
    if(verbose){
      cat("\nCalculating pathway activities............................................................\n")
    }
    pathways = quiet(compute.pathway.activity(counts.norm, gene_sets = gene_sets, paths = paths, return = return, file.name = file_name))
  }

  # 3. Deconvolution analysis / cell subgroups (skipped if precomputed `dt` is supplied)
  if(is.null(dt)){
    if(verbose){
      cat("\nPerforming deconvolution analysis............................................................\n")
    }

    dt = compute.deconvolution.analysis(deconv, corr = corr, corr_type = corr_type, batch = batch_vec, cells_extra = cells_extra, file_name = file_name, return = return, verbose = FALSE)
  }

  # 4. Cell groups construction and scores
  if(verbose){
    cat("\nCell groups identification............................................................\n")
  }
  cell.groups = quiet(construct_cell_groups(network, dt, batch = batch_vec, pval = pval, clustering.method = "ward.D2", dendrogram_file = file_name, return_dendrogram = return, return = return, cells_extra = cells_extra, min_cluster_size = min_cluster_size, deep_split = deep_split))

  # 5. Cell groups latent spaces
  if(verbose){
    cat("\nLatent spaces calculation............................................................\n")
  }
  latent_spaces <- quiet(compute.latent_factors(cell.groups[[1]], rank = NULL, seed = 123, file_name = file_name, return = return))

  # 6. Extract cell niches
  if(verbose){
    cat("\nExtract cell niches from latent spaces............................................................\n")
  }
  cells_niches = compute_cells_niches(latent_spaces, dt, cell.groups, enrich_thresh = enrich_thresh, quantile_cutoff = quantile_cutoff, cells_extra = cells_extra, return = return, file_name = file_name)

  # 7. Find TME states
  if(verbose){
    cat("\nCompute TME states............................................................\n")
  }
  gsea_results <- quiet(compute_factor_gsea(
    RNA.tpm      = counts.norm,          # genes x samples
    features_df  = latent_spaces$Z,      # samples x factors
    plot_dot     = return,
    top_n        = 10,
    file_name    = file_name,
    batch        = batch_vec             # cohort of each sample (NULL with a single cohort)
  ))

  if(!map_metaprograms){
    if(verbose){
      cat("\nNo cancer_type with TCGA meta-programs given: skipping mapping to TCGA meta-programs.\n")
    }
    metaprograms_mapping <- list(factor_mapping = NULL, reference = NULL)
  }else{
    if(verbose){
      cat("\nMapping to metaprograms............................................................\n")
    }
    metaprograms_mapping <- quiet(map_factors_to_metaprograms(
      gsea_study   = gsea_results,
      cancer_type  = cancer_type,
      plot = return, file_name = file_name
    ))
  }

  if(verbose){
    cat("\nEverything done!", if(return) "Results are saved in Results/ folder", "............................................................\n")
  }

  res = list("Deconvolution" = deconv, "TFs_matrix" = tfs, "TF_network" = network, "Pathways_scores" = pathways,
             "Processed_deconvolution" = dt, "Cell_groups" = cell.groups, "Latent_spaces" = latent_spaces, 
             "Cells_niches" = cells_niches, "TME_states" = metaprograms_mapping$factor_mapping, "Metaprograms_reference" = metaprograms_mapping$reference)

  return(res)

}

#' Compute cell group scores from deconvolution and TF module network
#'
#' This function identifies cell groups based on dendrogram cuts, computes composite scores for each group
#' using deconvolution features and TF module networks, and optionally exports the results.
#'
#' @param deconvolution A data frame with deconvolution features (samples as rows, cell-type features as columns).
#' This is usually the first element returned by \code{multideconv::compute.deconvolution.analysis()}. It must
#' contain all the features of \code{cell.dendrograms}.
#'
#' @param cell.dendrograms A named list of dendrogram (\code{hclust}) objects, each corresponding to a TF module, typically
#' returned by \code{identify.cell.groups()}.
#'
#' @param tfs.module.network A list containing network information of transcription factor (TF) modules, as
#' obtained from \code{compute.WTCNA()} (the elements \code{TFs module matrix}, \code{TFs per module} and
#' \code{TFs_matrix} are used).
#' @param batch Optional vector indicating batch assignment for samples.
#' @param return Logical; if TRUE (default), writes CSV files with cell group compositions and scores
#' to the "Results/" folder (named with \code{dendrogram_file}; the folder is created if needed).
#' @param dendrogram_file Optional character. Suffix of the output file names: the dendrogram PDF (see
#'   \code{return_dendrogram}) and the CSV files (see \code{return}).
#' @param return_dendrogram Logical. If TRUE, saves a PDF of the colored cell-group dendrograms to
#'   \code{Results/Dendrogram_color_clusters_<dendrogram_file>.pdf} (only when \code{dendrogram_file} is set). Default FALSE.
#' @param cells_extra Optional character vector with the cell names that are not in the nomenclature of multideconv,
#'   used to list the cell types of each cell group.
#' @param min_cluster_size Integer. Minimum number of features of a cell group (minimum cluster size of the
#'   dynamic tree cut of each dendrogram). Default is 3.
#' @param deep_split Logical. If TRUE (default), the dynamic tree cut splits the clusters more finely.
#'
#' @return An unnamed list of four elements:
#' \itemize{
#'   \item \code{[[1]]}: A data frame with the composite scores of all identified cell groups across samples. A
#'     positive score means that the sample is on the active side of the TF module of that cell group.
#'   \item \code{[[2]]}: A named list of vectors indicating the composition (original features) of each cell group.
#'   \item \code{[[3]]}: A named list of CCA projection parameters (\code{xcoef}, \code{train_means}, \code{train_sds}) for each cell group.
#'   \item \code{[[4]]}: A data frame with the relation of each cell group to its TF module (see
#'     \code{construct_cell_groups()}).
#' }
#' Cell groups are named \code{Dendrogram_<module>_group_<n>}, where \code{n} is the cluster number of the
#' dendrogram of that module (the number shown in the dendrogram PDF); the first three elements use the same names.
#' An error is raised if no cluster is found in any module. Features that the tree cut does not assign to any cluster
#' are not part of any cell group: a message reports how many there are per module.
#' If \code{return=TRUE}, two CSV files will be created:
#' \itemize{
#'   \item \code{Results/Cell.groups.composition_<dendrogram_file>.csv}: A table showing the composition of each cell
#'     group, its TF module, its relation to it (positive or negative, and the correlation) and its cell types.
#'   \item \code{Results/Cell.groups.scores_<dendrogram_file>.csv}: A matrix of cell group scores across samples.
#' }
#' Without \code{dendrogram_file}, they are named \code{Cell.groups.composition.csv} and \code{Cell.groups.scores.csv}.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' deconv_results <- multideconv::compute.deconvolution.analysis(...)
#' tf_network <- compute.WTCNA(...)
#' dendrograms <- identify.cell.groups(...)
#'
#' cell.groups <- cell.groups.computation(
#'   deconvolution = deconv_results[[1]],
#'   cell.dendrograms = dendrograms,
#'   tfs.module.network = tf_network,
#'   return = TRUE
#' )
#' }
cell.groups.computation = function(deconvolution, cell.dendrograms, tfs.module.network, batch = NULL, return = T, dendrogram_file = NULL, return_dendrogram = FALSE, cells_extra = NULL, min_cluster_size = 3, deep_split = TRUE){

  if(length(cell.dendrograms) == 0){
    stop("No dendrograms given: no TF module has deconvolution features to cluster")
  }
  missing_features = setdiff(unlist(lapply(cell.dendrograms, function(x) as.hclust(x)$labels)), colnames(deconvolution))
  if(length(missing_features) > 0){
    stop("Features of the dendrograms that are not columns of deconvolution: ", paste(utils::head(missing_features, 5), collapse = ", "),
         if(length(missing_features) > 5) paste0(", ... (", length(missing_features), " in total)"))
  }
  if(return || return_dendrogram){
    dir.create("Results", showWarnings = FALSE)
  }

  cell.groups = list()
  unassigned = c()
  cuts_per_dendrogram = calculate_dendrogram_cuts(cell.dendrograms, deep_split = deep_split, min_cluster_size = min_cluster_size)

  if(return_dendrogram){
    plot_dendrogram_clusters(cell.dendrograms, cuts_per_dendrogram, file_name = dendrogram_file)
  }

  for(k in 1:length(cell.dendrograms)){
    clusters = cuts_per_dendrogram[[k]]
    names(clusters) = as.hclust(cell.dendrograms[[k]])$labels # Same order as the cluster labels (original feature order, not the leaf order of the plot)

    cluster_ids = sort(unique(clusters[clusters > 0]))  # exclude unassigned (label 0)
    unassigned[names(cell.dendrograms)[k]] = sum(clusters == 0)

    y = list()
    x = list()
    z = list()
    r = c()

    for (j in seq_along(cluster_ids)) {
      cells = names(clusters)[clusters == cluster_ids[j]]
      y[[j]] = cells
      ###################################################Compute score for each cell group
      pca_group = deconvolution[, colnames(deconvolution) %in% cells, drop = F]
      color = names(cell.dendrograms)[k]
      score = compute_composite_score(pca_group, color, tfs.module.network, batch = batch)
      x[[j]] = score[[1]]
      z[[j]] = score[[2]]
      r[j] = score[[3]]
      y[[j]] = rownames(score[[2]]$xcoef) # Features with a weight (the CCA can drop collinear features)
    }

    if(length(x) != 0){
      ###################################################Naming of cell groups: cluster number of the dendrogram (same number as in the dendrogram plot)
      names(x) = names(y) = names(z) = paste0("Dendrogram_", names(cell.dendrograms)[k], "_group_", cluster_ids)
      cell.groups[[k]] = list(x, y, z, r)
    }else{
      cell.groups[[k]] <- NA #To avoid the problem when putting NULL (list ignore it and not consider it an element)
    }
  }

  names(cell.groups) = names(cell.dendrograms) #Named before removing NA to avoid omitting dendrogram names
  cell.groups <- cell.groups[!sapply(cell.groups, function(x) all(is.na(x)))]

  if(length(cell.groups)==0){
    stop("No cell groups found with significant features")
  }

  if(sum(unassigned) > 0){
    message(sum(unassigned), " of ", sum(lengths(cuts_per_dendrogram)), " features were not assigned to any cluster by the tree cut and are not part of any cell group (",
            paste0(names(unassigned)[unassigned > 0], ": ", unassigned[unassigned > 0], collapse = ", "), ")")
  }

  cell.groups.values <- lapply(cell.groups, function(x) x[[1]])
  cell.groups.names <- lapply(cell.groups, function(x) x[[2]])
  cell.groups.loadings <- lapply(cell.groups, function(x) x[[3]])

  ###################################################Relation of each cell group to its TF module
  cell.groups.relation <- unlist(lapply(unname(cell.groups), function(x) x[[4]]))
  module.relation = data.frame(cell_group = unlist(lapply(unname(cell.groups.names), names)),
                               module = rep(names(cell.groups), times = sapply(cell.groups, function(x) length(x[[2]]))),
                               relation = ifelse(cell.groups.relation < 0, "negative", "positive"),
                               correlation = cell.groups.relation,
                               n_features = unlist(lapply(unname(cell.groups.names), lengths)),
                               stringsAsFactors = F, row.names = NULL)
  # Cell types of each cell group, the most frequent among its features first
  module.relation$cell_types = sapply(unlist(unname(cell.groups.names), recursive = F), function(features){
    paste(names(sort(table(unlist(lapply(features, extract_cells, cells_extra = cells_extra))), decreasing = T)), collapse = ", ")
  })

  cell.groups = list(unlist(unname(cell.groups.values), recursive = F), unlist(unname(cell.groups.names), recursive = F),  unlist(unname(cell.groups.loadings), recursive = F), module.relation)
  cell.groups[[1]] = cell.groups[[1]] %>%
    data.frame()

  ###################################################Export clusters in a table and save
  clusters = data.frame(matrix(nrow = length(cell.groups[[2]]), ncol = 2))
  for (i in 1:length(cell.groups[[2]])) {
    clusters[i,1] = names(cell.groups[[2]])[[i]]
    clusters[i,2] = paste(cell.groups[[2]][[i]], collapse ="\n")
  }
  colnames(clusters) = c("Cell groups", "Methods-signatures")
  clusters = cbind(clusters, "Module" = module.relation$module, "Relation" = module.relation$relation, "Correlation" = module.relation$correlation, "Cell types" = module.relation$cell_types)

  if(return){
    file_suffix = if(!is.null(dendrogram_file)) paste0("_", dendrogram_file)
    utils::write.csv(clusters, paste0("Results/Cell.groups.composition", file_suffix, ".csv"), row.names = F)
    utils::write.csv(cell.groups[[1]], paste0("Results/Cell.groups.scores", file_suffix, ".csv"))
  }

  return(cell.groups)
}

#' Compute associations between TF module scores and clinical metadata
#'
#' This function tests for associations between transcription factor (TF) module scores
#' and available clinical traits. Numeric traits are correlated with each module (Pearson or
#' Spearman) and shown as a heatmap. If \code{plot_grid = TRUE}, categorical traits (text, factor and
#' TRUE/FALSE columns) are additionally tested with one-way ANOVA (Tukey HSD post-hoc) and shown as
#' boxplot grids. All plots are saved in the `Results/` directory, and the tables with the results are returned.
#'
#' @param tfs.modules A numeric matrix or data frame of TF module scores across samples.
#'        Typically the first element (`TFs module matrix`) of the output from `compute.WTCNA()`. Rows represent samples, columns represent TF modules.
#' @param coldata A data frame containing clinical traits (both categorical and numerical) for the same samples.
#'        Rows must be the same samples, in the same order, as the rows of `tfs.modules` (checked by name; an
#'        error is raised otherwise). Numeric columns are treated as continuous: convert number-coded categories
#'        (e.g. a stage coded 1 to 4) to factor to test them as groups.
#' @param pval A numeric threshold (default = 0.05) to determine statistical significance (p-value, or adjusted p-value,
#'   \code{<= pval}).
#'        Only associations with (adjusted) p-values below this threshold are annotated in the heatmap or drawn as boxplots.
#' @param corr_method Character string specifying the correlation method to use for continuous variables.
#'        Options are `"p"` for pearson and `"s"` for spearman. Default is `"p"`.
#' @param file.name Character. Base file name for saving plots of results.
#' @param ncol Integer. Number of columns in the grid of association boxplots. Default is 5.
#' @param plot_grid Logical; if TRUE, tests categorical traits with ANOVA and saves boxplot grids. Default is FALSE.
#' @param p_adjust Character. Method to correct the p-values for multiple testing, applied across the modules for
#'        each trait (any method of `stats::p.adjust()`, e.g. `"BH"` or `"bonferroni"`). Default is `"none"`.
#'
#' @details
#' Traits that cannot be tested are skipped with a message: numeric traits without variation, and categorical
#' traits with a single level or with as many levels as samples (e.g. a sample identifier). Categorical traits are
#' only analysed when `plot_grid = TRUE`; a message lists them otherwise.
#'
#' The size of the figures is computed from their content (number of modules and traits for the heatmap, number of
#' panels for the boxplot grids), and the y-axis of each boxplot panel from its data.
#'
#' @return Invisibly, a list with the results:
#' \itemize{
#'   \item \code{correlation}: matrix of module-trait correlations for the numeric traits (modules x traits).
#'   \item \code{p_value}: matrix of their p-values.
#'   \item \code{p_adjusted}: matrix of their p-values adjusted with \code{p_adjust} (equal to \code{p_value} with \code{"none"}).
#'   \item \code{anova}: data frame with the ANOVA of every module for every categorical trait (\code{trait},
#'     \code{module}, \code{F}, \code{p}, \code{p_adjusted}); \code{NULL} if \code{plot_grid = FALSE}.
#' }
#' The elements are \code{NULL} when there is no trait of that type. The following files are saved to the `Results/` directory:
#' \itemize{
#'   \item \code{TF.modules_metadata_<file.name>.pdf}: a heatmap of module-trait correlations
#'     for numeric traits (only associations with p-value <= \code{pval} are annotated, with the correlation and
#'     its significance: * < 0.05, ** < 0.01, *** < 0.001, **** < 0.0001).
#'   \item \code{ANOVA_boxplot_summary_<file.name>_<trait>.svg} (if \code{plot_grid = TRUE}): one boxplot
#'     grid per categorical trait with the significant modules.
#' }
#'
#' @examples
#'
#' data("network.tuto")
#' data("traitdata.tuto")
#'
#' # The clinical traits of the tutorial data are categorical: they are tested with plot_grid = TRUE
#' res <- compute.metadata.association(
#'   tfs.modules = network.tuto[[1]],
#'   coldata = traitdata.tuto,
#'   pval = 0.05,
#'   file.name = 'Tutorial',
#'   plot_grid = TRUE
#' )
#' head(res$anova)
#'
compute.metadata.association <- function(
    tfs.modules,
    coldata,
    pval = 0.05,
    corr_method = "p",
    file.name,
    ncol = 5,
    plot_grid = F,
    p_adjust = "none"
) {
  if(!identical(make.names(rownames(tfs.modules)), make.names(rownames(coldata)))){
    stop("Samples must be the same and in the same order in tfs.modules and coldata rows")
  }

  ### Association with categorical variables (text, factor and TRUE/FALSE columns)
  coldata_categorical <- coldata %>%
    dplyr::select(dplyr::where(is.character) | dplyr::where(is.factor) | dplyr::where(is.logical))

  anova_res <- NULL
  if(ncol(coldata_categorical) != 0) {
    if(plot_grid){
      # Use the new wrapper function to generate boxplot summaries
      anova_res <- compute.metadata.association.boxplot_summary(
        tfs.modules = tfs.modules,
        coldata = coldata_categorical,
        pval = pval,
        file.name = file.name,
        ncol = ncol,
        p_adjust = p_adjust
      )
    }else{
      message(ncol(coldata_categorical), " categorical trait(s) not analysed (", paste(colnames(coldata_categorical), collapse = ", "),
              "): set plot_grid = TRUE to test them with ANOVA")
    }
  }

  ### Association with quantitative variables
  coldata_quantitative <- coldata %>%
    dplyr::select(dplyr::where(is.numeric))

  constant <- vapply(coldata_quantitative, function(x) length(unique(stats::na.omit(x))) < 2, logical(1))
  if(any(constant)){
    message("Numeric trait(s) without variation skipped: ", paste(colnames(coldata_quantitative)[constant], collapse = ", "))
    coldata_quantitative <- coldata_quantitative[, !constant, drop = FALSE]
  }

  moduleTraitCor <- NULL
  moduleTraitPvalue <- NULL
  moduleTraitPadj <- NULL

  if(ncol(coldata_quantitative) != 0){
    moduleTraitCor <- WGCNA::cor(
      tfs.modules,
      coldata_quantitative,
      method = corr_method,
      use = "pairwise.complete.obs"
    )
    nObs <- t(!is.na(tfs.modules)) %*% (!is.na(coldata_quantitative)) # Samples with both values, per module-trait pair
    moduleTraitPvalue <- WGCNA::corPvalueStudent(moduleTraitCor, nObs)

    # Multiple testing correction across the modules, for each trait
    moduleTraitPadj <- moduleTraitPvalue
    for (i in seq_len(ncol(moduleTraitPadj))) {
      moduleTraitPadj[, i] <- stats::p.adjust(moduleTraitPvalue[, i], method = p_adjust)
    }

    #### Text of the significant cells: correlation and significance label
    breaks <- c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf)
    labels <- c("****", "***", "**", "*", "")
    stars <- as.character(cut(moduleTraitPadj, breaks = breaks, labels = labels, right = FALSE))
    significant <- !is.na(moduleTraitPadj) & moduleTraitPadj <= pval

    textMatrix <- ifelse(significant, paste0(signif(moduleTraitCor, 2), ifelse(stars != "", paste0("\n", stars), "")), "")
    dim(textMatrix) <- dim(moduleTraitCor)

    p_label <- if(p_adjust == "none") "p-value" else paste0(p_adjust, " adjusted p-value")
    heatmap_title <- c(paste0("Clinical associations ", file.name),
                       paste0(if(corr_method == "s") "Spearman" else "Pearson", " correlation"),
                       paste0("shown for ", p_label, " <= ", pval),
                       "* < 0.05, ** < 0.01, *** < 0.001, **** < 0.0001")
    # The figure size is computed from the number of modules and traits (fixed cell size); the title is centred on the cells,
    # so they are made wide enough for it when there are few traits
    cell_width <- max(60, ceiling(6.3 * max(nchar(heatmap_title)) / ncol(moduleTraitCor)))
    pheatmap::pheatmap(
      moduleTraitCor,
      color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
      breaks = seq(-1, 1, length.out = 101),
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      display_numbers = textMatrix,
      number_color = "black",
      fontsize_number = 9,
      fontsize = 9,
      fontsize_row = 11,
      fontsize_col = 11,
      angle_col = 45,
      border_color = "grey85",
      cellwidth = cell_width,
      cellheight = 30,
      main = paste(heatmap_title, collapse = "\n"),
      filename = paste0("Results/TF.modules_metadata_", file.name, ".pdf")
    )
  }else{
    message("No numeric trait to correlate: the correlation heatmap is not drawn")
  }

  invisible(list(correlation = moduleTraitCor, p_value = moduleTraitPvalue, p_adjusted = moduleTraitPadj, anova = anova_res))
}

#' Compute TF module enrichment using directed target genes
#'
#' This function performs enrichment analysis for transcription factor (TF) modules
#' using known TF-target interactions. For each module, it takes the hub TFs and retrieves
#' their known targets from the TF-target network (CollecTRI by default). It then conducts an over-representation
#' analysis (ORA) against the Reactome pathway database. To reduce redundancy, only unique
#' pathways per module are retained by filtering out overlaps between modules.
#'
#' @param RNA.tpm A numeric matrix of normalized gene expression values.
#'        Rows are gene symbols, columns are samples. This matrix is used to restrict the universe
#'        of genes used in enrichment and validate TF-target gene presence.
#' @param hub_tfs A list of hub TFs per module, typically obtained using `identify_hub_TFs()`.
#'        This should be a list of two elements: one with named hub TFs per module, and another
#'        with their module membership or additional metadata. Only the modules in this list (those
#'        with hub TFs) are analysed.
#' @param TF.collection Character. The TF-target network from which the targets of the hub TFs are taken:
#'        `"CollecTRI"` (default), `"Dorothea"` (confidence A and B) or `"ARACNE"`, read as in
#'        `compute.TFs.activity()`. Use the collection the TF activity was computed with.
#' @param universe Optional. A user-specified data frame of TF-target interactions (columns `source` and
#'        `target`), used instead of `TF.collection`.
#' @param cancer.type Optional character. Cancer type used to locate the ARACNe network
#'        (`input/ARACNE/<cancer.type>/network/network.txt`), only used when `TF.collection = "ARACNE"`. If `NULL`,
#'        the network is auto-detected when only one exists under `input/ARACNE/`.
#' @param top_n Integer. Number of pathways shown in the dot plot of each module. Default is 10. All the
#'        pathways are in the returned tables.
#' @param file.name Optional character suffix for the names of the saved plots.
#'
#' @return Invisibly, a list with two elements, each a named list with one table per module (the result table
#'         of `ReactomePA::enrichPathway()`, or `NULL` for a module without enrichment):
#' \itemize{
#'   \item \code{enrichment}: the module-specific pathways (adjusted p-value < 0.05, and not significant in any other module).
#'   \item \code{enrichment_all}: all the significant pathways of each module, before removing those shared between modules.
#' }
#'         For each module with module-specific pathways, a dot plot is also saved in the `Results/` directory as a PDF file named
#'         `Enrichment_Reactome_Module_<color>_<file.name>.pdf` (`Enrichment_Reactome_Module_<color>.pdf` without `file.name`).
#'         If no enrichment is found for a module, a message is given and no file is saved for that module.
#'         An error is raised if none of the hub TFs is in the TF-target network, or if none of their targets is
#'         found in the row names of `RNA.tpm` (they must be gene symbols).
#'
#' @details
#' The targets of the hub TFs are taken from the TF-target network (for CollecTRI and Dorothea, cached in
#' `Results/TF_target_collection_<TF.collection>.csv`), and `ReactomePA::enrichPathway()` is used for ORA
#' with the Reactome database. The ORA input is the 20% most variable targets of the module's hub TFs; a module with
#' fewer than 5 such genes is skipped with a message. The background is all the genes of `RNA.tpm`. Any set of TFs has
#' targets enriched in the pathways those TFs regulate, so the results describe what the hub TFs of each module
#' regulate, not evidence that the module itself is special.
#' Pathways enriched in more than one module are removed, so only module-specific results are kept.
#' TF names are matched in their syntactic form too (e.g. `NKX2.1`, as returned by `compute.TFs.activity()`, matches `NKX2-1`).
#'
#'
#' @examples
#' \dontrun{
#' # Identify the hub TFs of a TF network from compute.WTCNA()
#' hub_tfs <- identify_hub_TFs(wgcna_network, MM_thresh = 0.8, degree_thresh = 0.9)
#'
#' # Perform module-level Reactome enrichment
#' enr <- compute.modules.enrichment(counts.norm, hub_tfs)
#' head(enr$enrichment[[1]])
#' }
compute.modules.enrichment <- function(RNA.tpm, hub_tfs, TF.collection = "CollecTRI", universe = NULL, cancer.type = NULL, top_n = 10, file.name = NULL){
  if(length(hub_tfs[[1]]) == 0){
    stop("hub_tfs contains no module with hub TFs (see identify_hub_TFs())")
  }

  #TF-target network: the collection used to compute the TF activity (as in compute.TFs.activity())
  if(!is.null(universe)){
    net = universe
  }else if(TF.collection == "ARACNE"){
    if(is.null(cancer.type)){
      # auto-discover when only one network exists
      candidates <- list.files("input/ARACNE", pattern = "^network\\.txt$", recursive = TRUE, full.names = TRUE)
      if(length(candidates) != 1)
        stop("TF.collection = 'ARACNE': specify 'cancer.type' (", length(candidates), " network.txt files found under input/ARACNE/)")
      aracne.network <- candidates[1]
    } else {
      aracne.network <- file.path("input/ARACNE", cancer.type, "network", "network.txt")
      if(!file.exists(aracne.network))
        stop("ARACNe network not found for cancer type '", cancer.type, "': ", aracne.network)
    }
    net = utils::read.table(aracne.network, header = TRUE, sep = "\t") %>%
      dplyr::select(source = Regulator, target = Target)
  }else if(TF.collection %in% c("CollecTRI", "Dorothea")){
    tf_cache_file <- paste0("Results/TF_target_collection_", TF.collection, ".csv")
    if (file.exists(tf_cache_file)) {
      net = utils::read.csv(tf_cache_file, row.names = 1)
    } else if(TF.collection == "Dorothea"){
      net = dplyr::filter(dorothea::dorothea_hs, .data$confidence %in% c("A", "B")) %>%
        dplyr::mutate(source = .data$tf) %>%
        dplyr::select(-tf)
      utils::write.csv(net, tf_cache_file)
    } else {
      net = decoupleR::get_collectri(organism = 'human', split_complexes = F)
      utils::write.csv(net, tf_cache_file)
    }
  }else{
    stop("TF.collection must be 'CollecTRI', 'Dorothea' or 'ARACNE' (or give the TF-target network in 'universe')")
  }

  if(!all(c("source", "target") %in% colnames(net))){
    stop("The TF-target network must have the columns 'source' (TF) and 'target' (gene)")
  }

  #Hub TFs found in the network and target genes found in the expression matrix, per module
  tfs_found = sapply(hub_tfs[[1]], function(tfs) sum(make.names(tfs) %in% make.names(net$source)))
  targets_found = sapply(hub_tfs[[1]], function(tfs) sum(rownames(RNA.tpm) %in% net$target[make.names(net$source) %in% make.names(tfs)]))
  if(all(tfs_found == 0)){
    stop("None of the hub TFs is in the TF-target network: use the TF.collection (or universe) the TF activity was computed with")
  }
  if(all(targets_found == 0)){
    stop("No target gene of the hub TFs was found in the row names of RNA.tpm: they must be gene symbols")
  }
  no_targets = names(targets_found)[targets_found == 0]
  if(length(no_targets) > 0){
    message("No target gene of the hub TFs found in RNA.tpm for module(s): ", paste(no_targets, collapse = ", "))
  }

  message("Enrichment computed for the ", length(hub_tfs[[1]]), " module(s) with hub TFs: ", paste(names(hub_tfs[[1]]), collapse = ", "))

  res = list()
  pathways = list()

  #Pathway enrichment using target genes
  for (i in 1:length(hub_tfs[[1]])) {
    color = names(hub_tfs[[1]])[i]
    res[i] = list(module_enrich(as.matrix(RNA.tpm), color, hub_tfs, net)) #list() keeps modules without enrichment (NULL)
    names(res)[i] = color
    pathways[[i]] = if(is.null(res[[i]])) character(0) else res[[i]]@result[["Description"]]
    names(pathways)[i] = color
  }

  #All the significant pathways of each module (before removing those shared between modules)
  enrichment_all = lapply(res, function(x) if(is.null(x)) NULL else x@result)

  if(length(unlist(pathways)) == 0){ #gplots::venn() fails when all sets are empty
    message("No enrichment for any module")
    return(invisible(list(enrichment = enrichment_all, enrichment_all = enrichment_all)))
  }

  ItemsList <- gplots::venn(pathways, show.plot = FALSE)
  x = attributes(ItemsList)

  #Keep only unique pathways
  for (i in 1:length(hub_tfs[[1]])) {
    color = names(hub_tfs[[1]])[i]
    if(is.null(res[[i]]) == T){
      if(!color %in% no_targets) message("No enrichment for module ", color)
    }else{
      res[[i]]@result = res[[i]]@result[which(res[[i]]@result[["Description"]]%in%x$intersections[[color]]),] #take unique pathways per module
      if(nrow(res[[i]]@result)==0){
        message("No enrichment for module ", color)
        res[i] = list(NULL)
      }else{
        pdf(paste0("Results/Enrichment_Reactome_Module_", color, if(!is.null(file.name)) paste0("_", file.name), ".pdf"))
        suppressMessages(print(enrichplot::dotplot(res[[i]], showCategory = top_n, title=paste0("Enrichment by Reactome\nModule ", color))))
        dev.off()
      }
    }
  }

  invisible(list(enrichment = lapply(res, function(x) if(is.null(x)) NULL else x@result), enrichment_all = enrichment_all))
}

#' Compute modules relationship
#'
#' Correlates every column of one matrix with every column of another (e.g., TF modules and external features such as
#' deconvolution estimates or pathway activities) across the same set of samples, and visualizes the significant associations.
#' It is meant for continuous variables only: all the columns of both matrices must be numeric. For categorical
#' variables (e.g. clinical groups), use \code{compute.metadata.association()}.
#'
#' @param matA A numeric matrix or data frame of continuous features (samples x features), e.g. TF module scores.
#' @param matB A numeric matrix or data frame of continuous features to correlate with (samples x features).
#' @param file_name A string indicating the base name (without extension) of the figures to be saved in the "Results/" folder.
#'   If \code{NULL} (default), no figure is saved.
#' @param batch Optional vector indicating batch assignment for samples (no missing values). If it has two or more
#'   levels, partial correlations (\code{ppcor::pcor.test()}, Pearson or Spearman according to \code{cor_type}) are
#'   computed controlling for batch. Batch is treated as categorical (one dummy variable per level beyond the first).
#' @param pval A numeric threshold for the p-value to define significance (p-value <= \code{pval}). Default is 0.05.
#' @param cor_type Type of correlation to compute: "p" (Pearson) or "s" (Spearman). Default is "p". It also sets
#'   the type of the partial correlations computed with \code{batch}.
#' @param vertical Logical; if TRUE, modules are on the x-axis and the \code{matB} features on the y-axis.
#'   Otherwise (default), the \code{matB} features are on the x-axis and modules on the y-axis.
#' @param plot Logical; if TRUE, saves the heatmap plot as a PDF (\code{file_name} is needed). Default is TRUE.
#' @param plot.grid Logical; if TRUE, generates per-pair scatter grid plots for significant associations, saved to
#'   "Results/<file_name>_scatter_grid.svg" (the 100 strongest pairs at most). With Spearman, the ranks of the values
#'   are shown, as the correlation is computed on them. With \code{batch}, the values shown are
#'   corrected for batch, as the correlations are. Default is FALSE.
#' @param ncol.grid Integer number of columns used in scatter grid layout. If NULL (default), it is set to
#'   the ceiling of the square root of the number of panels.
#'
#' @return Invisibly, a list with:
#' \itemize{
#'   \item \code{correlation}: the correlation matrix (columns of \code{matA} as rows, columns of \code{matB} as columns).
#'   \item \code{significant}: a named list with, for each column of \code{matA}, the columns of \code{matB} significantly
#'     correlated with it (p-value <= \code{pval}).
#'   \item \code{p_value}: the matrix of p-values.
#' }
#' If `plot = TRUE` and `file_name` is given, a heatmap of the correlations is saved to
#' "Results/<file_name>.pdf" (the ".pdf" extension is added if missing). Its size is computed from the number of
#' features. With up to 40 features per axis, the significant cells show the correlation and its significance
#' (* < 0.05, ** < 0.01, *** < 0.001, **** < 0.0001); with more, they are marked with an asterisk.
#'
#'
#' @details
#' `matA` and `matB` must have the same samples in the same order (same rownames);
#' names are compared after `make.names()`, so e.g. "TCGA-XX" and "TCGA.XX" match. Correlations and p-values use,
#' for each pair of columns, the samples with a value in both (missing values are allowed); a pair with fewer
#' than three such samples, or without variation among them, gets \code{NA}. Columns without variation are skipped with a message. Column names
#' must be unique within each matrix.
#'
#' @examples
#' data("network.tuto")
#' data("deconv_subgroups.tuto")
#' corr = compute.modules.relationship(network.tuto[[1]],
#'                                     deconv_subgroups.tuto[[1]],
#'                                     plot = FALSE,
#'                                     pval = 0.01)
#' corr$correlation[1:3, 1:3]
#' corr$significant[[1]]
#'
#' \dontrun{
#' # TF modules vs PROGENy pathways (downloads the PROGENy model from OmniPath)
#' data("counts.norm.tuto")
#' pathways <- compute.pathway.activity(counts.norm.tuto)
#' compute.modules.relationship(network.tuto[[1]],
#'                              pathways,
#'                              "Pathways_Progeny-TFs_Modules")
#' }
#'
compute.modules.relationship <- function(matA, matB, file_name = NULL, batch = NULL, pval = 0.05, cor_type = "p", vertical = F, plot = T, plot.grid = F, ncol.grid = NULL){

  matA = data.frame(matA, check.names = FALSE)
  matB = data.frame(matB, check.names = FALSE) # keep feature names unchanged (e.g. "-") so they match the deconvolution columns

  if(!cor_type %in% c("p", "s")){
    stop("cor_type must be 'p' (Pearson) or 's' (Spearman)")
  }
  if(!is.null(batch) && anyNA(batch)){
    stop("batch has missing values: every sample needs a batch")
  }
  duplicated_names = unique(c(colnames(matA)[duplicated(colnames(matA))], colnames(matB)[duplicated(colnames(matB))]))
  if(length(duplicated_names) > 0){
    stop("Column names must be unique within matA and within matB. Duplicated: ", paste(utils::head(duplicated_names, 5), collapse = ", "))
  }
  not_numeric = c(colnames(matA)[!sapply(matA, is.numeric)], colnames(matB)[!sapply(matB, is.numeric)])
  if(length(not_numeric) > 0){
    stop("compute.modules.relationship() is for continuous variables only. Columns that are not numeric: ", paste(utils::head(not_numeric, 5), collapse = ", "),
         if(length(not_numeric) > 5) paste0(", ... (", length(not_numeric), " in total)"), ". Use compute.metadata.association() for categorical variables")
  }
  if(nrow(matA) != nrow(matB)){
    stop("matA and matB must have the same samples: they have ", nrow(matA), " and ", nrow(matB), " rows")
  }
  if(length(rownames(matA)) == 0 || !all(make.names(rownames(matA)) == make.names(rownames(matB)))){
    stop("matA and matB must have the same samples in the same order: their row names differ")
  }

  # Columns without variation cannot be correlated
  constantA = sapply(matA, function(x) length(unique(stats::na.omit(x))) < 2)
  constantB = sapply(matB, function(x) length(unique(stats::na.omit(x))) < 2)
  if(any(constantA) || any(constantB)){
    message("Columns without variation skipped: ", paste(c(colnames(matA)[constantA], colnames(matB)[constantB]), collapse = ", "))
    matA = matA[, !constantA, drop = FALSE]
    matB = matB[, !constantB, drop = FALSE]
  }
  if(ncol(matA) == 0 || ncol(matB) == 0){
    stop("No column with variation left to correlate")
  }

  # Values shown in the scatter plots: for Spearman the ranks, which the correlation is computed on
  matA_plot = matA
  matB_plot = matB
  if(cor_type == "s"){
    matA_plot[] = lapply(matA, rank, na.last = "keep")
    matB_plot[] = lapply(matB, rank, na.last = "keep")
  }
  cor_method = if(cor_type == "s") "spearman" else "pearson"

  # ---------- PARTIAL CORRELATION IF BATCH PROVIDED ----------
  if(!is.null(batch) && length(unique(batch)) > 1){
    if(length(batch) != nrow(matA)) stop("Length of batch must match number of samples")
    cohort = factor(batch)

    # Batch is categorical: one dummy variable per batch level beyond the first
    batch <- stats::model.matrix(~ cohort)[, -1, drop = FALSE]

    moduleTraitCor <- matrix(NA, nrow = ncol(matA), ncol = ncol(matB))
    moduleTraitPvalue <- matrix(NA, nrow = ncol(matA), ncol = ncol(matB))
    rownames(moduleTraitCor) <- colnames(matA)
    colnames(moduleTraitCor) <- colnames(matB)
    rownames(moduleTraitPvalue) <- colnames(matA)
    colnames(moduleTraitPvalue) <- colnames(matB)

    for(i in 1:ncol(matA)){
      for(j in 1:ncol(matB)){
        ok <- !is.na(matA[,i]) & !is.na(matB[,j]) # Samples with both values
        if(sum(ok) < 3 || length(unique(matA[ok,i])) < 2 || length(unique(matB[ok,j])) < 2) next # Not enough samples, or no variation among them: left as NA
        covar <- batch[ok, , drop = FALSE]
        covar <- covar[, apply(covar, 2, function(x) length(unique(x)) > 1), drop = FALSE] # Batches still present in these samples
        pc <- if(ncol(covar) > 0) ppcor::pcor.test(matA[ok,i], matB[ok,j], covar, method = cor_method) else stats::cor.test(matA[ok,i], matB[ok,j], method = cor_method, exact = FALSE)
        moduleTraitCor[i,j] <- pc$estimate
        moduleTraitPvalue[i,j] <- pc$p.value
      }
    }

    # Values corrected for batch, shown in the scatter plots
    matA_plot = as.data.frame(lapply(matA_plot, function(x) x - stats::ave(x, cohort, FUN = function(v) mean(v, na.rm = TRUE))), check.names = FALSE, row.names = rownames(matA))
    matB_plot = as.data.frame(lapply(matB_plot, function(x) x - stats::ave(x, cohort, FUN = function(v) mean(v, na.rm = TRUE))), check.names = FALSE, row.names = rownames(matB))
    cor_label = paste0(if(cor_type == "s") "Spearman" else "Pearson", " partial correlation (corrected for batch)")
  } else {
    # ---------- STANDARD CORRELATION (samples with both values, for each pair of columns) ----------
    moduleTraitCor = WGCNA::cor(matA, matB, method = cor_type, use = "pairwise.complete.obs")
    nObs = t(!is.na(matA)) %*% (!is.na(matB))
    moduleTraitPvalue = WGCNA::corPvalueStudent(moduleTraitCor, nObs)
    dimnames(moduleTraitPvalue) = dimnames(moduleTraitCor)
    cor_label = paste0(if(cor_type == "s") "Spearman" else "Pearson", " correlation")
  }

  #Significant features of matB for each column of matA
  sig = lapply(seq_len(nrow(moduleTraitPvalue)), function(i) colnames(moduleTraitPvalue)[which(moduleTraitPvalue[i,] <= pval)])
  names(sig) = rownames(moduleTraitPvalue)

  if((plot || plot.grid) && is.null(file_name)){
    message("No file_name given: the figures are not saved")
  }

  if (plot.grid && !is.null(file_name)) {
    plot_module_scatter_grid(
      matA = matA_plot,
      matB = matB_plot,
      cor_mat = moduleTraitCor,
      p_mat = moduleTraitPvalue,
      file_name = file_name,
      pval = pval,
      ncol = ncol.grid,
      only_sig = TRUE,
      axis_suffix = if(cor_type == "s") " (rank)" else ""
    )
  }

  ###Heatmap of the correlations (size computed from the number of features)
  if(plot && !is.null(file_name)){
    significant = !is.na(moduleTraitPvalue) & moduleTraitPvalue <= pval
    stars = as.character(cut(moduleTraitPvalue, breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf), labels = c("****", "***", "**", "*", ""), right = FALSE))
    show_numbers = ncol(matA) <= 40 && ncol(matB) <= 40 # Otherwise the cells are too small for text: significant cells are marked with an asterisk
    textMatrix = if(show_numbers) ifelse(significant, paste0(signif(moduleTraitCor, 2), ifelse(stars != "", paste0("\n", stars), "")), "") else ifelse(significant, "*", "")
    dim(textMatrix) = dim(moduleTraitCor)
    dimnames(textMatrix) = dimnames(moduleTraitCor)

    #matB features ordered by the similarity of their correlations
    feature_order = if(ncol(moduleTraitCor) > 2 && !anyNA(moduleTraitCor)) stats::hclust(stats::dist(t(moduleTraitCor), method = "euclidean"), method = "ward.D2")$order else seq_len(ncol(moduleTraitCor))
    heat = moduleTraitCor[, feature_order, drop = FALSE]
    heat_text = textMatrix[, feature_order, drop = FALSE]
    if(vertical){
      heat = t(heat)
      heat_text = t(heat_text)
    }

    heatmap_title = c(cor_label, paste0("shown for p-value <= ", pval), if(show_numbers) "* < 0.05, ** < 0.01, *** < 0.001, **** < 0.0001")
    cell_width = if(show_numbers) 36 else if(vertical) 14 else 11
    cell_width = max(cell_width, ceiling(6.3 * max(nchar(heatmap_title)) / ncol(heat))) # The title is centred on the cells: wide enough for it with few columns
    pheatmap::pheatmap(
      heat,
      color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
      breaks = seq(-1, 1, length.out = 101),
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      display_numbers = heat_text,
      number_color = "black",
      fontsize_number = if(show_numbers) 7.5 else 7,
      fontsize = 9,
      fontsize_row = if(nrow(heat) > 40) 7 else 10,
      fontsize_col = if(ncol(heat) > 40) 7 else 10,
      angle_col = 90, # vertical labels: long feature names stay inside the page
      border_color = "grey85",
      na_col = "grey90",
      cellwidth = cell_width,
      cellheight = if(show_numbers) 26 else if(vertical) 11 else 14,
      main = paste(heatmap_title, collapse = "\n"),
      filename = paste0("Results/", file_name, if (!grepl("\\.pdf$", file_name)) ".pdf")
    )
  }

  invisible(list(correlation = moduleTraitCor, significant = sig, p_value = moduleTraitPvalue))
}

#' Compute pathway activity scores
#'
#' This function computes pathway activity scores from log-normalized gene expression data with \code{decoupleR}
#' (Badia-i-Mompel et al., 2022), based on the PROGENy resource (Schubert et al., 2018). By default the univariate
#' linear model (\code{"ulm"}) is used.
#' Optionally, it also performs Gene Set Variation Analysis (GSVA) using hallmark signatures or any user-provided gene sets.
#'
#' @param RNA.tpm A numeric matrix of log-normalized gene expression values (e.g. log2(TPM + 1)) with genes (symbols) as
#'   rows and samples as columns. A warning is given if it has values above 20, which log2(TPM + 1) cannot reach.
#' @param gene_sets A named list of gene sets (character vectors of gene symbols; e.g., hallmark signatures or user-defined
#'   sets). If provided, GSVA scores will be computed for these sets. Default is \code{NULL}.
#' @param paths A data frame describing the pathway-gene interactions (columns \code{source}, \code{target}, \code{weight}) for use with PROGENy.
#'   If \code{NULL}, the human PROGENy resource (top 500 genes) is used and cached in \code{Results/Pathways_collection_PROGENy.csv}
#'   (also when \code{return = FALSE}, so that it is downloaded only once). A message gives the number of pathways and
#'   interactions of the cached file and the date it was saved; a cached file that is empty or lacks the needed columns
#'   is not used (the collection is downloaded again, with a warning). Delete the file to download the collection again.
#' @param method Character. The \code{decoupleR} statistic used to score the PROGENy pathways: \code{"ulm"} (univariate
#'   linear model, default), \code{"mlm"} (multivariate linear model), \code{"wsum"} (weighted sum) or \code{"consensus"}
#'   (consensus of the three). The gene sets of \code{gene_sets} are always scored with GSVA.
#' @param scale Logical. If TRUE (default), each pathway is z-scored across the samples. The scaled scores of a
#'   sample then depend on the other samples given in the same call; use \code{scale = FALSE} to get scores that can be
#'   compared between separate calls (e.g. a training and a test set). With fewer than three samples the scores are
#'   not scaled (a message is given).
#' @param return Logical; if TRUE, saves the score matrices in Results/ folder. Default is TRUE.
#' @param file.name Optional character suffix used when writing output CSV files
#'   (\code{Results/Pathway_matrix_PROGENy_<file.name>.csv} and \code{Results/Pathway_matrix_GSVA_<file.name>.csv};
#'   without it, \code{Results/Pathway_matrix_PROGENy.csv} and \code{Results/Pathway_matrix_GSVA.csv}).
#'
#' @return If \code{gene_sets} is \code{NULL}, a data frame of PROGENy pathway activity scores (samples as rows, pathways as columns),
#' scaled if \code{scale = TRUE}.
#' If \code{gene_sets} is provided, a list with two elements:
#' \itemize{
#'   \item \code{PROGENy}: A data frame of PROGENy pathway activity scores.
#'   \item \code{GSVA}: A data frame of GSVA scores based on the provided gene sets.
#' }
#' Column names are made syntactically valid with \code{make.names()} (e.g. \code{JAK.STAT}), in the returned tables and
#' in the saved files. Pathways with fewer than 5 of their genes in the data, and gene sets that GSVA cannot score, are
#' left out with a message.
#'
#'
#' @references
#' Schubert M, Klinger B, Kluenemann M, Sieber A, Uhlitz F, Sauer S, Garnett MJ, Bluethgen N, Saez-Rodriguez J.
#' Perturbation-response genes reveal signaling footprints in cancer gene expression. Nature Communications. 2018. \doi{10.1038/s41467-017-02391-6}
#'
#' Badia-i-Mompel, P. et al. (2022). decoupleR: ensemble of computational methods to infer biological activities from omics data. *Bioinformatics Advances*, 2(1), vbac016. https://doi.org/10.1093/bioadv/vbac016
#'
#' @examples
#' \dontrun{
#' # Compute only PROGENy activities (downloads the PROGENy model from OmniPath)
#' data("counts.norm.tuto")
#' pathways <- compute.pathway.activity(counts.norm.tuto)
#' }
#'
compute.pathway.activity <- function(RNA.tpm, gene_sets = NULL, paths = NULL, method = "ulm", scale = TRUE, return = TRUE, file.name = NULL) {

  if (!method %in% c("ulm", "mlm", "wsum", "consensus")) {
    stop("method must be 'ulm', 'mlm', 'wsum' or 'consensus'")
  }
  if (!is.null(gene_sets) && (!is.list(gene_sets) || is.null(names(gene_sets)) || any(is.na(names(gene_sets)) | names(gene_sets) == "") || !all(sapply(gene_sets, is.character)))) {
    stop("gene_sets must be a named list of character vectors (the genes of each gene set)")
  }

  rn <- rownames(RNA.tpm)
  RNA.tpm <- apply(as.matrix(RNA.tpm), 2, as.numeric)
  rownames(RNA.tpm) <- rn
  results_list <- list()

  if (max(RNA.tpm, na.rm = TRUE) > 20) { # log2(TPM + 1) cannot exceed 20
    warning("RNA.tpm has values above 20, so it does not look log-transformed: log-normalized expression (e.g. log2(TPM + 1)) is expected")
  }
  if (scale && ncol(RNA.tpm) < 3) {
    message("Fewer than three samples: the pathway scores are not scaled")
    scale <- FALSE
  }

  ###### PROGENy
  progeny_cache_file <- "Results/Pathways_collection_PROGENy.csv"
  # Cached collection: checked before it is used
  if (is.null(paths) && file.exists(progeny_cache_file)) {
    paths <- tryCatch(utils::read.csv(progeny_cache_file, row.names = 1), error = function(e) NULL)
    if (is.null(paths) || nrow(paths) == 0 || !all(c("source", "target", "weight") %in% colnames(paths))) {
      warning("The cached PROGENy collection ", progeny_cache_file, " cannot be used (it is empty or lacks the columns 'source', 'target' and 'weight'): the collection is downloaded again")
      paths <- NULL
    } else {
      message("Using cached PROGENy pathways collection from ", progeny_cache_file, " (", length(unique(paths$source)), " pathways, ", nrow(paths), " interactions, saved on ",
              format(file.mtime(progeny_cache_file), "%Y-%m-%d"), "). Delete this file to download the collection again")
    }
  }
  if (is.null(paths)) {
    paths <- decoupleR::get_progeny(organism = "human", top = 500)
    if (!is.data.frame(paths) || nrow(paths) == 0 || !all(c("source", "target", "weight") %in% colnames(paths))) {
      stop("The PROGENy collection could not be downloaded (empty, or without the columns 'source', 'target' and 'weight')")
    }
    # Written to a temporary file first, so that an interrupted or simultaneous run never leaves or reads an incomplete cache
    utils::write.csv(paths, paste0(progeny_cache_file, ".tmp"))
    file.rename(paste0(progeny_cache_file, ".tmp"), progeny_cache_file)
  }
  if (!all(c("source", "target", "weight") %in% colnames(paths))) {
    stop("paths must have the columns 'source', 'target' and 'weight'")
  }

  progeny <- decoupleR::decouple(
    mat             = RNA.tpm,
    network         = dplyr::mutate(as.data.frame(paths), mor = .data$weight),
    .source         = "source",
    .target         = "target",
    statistics      = if (method == "consensus") NULL else method,
    consensus_score = (method == "consensus"),
    minsize         = 5
  ) %>%
    dplyr::filter(.data$statistic == method)

  sample_acts_progeny <- progeny %>%
    tidyr::pivot_wider(id_cols = "condition", names_from = "source", values_from = "score") %>%
    tibble::column_to_rownames("condition") %>%
    as.matrix()
  sample_acts_progeny <- sample_acts_progeny[colnames(RNA.tpm), , drop = FALSE] # samples in the order of the input

  not_scored <- setdiff(unique(paths$source), colnames(sample_acts_progeny))
  if (length(not_scored) > 0) {
    message("Pathways left out (fewer than 5 of their genes in the data): ", paste(not_scored, collapse = ", "))
  }
  if (scale) sample_acts_progeny <- base::scale(sample_acts_progeny)

  results_list$PROGENy <- as.data.frame(sample_acts_progeny)

  ###### GSVA (optional)
  if (!is.null(gene_sets)) {
    gsva_results <- GSVA::gsva(
      GSVA::gsvaParam(RNA.tpm, gene_sets, kcdf = "Gaussian", minSize = 1, maxDiff = TRUE),
      verbose = FALSE
    )

    not_scored <- setdiff(names(gene_sets), rownames(gsva_results))
    if (length(not_scored) > 0) {
      message("Gene sets left out by GSVA (not enough of their genes in the data): ", paste(not_scored, collapse = ", "))
    }

    sample_acts_gsva <- t(gsva_results)
    if (scale) sample_acts_gsva <- base::scale(sample_acts_gsva)

    results_list$GSVA <- as.data.frame(sample_acts_gsva)
  }

  # standardise column names so they are valid R identifiers (consistent with caret internals)
  results_list <- lapply(results_list, function(df) {
    colnames(df) <- make.names(colnames(df))
    df
  })

  ###### Save outputs if requested
  if (return) {
    file_suffix <- if (!is.null(file.name)) paste0("_", file.name) else ""
    if (!is.null(results_list$PROGENy)) {
      utils::write.csv(results_list$PROGENy, paste0("Results/Pathway_matrix_PROGENy", file_suffix, ".csv"))
    }
    if (!is.null(results_list$GSVA)) {
      utils::write.csv(results_list$GSVA, paste0("Results/Pathway_matrix_GSVA", file_suffix, ".csv"))
    }
  }

  ###### RETURN FORMAT LOGIC
  # If only one element -> return it directly (matrix)
  if (length(results_list) == 1) {
    return(results_list[[1]])
  }

  return(results_list)
}

#' Compute Transcription Factor (TF) activity
#'
#' Infers transcription factor (TF) activity from a gene expression matrix with \code{decoupleR::decouple()}
#' (Badia-i-Mompel et al., 2022). By default the univariate linear model (\code{"ulm"}) is used.
#' The TF-target network can be provided by the user, obtained from prebuilt collections (CollecTRI from
#' OmnipathR or Dorothea from the \code{dorothea} package), or read from an ARACNe-inferred network.
#'
#' @param RNA.counts A gene expression matrix with genes (symbols) as rows and samples as columns, with unique sample
#'   names. The matrix should be log-normalized (e.g., log2(TPM + 1), log2CPM, etc.). A warning is given if it has
#'   values above 20, which log2(TPM + 1) cannot reach.
#' @param TF.collection Character. The source of the TF-target network. Options are `"CollecTRI"` (default), `"Dorothea"`, or `"ARACNE"`.
#' - `"CollecTRI"` uses the prebuilt collection from OmnipathR and `"Dorothea"` (confidence A and B) the one from the
#'   `dorothea` package. Each collection is cached in its own file,
#'   `Results/TF_target_collection_<TF.collection>.csv`, and reused on later calls: a message gives the number of
#'   TFs and interactions of the cached file and the date it was saved. A cached file that is empty or lacks the
#'   columns `source`, `target` and `mor` is not used (the collection is obtained again, with a warning). Delete the
#'   file to obtain the collection again.
#' - `"ARACNE"` reads a tab-separated network file with `Regulator` and `Target` columns from
#'   `input/ARACNE/<cancer.type>/network/network.txt` (relative to the working directory). The mode of
#'   regulation of each edge is read from a `mor` column of the file (1 activation, -1 repression), which should be
#'   the sign of the Spearman correlation between TF and target computed on the cohort used to build the network.
#'   If the file has no `mor` column, a warning is given and the signs are computed as the sign of the Spearman
#'   correlation between TF and target in the samples of `RNA.counts`; they then depend on the samples given (e.g.
#'   they can differ between cohorts, or between a training and a test set).
#' @param min_targets_size Integer. Minimum number of target genes per regulon (passed to
#'   \code{decoupleR::decouple()} as \code{minsize}). Default is 5.
#' @param universe Optional. A user-specified data frame of TF-target interactions (columns \code{source},
#'   \code{target}, \code{mor}). The three columns are required: \code{mor} is the mode of regulation
#'   (1 activation, -1 repression; use 1 for every edge if the network has no sign). If not provided, the
#'   network is fetched based on `TF.collection`. Ignored when \code{TF.collection = "ARACNE"}.
#' @param cancer.type Optional character. TCGA cancer type abbreviation used to locate the ARACNe network
#'   (only used when \code{TF.collection = "ARACNE"}). If \code{NULL}, the network is auto-detected when only one
#'   \code{network.txt} exists under \code{input/ARACNE/}.
#' @param method Character. The \code{decoupleR} statistic used to score the TFs: \code{"ulm"} (univariate linear
#'   model, default), \code{"mlm"} (multivariate linear model), \code{"wsum"} (weighted sum) or \code{"consensus"}
#'   (consensus of the three).
#' @param scale Logical. If TRUE (default), z-score scales the TF activity matrix across samples. The scaled scores
#'   of a sample then depend on the other samples given in the same call; use \code{scale = FALSE} to get scores that
#'   can be compared between separate calls (e.g. a training and a test set). With fewer than three samples the
#'   scores are not scaled (a message is given).
#' @param return Logical; if TRUE, saves matrix in Results/ folder. Default is TRUE.
#' @param file.name Optional character suffix used when writing the TF activity matrix to disk
#'   (\code{Results/TF_matrix_<file.name>.csv}; without it, \code{Results/TF_matrix.csv}).
#'
#' @return A data frame of inferred (and, if \code{scale = TRUE}, scaled) TF activity scores, with samples as rows
#'   and TFs as columns. Column names are made syntactically valid with \code{make.names()} (e.g. \code{NKX2.1}),
#'   in the returned table and in the saved file.
#'
#' @references
#' Badia-i-Mompel, P. et al. (2022). decoupleR: ensemble of computational methods to infer biological activities from omics data. *Bioinformatics Advances*, 2(1), vbac016. https://doi.org/10.1093/bioadv/vbac016
#'
#' Tuerei, D., Korcsmaros, T., & Saez-Rodriguez, J. (2016). OmniPath: guidelines and gateway for literature-curated signaling pathway resources. *Nature Methods*, 13(12), 966-967. https://doi.org/10.1038/nmeth.4077
#'
#' Garcia-Alonso, L. et al. (2019). Benchmark and integration of resources for the estimation of human transcription factor activities. *Genome Research*. https://doi.org/10.1101/gr.240663.118
#'
#' Lachmann, A. et al. (2016). ARACNe-AP: gene network reverse engineering through adaptive partitioning inference of mutual information. *Bioinformatics*, 32(14), 2233-2235. https://doi.org/10.1093/bioinformatics/btw216
#'
#' Margolin, A.A. et al. (2006). ARACNE: an algorithm for the reconstruction of gene regulatory networks in a mammalian cellular context. *BMC Bioinformatics*, 7(Suppl 1), S7. https://doi.org/10.1186/1471-2105-7-S1-S7
#'
#' @examples
#' \dontrun{
#' # Downloads the CollecTRI network from OmniPath
#' data("counts.norm.tuto")
#' tfs_activity <- compute.TFs.activity(counts.norm.tuto)
#' }
#'
compute.TFs.activity <- function(RNA.counts, TF.collection = "CollecTRI", min_targets_size = 5, universe = NULL, cancer.type = NULL, method = "ulm", scale = TRUE, return = TRUE, file.name = NULL){

  if(!method %in% c("ulm", "mlm", "wsum", "consensus")){
    stop("method must be 'ulm', 'mlm', 'wsum' or 'consensus'")
  }
  if(!TF.collection %in% c("CollecTRI", "Dorothea", "ARACNE") && is.null(universe)){
    stop("TF.collection must be 'CollecTRI', 'Dorothea' or 'ARACNE'")
  }
  if(!is.null(universe) && TF.collection != "ARACNE" && !all(c("source", "target", "mor") %in% colnames(universe))){
    stop("universe must have the columns 'source', 'target' and 'mor' (mode of regulation: 1 activation, -1 repression; use 1 for every edge if the network has no sign)")
  }
  if(is.null(colnames(RNA.counts)) || anyDuplicated(colnames(RNA.counts)) > 0){
    stop("RNA.counts must have unique sample names as column names")
  }
  if(max(RNA.counts, na.rm = TRUE) > 20){ # log2(TPM + 1) cannot exceed 20
    warning("RNA.counts has values above 20, so it does not look log-transformed: log-normalized expression (e.g. log2(TPM + 1)) is expected")
  }
  if(scale && ncol(RNA.counts) < 3){
    message("Fewer than three samples: the TF activity scores are not scaled")
    scale <- FALSE
  }

  tf_cache_file <- paste0("Results/TF_target_collection_", TF.collection, ".csv")

  if(TF.collection == "ARACNE"){

    if(is.null(cancer.type)){
      # auto-discover when only one network exists
      candidates <- list.files("input/ARACNE", pattern = "^network\\.txt$",
                               recursive = TRUE, full.names = TRUE)
      if(length(candidates) == 0)
        stop("TF.collection = 'ARACNE' requires a 'cancer.type' or a network.txt under input/ARACNE/")
      if(length(candidates) > 1)
        stop("Multiple ARACNe networks found. Specify 'cancer.type' (e.g. cancer.type = 'skcm'):\n",
             paste(dirname(dirname(candidates)), collapse = "\n"))
      aracne.network <- candidates[1]
      message("Auto-detected ARACNe network: ", aracne.network)
    } else {
      aracne.network <- file.path("input/ARACNE", cancer.type, "network", "network.txt")
      if(!file.exists(aracne.network))
        stop("ARACNe network not found for cancer type '", cancer.type, "': ", aracne.network)
    }

    message("Loading ARACNe network from: ", aracne.network)
    # Read network edges, filter to genes present in expression matrix
    aracne_net <- utils::read.table(aracne.network, header = TRUE, sep = "\t")
    if(!all(c("Regulator", "Target") %in% colnames(aracne_net))){
      stop("The ARACNe network file must have the columns 'Regulator' and 'Target': ", aracne.network)
    }
    if("mor" %in% colnames(aracne_net)){
      # Signs given in the network file (computed on the cohort used to build the network)
      if(!is.numeric(aracne_net$mor)){
        stop("The 'mor' column of the ARACNe network file must be numeric (1 activation, -1 repression): ", aracne.network)
      }
      universe <- aracne_net %>%
        dplyr::select(source = Regulator, target = Target, mor) %>%
        dplyr::filter(target %in% rownames(RNA.counts) & !is.na(mor) & mor != 0)
    } else {
      warning("The ARACNe network file has no 'mor' column (sign of each interaction): the signs are computed as the sign of the Spearman correlation ",
              "between each TF and its target in the samples given, so they depend on these samples. Add a 'mor' column (1 activation, -1 repression) ",
              "computed on the cohort used to build the network: ", aracne.network)
      aracne_net <- aracne_net %>%
        dplyr::select(source = Regulator, target = Target) %>%
        dplyr::filter(source %in% rownames(RNA.counts) & target %in% rownames(RNA.counts))

      # Compute Spearman correlation between every TF and its targets in one matrix op
      # (this is exactly what TFmode1 does internally)
      all_tfs     <- unique(aracne_net$source)
      all_targets <- unique(aracne_net$target)
      cor_mat <- suppressWarnings(
        stats::cor(t(RNA.counts[all_tfs, , drop = FALSE]),
                   t(RNA.counts[all_targets, , drop = FALSE]),
                   method = "spearman")
      )

      # mor = sign of Spearman correlation - +1 activation, -1 repression
      universe <- aracne_net %>%
        dplyr::mutate(mor = cor_mat[cbind(source, target)]) %>%
        dplyr::mutate(mor = sign(mor)) %>%
        dplyr::filter(!is.na(mor) & mor != 0)
    }

    message("Computing TF activities...")

    sample_acts <- decoupleR::decouple( mat     = RNA.counts,
                                        network = universe,
                                        .source = "source",
                                        .target = "target",
                                        statistics = if(method == "consensus") NULL else method,
                                        consensus_score = (method == "consensus"),
                                        minsize = min_targets_size
                                      ) %>%
      dplyr::filter(.data$statistic == method) %>%
      decoupleR::pivot_wider_profile(id_cols     = source,
                                     names_from  = condition,
                                     values_from = score) %>%
      as.matrix() %>%
      t()

  } else {

    # Cached collection: checked before it is used
    if(is.null(universe) && file.exists(tf_cache_file)){
      universe = tryCatch(utils::read.csv(tf_cache_file, row.names = 1), error = function(e) NULL)
      if(is.null(universe) || nrow(universe) == 0 || !all(c("source", "target", "mor") %in% colnames(universe))){
        warning("The cached TF-target collection ", tf_cache_file, " cannot be used (it is empty or lacks the columns 'source', 'target' and 'mor'): the collection is obtained again")
        universe = NULL
      } else {
        message("Using cached TF-target collection from ", tf_cache_file, " (", length(unique(universe$source)), " TFs, ", nrow(universe), " interactions, saved on ",
                format(file.mtime(tf_cache_file), "%Y-%m-%d"), "). Delete this file to obtain the collection again")
      }
    }

    if(is.null(universe)){
      if(TF.collection == "CollecTRI"){
        universe = decoupleR::get_collectri(organism = 'human', split_complexes = F)
      } else {
        universe = dplyr::filter(dorothea::dorothea_hs, .data$confidence %in% c("A", "B")) %>%
          dplyr::mutate(source = .data$tf) %>%
          dplyr::select(-tf)
      }
      if(!is.data.frame(universe) || nrow(universe) == 0 || !all(c("source", "target", "mor") %in% colnames(universe))){
        stop("The ", TF.collection, " collection could not be obtained (empty, or without the columns 'source', 'target' and 'mor')")
      }
      # Written to a temporary file first, so that an interrupted or simultaneous run never leaves or reads an incomplete cache
      utils::write.csv(universe, paste0(tf_cache_file, ".tmp"))
      file.rename(paste0(tf_cache_file, ".tmp"), tf_cache_file)
    }

    sample_acts <- decoupleR::decouple(mat     = RNA.counts,
                                       network = universe,
                                      .source = "source",
                                      .target = "target",
                                      statistics = if(method == "consensus") NULL else method,
                                      consensus_score = (method == "consensus"),
                                      minsize = min_targets_size
                                    ) %>%
                                      dplyr::filter(.data$statistic == method) %>%
                                      decoupleR::pivot_wider_profile(id_cols     = source,
                                                                     names_from  = condition,
                                                                     values_from = score) %>%
                                      as.matrix() %>%
                                      t()

  }

  sample_acts <- sample_acts[colnames(RNA.counts), , drop = FALSE]
  sample_acts <- if (scale) base::scale(sample_acts) else sample_acts

  result <- data.frame(sample_acts)
  colnames(result) <- make.names(colnames(result))

  if(return){
    utils::write.csv(result, paste0("Results/TF_matrix", if(!is.null(file.name)) paste0("_", file.name), ".csv"))
  }

  return(result)

}

#' Compute Weighted TF-coactivity Network Analysis (WTCNA)
#'
#' Construct a weighted signed or unsigned network using TF activity to cluster protein regulators
#' into modules that share similar activity patterns. Each TF module will have a sample-level score
#' represented by the eigengene of the module.
#'
#' @param TFs.matrix Matrix of TF activity (samples x TFs). When \code{batch = TRUE}, a list of such matrices,
#'   one per cohort.
#' @param batch Logical; if TRUE, performs consensus WGCNA (\code{WGCNA::blockwiseConsensusModules()}) across cohorts
#'   provided as a list of matrices. In this mode \code{clustering.method} is not used (WGCNA builds the consensus tree
#'   itself), modules are merged by WGCNA with \code{mergeCutHeight = 1 - corr_mod}, and module eigengenes are scaled
#'   within each cohort.
#' @param network.type Network type: "signed" (default) or "unsigned".
#' @param clustering.method Linkage method of the hierarchical clustering of the TFs (\code{stats::hclust()}), applied to the
#'   TOM dissimilarity (1 - topological overlap) between TFs (single-cohort mode only). Default is "average", as in WGCNA:
#'   two groups of TFs are joined according to the average dissimilarity between their TFs.
#' @param minMod Minimum number of TFs per module. Default is 15.
#' @param corr_mod Correlation threshold (0-1): modules whose eigengenes have a correlation (of type \code{cor_type})
#'   above it are merged. Default is 0.9.
#' @param cor_type Correlation used to pick the soft-threshold power, build the adjacency matrix, merge modules and
#'   compute the proportion of variance explained: "p" (Pearson) or "s" (Spearman). The module eigengene is always the
#'   first principal component of the module's TF activities. Spearman is only available when \code{batch = FALSE}, because WGCNA
#'   consensus modules only support Pearson correlation. Default is "p".
#' @param softPower Optional numeric value specifying the soft-thresholding power used to build the adjacency
#'   matrix (one value per cohort when \code{batch = TRUE}). If \code{NULL}, the lowest power whose scale-free fit
#'   \eqn{R^2} reaches 0.9 is chosen automatically (the power with the best fit if none reaches 0.9, with a message).
#' @param verbose Logical; whether to print the function messages. Default is FALSE.
#' @param file.name Optional character suffix used when writing WTCNA outputs (CSV files and, in single-cohort
#'   mode, the diagnostic plots \code{Soft_Threshold_<file.name>.pdf},
#'   \code{Gene_dendrogram_and_module_colors_<file.name>.pdf} and
#'   \code{Gene_dendrogram_and_module_colors_after_merging_<file.name>.pdf}).
#' @param return Logical, whether to save output plots, the module list (\code{Results/TFs_modules_<file.name>.csv})
#'   and the module matrix (\code{Results/TF_module_matrix_<file.name>.csv}, with the sample names as row names)
#'   to "Results/" (without \code{file.name}, \code{TFs_modules.csv} and \code{TF_module_matrix.csv}). Default is TRUE.
#'
#' @return A named list with:
#' \itemize{
#'   \item \code{TFs module matrix}: Scaled module eigengenes (samples x modules). In batch mode, samples are
#'     concatenated in cohort order (use the row names to match them to other tables).
#'   \item \code{TFs colors}: Vector of module colors assigned to each TF (\code{"grey"} for TFs not assigned to any module).
#'   \item \code{TFs per module}: List of TF names in each module.
#'   \item \code{Proportion of variance}: Variance explained per module (single-cohort mode only).
#'   \item \code{TFs_matrix}: The TF activity matrix used (a list of per-cohort matrices restricted to shared TFs in batch mode).
#'   \item \code{Network parameters}: A list with the soft-thresholding power(s) used (\code{softPower}, one per cohort in
#'     batch mode), \code{network.type} and \code{cor_type}, used by \code{identify_hub_TFs()}.
#' }
#'
#' @details
#' TFs that WGCNA does not assign to any module (color \code{"grey"}) are not a module: they are left out of the
#' module matrix and of \code{TFs per module}, with a message.
#'
#' In batch mode, the soft-thresholding powers are chosen on the TFs shared by all cohorts, the ones used to build
#' the network.
#'
#' In batch mode, \code{WGCNA::blockwiseConsensusModules()} can fail with
#' \code{"dim(X) must have a positive length"} when a consensus module has exactly one TF to reassign to
#' another module (a known WGCNA issue). The problem goes away when the modules change. If the powers were
#' chosen automatically (\code{softPower = NULL}), \code{compute.WTCNA()} retries moving the power of one cohort
#' one step down or up, and gives a warning with the powers that were finally used. If \code{softPower} was
#' given, or no neighbouring power works, it stops with a message reporting the powers: set \code{softPower}
#' to slightly different values or change \code{minMod} (a small change may not be enough).
#'
#' @references
#' Langfelder, P., & Horvath, S. (2008). WGCNA: an R package for weighted correlation network analysis.
#' BMC Bioinformatics, 9, 559. https://doi.org/10.1186/1471-2105-9-559
#'
#'
#' @examples
#'
#' data("tfs.tuto")
#' network <- compute.WTCNA(tfs.tuto, corr_mod = 0.9, return = FALSE)
compute.WTCNA <- function(TFs.matrix, batch = FALSE, network.type = "signed", clustering.method = "average",
                          minMod = 15, corr_mod = 0.9, cor_type = "p", verbose = F, file.name = NULL,
                          softPower = NULL, return = T) {

  if(verbose){
    cat("Creating weighted TF-coactivity network......................................................................\n\n")
  }

  if(!network.type %in% c("signed", "unsigned")){
    stop("network.type must be 'signed' or 'unsigned'")
  }
  if(!cor_type %in% c("p", "s")){
    stop("cor_type must be 'p' (Pearson) or 's' (Spearman)")
  }
  if(batch && (!is.list(TFs.matrix) || is.data.frame(TFs.matrix))){
    stop("If batch = TRUE, TFs.matrix must be a list of matrices, one per cohort")
  }
  if(!all(sapply(if(batch) TFs.matrix else list(TFs.matrix), function(x) all(sapply(as.data.frame(x), is.numeric))))){
    stop("TFs.matrix must contain only numeric TF activity values")
  }

  if(batch){
    if(cor_type != "p"){
      stop("cor_type = '", cor_type, "' is only available when batch = FALSE: WGCNA consensus modules only support Pearson correlation (cor_type = 'p')")
    }
    nCohorts <- length(TFs.matrix)
    if(verbose) cat("Running consensus WGCNA across", nCohorts, "cohorts...\n")

    # Step 1: intersect TFs across cohorts so all matrices share the same columns (powers are chosen on these TFs)
    common_tfs <- Reduce(intersect, lapply(TFs.matrix, colnames))
    if (length(common_tfs) == 0) stop("No common TFs found across cohorts")
    TFs.matrix <- lapply(TFs.matrix, function(x) x[, common_tfs, drop = FALSE])

    # Step 2: pick soft-threshold per cohort if not provided
    powers = c(1:10, 12:20)
    auto_power = is.null(softPower)
    if(auto_power){
      softPower <- numeric(nCohorts)
      for(i in 1:nCohorts){
        if(verbose) cat("Picking soft threshold for cohort", i, "...\n")
        utils::capture.output(sft <- WGCNA::pickSoftThreshold(TFs.matrix[[i]], powerVector = powers, verbose = 0, networkType = network.type)) # it prints its table otherwise
        target = 0.9
        fit = -sign(sft$fitIndices[,3]) * sft$fitIndices[,2] # signed scale-free fit R^2
        reached = which(fit >= target)
        softPower[i] <- if(length(reached) > 0) powers[reached[1]] else powers[which.max(fit)] # lowest power reaching the target
        if(length(reached) == 0) message("Cohort ", i, ": the scale-free fit never reaches ", target, " (best ", round(max(fit, na.rm = TRUE), 2), " at power ", softPower[i], "), this power is used")
        if(verbose) cat("Cohort", i, "softPower =", softPower[i], "\n")
      }
    } else {
      if(length(softPower) != nCohorts){
        stop("softPower must be a vector with one value per cohort")
      }
    }

    multiExpr <- lapply(TFs.matrix, function(x) list(data = x))

    # Step 3: run blockwiseConsensusModules
    # Known WGCNA issue: it fails ("dim(X) must have a positive length") when exactly one TF of a module is a candidate
    # for reassignment. NULL is returned in that case; any other error is raised as it is
    run_consensus <- function(power) tryCatch(WGCNA::blockwiseConsensusModules(
      multiExpr,
      power = power,                      # vector of powers per cohort
      networkType = network.type,
      minModuleSize = minMod,
      deepSplit = 2,
      pamRespectsDendro = FALSE,
      mergeCutHeight = 1 - corr_mod,       # merge modules whose eigengenes correlate above corr_mod, as in single-cohort mode
      corType = ifelse(cor_type == "p", "pearson", "spearman"),
      verbose = ifelse(verbose, 3, 0)
    ), error = function(e) if(grepl("dim(X) must have a positive length", conditionMessage(e), fixed = TRUE)) NULL else stop(e))

    consensusModules <- run_consensus(softPower)
    power_used <- softPower

    # The issue goes away when the modules change: with automatic powers, move the power of one cohort one step down or up
    if(is.null(consensusModules) && auto_power){
      for(step in c(-1, 1)){
        for(i in seq_len(nCohorts)){
          pos <- match(softPower[i], powers) + step
          if(is.null(consensusModules) && pos >= 1 && pos <= length(powers)){
            candidate <- softPower
            candidate[i] <- powers[pos]
            consensusModules <- run_consensus(candidate)
            if(!is.null(consensusModules)){
              power_used <- candidate
              warning("WGCNA::blockwiseConsensusModules() failed with the automatically chosen soft-thresholding powers (", paste(softPower, collapse = ", "),
                      ") because of a known WGCNA issue (a consensus module with exactly one TF to reassign): the powers ", paste(candidate, collapse = ", "),
                      " were used instead. Set 'softPower' to choose them yourself.", call. = FALSE)
            }
          }
        }
      }
    }

    if(is.null(consensusModules)){
      stop("WGCNA::blockwiseConsensusModules() failed because of a known WGCNA issue (a consensus module with exactly one TF ",
           "to reassign). It goes away when the modules change: set 'softPower' to other values than the ones used here (",
           paste(softPower, collapse = ", "), "; one step up or down is usually enough), or change 'minMod' (currently ", minMod,
           "; a small change may not be enough).", call. = FALSE)
    }

    dynamicColors <- consensusModules$colors

    # Step 4: Extract module eigengenes per cohort and scale within cohort (grey = TFs without module, not a module)
    MEs_list <- consensusModules$multiMEs
    scaled_MEs_list <- lapply(MEs_list, function(x) {
      ME <- x$data[, colnames(x$data) != "MEgrey", drop = FALSE]              # extract the data.frame
      ME_scaled <- scale(ME)    # scale MEs scores to account for variabilities across cohorts even though they are already scaled at TF level (MEs can be different in scale)
      colnames(ME_scaled) <- colnames(ME)
      return(ME_scaled)
    })

    # Step 5: concatenate scaled MEs across cohorts
    MEs <- do.call(rbind, scaled_MEs_list)
    if(ncol(MEs) == 0) stop("No TF module found: no TF was assigned to a consensus module. Try a smaller minMod")
    colnames(MEs) <- gsub("ME", "", colnames(MEs))

    # Step 6: Map TFs to modules (same as before)
    modtfs <- list()
    modules <- setdiff(unique(dynamicColors), "grey")
    tfs <- colnames(TFs.matrix[[1]])
    for(i in seq_along(modules)){
      modtfs[[i]] <- tfs[dynamicColors == modules[i]]
    }
    names(modtfs) <- modules

    output <- list(
      "TFs module matrix" = MEs,
      "TFs colors" = dynamicColors,
      "TFs per module" = modtfs,
      "TFs_matrix" = TFs.matrix,
      "Network parameters" = list(softPower = power_used, network.type = network.type, cor_type = cor_type)
    )


  } else {


    #####Choose parameter for scale-free network topology
    powers = c(c(1:10), seq(from = 12, to=20, by=1))
    utils::capture.output(sft <- WGCNA::pickSoftThreshold(TFs.matrix, powerVector = powers, verbose = 0, networkType = network.type,
                                                          corOptions = list(use = "p", method = ifelse(cor_type == "p", "pearson", "spearman")))) # it prints its table otherwise

    if(return){
      pdf(paste0("Results/Soft_Threshold", if(!is.null(file.name)) paste0("_", file.name), ".pdf"))
      plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
           xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
           main = paste("Scale independence"))
      graphics::text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=0.9,col="red");
      graphics::abline(h=0.90,col="red")
      dev.off()
    }

    # Automatic or user-defined soft-threshold selection
    if (is.null(softPower)) {
      target = 0.9 # Target SFT.R.sq value
      fit = -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2] # signed scale-free fit R^2
      reached = which(fit >= target)
      softPower = if(length(reached) > 0) powers[reached[1]] else powers[which.max(fit)] # lowest power reaching the target
      if(length(reached) == 0) message("The scale-free fit never reaches ", target, " (best ", round(max(fit, na.rm = TRUE), 2), " at power ", softPower, "): this power is used")

      if (verbose) {
        cat("Automatically choosing", softPower, "as soft-threshold......................................................................\n\n")
      }
    } else {
      if (verbose) {
        cat("Using user-defined soft-threshold =", softPower, "......................................................................\n\n")
      }
    }

    if(verbose){
      #####Co-expression matrix using nodes adjacency and topological overlapping nodes
      cat("Calculating nodes adjacency and topological overlapping nodes.................................................\n\n")
    }

    adjacency = WGCNA::adjacency(TFs.matrix, power =softPower, type=network.type, corFnc = "cor", corOptions = list(use = "p", method = ifelse(cor_type == "p", "pearson", "spearman")))
    TOM = WGCNA::TOMsimilarity(adjacency, TOMType = network.type, verbose = 0)
    dissTOM = 1-TOM

    #####Unsupervised hierarchical clustering using dissimilarity matrix
    geneTree = stats::hclust(stats::as.dist(dissTOM), method = clustering.method) # Tree built on the TOM dissimilarity itself, as in WGCNA
    dynamicMods = dynamicTreeCut::cutreeDynamic(dendro = geneTree, distM = dissTOM,
                                                deepSplit = 2, pamRespectsDendro = FALSE,
                                                minClusterSize = minMod, verbose = 0);
    dynamicColors = WGCNA::labels2colors(dynamicMods)

    #####Remove variables and clean garbage
    rm(adjacency, TOM, dissTOM)
    gc()

    if(return){
      pdf(paste0("Results/Gene_dendrogram_and_module_colors", if(!is.null(file.name)) paste0("_", file.name), ".pdf"))
      WGCNA::plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                                 dendroLabels = FALSE, hang = 0.03,
                                 addGuide = TRUE, guideHang = 0.05,
                                 main = "Gene dendrogram and module colors")
      dev.off()
    }

    #####Calculate eigenvectors from modules
    if(verbose){
      cat("Calculating eigenvectors from modules.................................................\n\n")
    }

    MEList = WGCNA::moduleEigengenes(TFs.matrix, colors = dynamicColors, scale = F) ## turn off scaling cause the output of compute.TFs.activity() is TFs = scale(NES)
    MEs = scale(MEList$eigengenes) # Scale MEs scores (not the same as above) and this is to have standard ranges and same input whether batch = TRUE or = FALSE (for this case)
    MEs = WGCNA::orderMEs(MEs)
    MEs = MEs[, colnames(MEs) != "MEgrey", drop = FALSE] # grey = TFs without module, not a module
    if(ncol(MEs) == 0) stop("No TF module found: no TF was assigned to a module. Try a smaller minMod")

    if(verbose){
      print(paste0("Merging modules significantly correlated with ", corr_mod, "........"))
    }

    merge = mergeModules(MEs, dynamicColors, corr_mod, method = ifelse(cor_type == "p", "pearson", "spearman"))
    MEs = merge[[1]]
    dynamicColors = merge[[2]]
    modtfs = list()
    for(i in 1:ncol(MEs)){
      tfs = colnames(TFs.matrix)
      modules = c(substring(names(MEs), 3))
      inModule = is.finite(match(dynamicColors,modules[i]))
      modtfs[[i]] = tfs[inModule]
    }
    names(modtfs) = modules

    if(return){
      pdf(paste0("Results/Gene_dendrogram_and_module_colors_after_merging", if(!is.null(file.name)) paste0("_", file.name), ".pdf"))
      WGCNA::plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                                 dendroLabels = FALSE, hang = 0.03,
                                 addGuide = TRUE, guideHang = 0.05,
                                 main = "Gene dendrogram and module colors")
      dev.off()
    }

    TFspropVar = WGCNA::propVarExplained(TFs.matrix[, dynamicColors != "grey", drop = FALSE], dynamicColors[dynamicColors != "grey"], MEs, corFnc = "cor", corOptions = if(cor_type == "p") "use = 'p'" else "use = 'p', method = 'spearman'")

    colnames(MEs) <- gsub("ME", "", colnames(MEs))

    output = list(MEs, dynamicColors, modtfs, TFspropVar, TFs.matrix, list(softPower = softPower, network.type = network.type, cor_type = cor_type))
    names(output) = c("TFs module matrix", "TFs colors", "TFs per module", "Proportion of variance", "TFs_matrix", "Network parameters")

  }

  unassigned = colnames(if(batch) TFs.matrix[[1]] else TFs.matrix)[dynamicColors == "grey"]
  if(length(unassigned) > 0){
    message(length(unassigned), " TFs not assigned to any module (grey) are left out of the modules: ", paste(utils::head(unassigned, 10), collapse = ", "),
            if(length(unassigned) > 10) ", ...")
  }

  ## Retrieve TFs modules in .csv
  contador = 1
  tfs_modules = data.frame(matrix(nrow = length(modtfs), ncol = 2))
  colnames(tfs_modules) = c("TFs module", "Composition")
  for (i in 1:length(modtfs)) {
    tfs_modules[contador,1] = names(modtfs)[i]
    tfs_modules[contador,2] = paste(modtfs[[i]], collapse = ",")
    contador = contador + 1
  }

  if(return){
    file_suffix = if(!is.null(file.name)) paste0("_", file.name) else ""
    utils::write.csv(tfs_modules, paste0("Results/TFs_modules", file_suffix, ".csv"), row.names = F)
    utils::write.csv(MEs, file = paste0("Results/TF_module_matrix", file_suffix, ".csv")) # Sample names as row names (in batch mode samples are stacked cohort by cohort)
  }

  return(output)
}

#' Identify hub TFs
#'
#' Identifies hub TFs per module using values of module membership and degree.
#' TFs with high module membership (Pearson correlation of their activity with the module score > \code{MM_thresh})
#' and a degree at or above the \code{degree_thresh} quantile of their module are considered hub TFs.
#' The degree of a TF is its intramodular connectivity, as in WGCNA: the sum of its network adjacencies with the
#' other TFs of the same module. The adjacencies are computed as in \code{compute.WTCNA()}, with the same
#' soft-thresholding power, network type and correlation (for a \code{"signed"} network, \eqn{((1 + r)/2)^{power}}, so a
#' negative correlation adds almost nothing instead of being subtracted).
#'
#' The TF activities are taken from \code{TF.network} (element \code{TFs_matrix}) and matched to the module scores
#' by sample name. For a multi-cohort network (\code{compute.WTCNA(batch = TRUE)}), the per-cohort TF activities
#' are stacked to compute the membership, and the degree is the mean of the per-cohort degrees (each with the power
#' of its cohort).
#'
#' @param TF.network TF network obtained from \code{compute.WTCNA()}.
#' @param MM_thresh Threshold for module membership. Default is 0.8.
#' @param degree_thresh Quantile threshold for degree within each module. Default is 0.9 (top 10%).
#'
#' @return A list with two elements:
#' \describe{
#'   \item{hubGenes}{Named list of hub TFs per module. Modules without hub TFs are left out (a message lists them).}
#'   \item{detailedData}{Data frame with one row per hub TF (TF names as row names) and the columns \code{Module},
#'     \code{ModuleMembership} and \code{Degree}.}
#' }
#'
#' @details
#' Networks built before the soft-thresholding power was stored in \code{TF.network} (element \code{Network parameters})
#' have their power chosen again with the rule of \code{compute.WTCNA()} (with a message), assuming a signed network
#' and Pearson correlation.
#'
#' @export
#'
#' @examples
#'
#' data("network.tuto")
#'
#' hub_tfs <- identify_hub_TFs(network.tuto, MM_thresh = 0.8, degree_thresh = 0.9)
#'
identify_hub_TFs <- function(TF.network, MM_thresh = 0.8, degree_thresh = 0.9) {

  if(!is.list(TF.network) || is.data.frame(TF.network) || is.null(TF.network[["TFs module matrix"]]) || is.null(TF.network[["TFs_matrix"]])){
    stop("TF.network must be the output of compute.WTCNA(): the TF activities are taken from it (the datExpr argument was removed)")
  }
  moduleEigengenes_df = TF.network[["TFs module matrix"]]
  moduleColors = TF.network[["TFs colors"]]
  modules = intersect(unique(moduleColors), colnames(moduleEigengenes_df)) # TFs without module (grey) are not a module

  # TF activities: one matrix (samples x TFs), or one per cohort for a multi-cohort network
  tfs_list = TF.network[["TFs_matrix"]]
  if(is.data.frame(tfs_list) || is.matrix(tfs_list)) tfs_list = list(tfs_list)
  tfs_list = lapply(tfs_list, as.matrix)
  tf_names = colnames(tfs_list[[1]])
  if(length(moduleColors) != length(tf_names)){
    stop("TF.network is inconsistent: ", length(moduleColors), " module colors for ", length(tf_names), " TFs")
  }

  # Samples matched to the module scores by name (cohort matrices stacked in batch mode)
  tfs_all = do.call(rbind, unname(tfs_list))
  idx = match(make.names(rownames(moduleEigengenes_df)), make.names(rownames(tfs_all)))
  if(is.null(rownames(moduleEigengenes_df)) || is.null(rownames(tfs_all)) || anyNA(idx)){
    stop("The samples of the module matrix and of the TF activities in TF.network do not match (sample names are needed)")
  }
  tfs_all = tfs_all[idx, , drop = FALSE]

  # Network parameters of compute.WTCNA() (chosen again for networks built before they were stored)
  params = TF.network[["Network parameters"]]
  if(is.null(params)){
    powers = c(1:10, 12:20)
    softPower = sapply(tfs_list, function(x){
      utils::capture.output(sft <- WGCNA::pickSoftThreshold(x, powerVector = powers, verbose = 0, networkType = "signed"))
      fit = -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
      reached = which(fit >= 0.9)
      if(length(reached) > 0) powers[reached[1]] else powers[which.max(fit)]
    })
    params = list(softPower = softPower, network.type = "signed", cor_type = "p")
    message("TF.network does not store its soft-thresholding power (built with an older version of compute.WTCNA()): power ",
            paste(softPower, collapse = ", "), " chosen again, signed network and Pearson correlation assumed")
  }
  cor_method = ifelse(params$cor_type == "p", "pearson", "spearman")

  # Identify hub genes
  hubGenesList <- list()
  hubGenesData <- data.frame()  # Initialize an empty data frame

  for (module in modules) {
    inModule <- which(moduleColors == module)

    # Module membership: correlation of each TF with the module score
    memberships <- stats::cor(tfs_all[, inModule, drop = FALSE], moduleEigengenes_df[, module])

    # Degree: intramodular connectivity as in WGCNA (sum of the adjacencies with the other TFs of the module), mean over cohorts
    degrees_per_cohort <- sapply(seq_along(tfs_list), function(i){
      adjacency <- WGCNA::adjacency(tfs_list[[i]][, inModule, drop = FALSE], power = params$softPower[i], type = params$network.type,
                                    corFnc = "cor", corOptions = list(use = "p", method = cor_method))
      diag(adjacency) <- 0
      rowSums(adjacency)
    })
    degrees <- rowMeans(matrix(degrees_per_cohort, nrow = length(inModule)))
    names(degrees) <- tf_names[inModule]

    # Calculate the cutoff for the top 10% by degree
    degreeCutoff <- stats::quantile(degrees, degree_thresh)

    # Identify hub genes (top 10% by degree)
    hubGenes <- names(degrees[degrees >= degreeCutoff])
    if (length(hubGenes) == 0) next

    # Identify hub genes (>0.8 MM)
    highMMGenes <- which(memberships > MM_thresh)
    if (length(highMMGenes) == 0) next

    # Filter hub genes (degree) to only include those with high MM
    finalHubGenes <- intersect(hubGenes, rownames(memberships)[highMMGenes])
    if (length(finalHubGenes) == 0) next

    hubGenesList[[module]] = finalHubGenes

    # Create a dataframe with detailed information for hub genes only
    moduleData <- data.frame(
      Module = module,
      ModuleMembership = unname(memberships[finalHubGenes,]),
      Degree = degrees[finalHubGenes],
      stringsAsFactors = FALSE)

    # Append to the combined dataframe
    hubGenesData <- rbind(hubGenesData, moduleData)
  }

  no_hubs <- setdiff(modules, names(hubGenesList))
  if (length(no_hubs) > 0) {
    message("No hub TFs in ", length(no_hubs), " of ", length(modules), " modules (", paste(no_hubs, collapse = ", "), "): they are not in the output")
  }

  return(list(hubGenes = hubGenesList, detailedData = hubGenesData))
}

#' Identify cell groups
#'
#' Builds, for each TF module, a dendrogram of the deconvolution features significantly correlated with it.
#' The features are clustered on the value of that correlation: the distance between two features is the absolute
#' difference between their correlations with the module (e.g. 0.62 and 0.58 are at distance 0.04). Two features are
#' therefore similar when they are related to the module to the same degree (similar biological activity), not when
#' they vary together across samples. The dendrograms are cut into cell groups by \code{cell.groups.computation()}.
#'
#' @param features The list returned by \code{compute.modules.relationship()}. Its first two elements are used
#'   (accessed by position):
#'   \itemize{
#'     \item \code{[[1]]} (\code{correlation}): A matrix of correlations between TF modules (rows, named) and deconvolution features (columns).
#'     \item \code{[[2]]} (\code{significant}): A list, named by TF module, of the features significantly correlated with each module
#'       (p-value <= the \code{pval} given to \code{compute.modules.relationship()}).
#'   }
#'   The correlations of each module are taken from the row with its name.
#' @param clustering.method Linkage method of the hierarchical clustering (\code{stats::hclust()}) of the features
#'   on the distance described above. Default: "ward.D2".
#' @param file_name Optional character suffix for the names of the saved plots.
#' @param return Logical; whether to save one dendrogram plot per module to
#'   \code{Results/Dendrogram_cell_types_<module>_<file_name>.pdf} (\code{Results/Dendrogram_cell_types_<module>.pdf}
#'   without \code{file_name}). Each plot shows the correlation of every feature with the module next to its name, and
#'   its size follows the number of features. Default: TRUE.
#'
#' @return A named list of dendrograms (\code{hclust} objects, labelled with the feature names), one per TF module. Each
#'   one also stores the correlation of its features with the module (element \code{correlation}, used by
#'   \code{calculate_dendrogram_cuts()}). Modules with fewer than two significant features are discarded; \code{NULL} is
#'   returned if no module remains.
#'
#' @examples
#' data("network.tuto")
#' data("deconv_subgroups.tuto")
#' rel <- compute.modules.relationship(network.tuto[[1]], deconv_subgroups.tuto[[1]], plot = FALSE)
#' dendrograms <- identify.cell.groups(rel, return = FALSE)
#' names(dendrograms)
#'
identify.cell.groups = function(features, clustering.method = "ward.D2", file_name = NULL, return = T){

  if(!is.list(features) || length(features) < 2 || !is.matrix(features[[1]]) || !is.numeric(features[[1]]) || !is.list(features[[2]]) ||
     is.null(rownames(features[[1]])) || is.null(names(features[[2]]))){
    stop("features must be the list returned by compute.modules.relationship(): a correlation matrix with TF modules as row names, ",
         "and a list of significant features named by TF module")
  }
  moduleTraitCor = features[[1]]
  features_vec = features[[2]]
  missing_modules = setdiff(names(features_vec), rownames(moduleTraitCor))
  if(length(missing_modules) > 0){
    stop("No correlations for the TF modules: ", paste(missing_modules, collapse = ", "))
  }

  #Discard modules with fewer than two significant features
  features_vec = features_vec[lengths(features_vec) >= 2]

  if(length(features_vec) == 0){
    return(NULL)
  }
  lis.dendrogram = list()

  for (module in names(features_vec)){
    TFmoduleTraitcor = moduleTraitCor[module, colnames(moduleTraitCor) %in% features_vec[[module]], drop = F] # Correlations of this module (by name)
    ###Dendrogram by Module

    d <- stats::dist(t(TFmoduleTraitcor))
    dendrogram <- stats::hclust(d, method = clustering.method)
    dendrogram$correlation <- TFmoduleTraitcor[1, dendrogram$labels] # Kept for calculate_dendrogram_cuts(): each cell group must be one range of these values
    if(return){
      # Labels with the correlation of each feature; size from the number of features and the length of their labels
      dendrogram_plot <- dendrogram
      dendrogram_plot$labels <- paste0(dendrogram$labels, " (r = ", sprintf("%.2f", TFmoduleTraitcor[1, dendrogram$labels]), ")")
      label_lines <- max(nchar(dendrogram_plot$labels)) * 0.62
      pdf(paste0("Results/Dendrogram_cell_types_", module, if(!is.null(file_name)) paste0("_", file_name), ".pdf"),
          width = 4 + label_lines * 0.2 * 0.8, height = max(4, 1.5 + 0.16 * length(dendrogram$labels)))
      par(mar = c(5, 1, 4, label_lines), cex = 0.8) #bottom, left, top, right
      plot(as.dendrogram(dendrogram_plot), horiz = T, xlab = paste0("Merge height (", clustering.method, ") on differences in correlation"),
           main = paste0("TF module ", module, ": ", length(dendrogram$labels), " significantly correlated features"),
           adj = 0, cex.main = 1) # title starts at the left edge of the tree and runs over the labels, so it is not cut
      dev.off()
    }
    lis.dendrogram[[module]] = dendrogram
  }

  return(lis.dendrogram)

}

#' Compute a cell-type composition matrix from deconvolution subgroups
#'
#' Builds a binary presence matrix indicating which original cell types are present in each cell group. Each feature
#' of a cell group is opened down to its original deconvolution features (a subgroup is replaced by its members), and
#' the cell type of each feature is read from its name with the cell type vocabulary of multideconv
#' (\code{multideconv::get_cell_type_nomenclature()}). Cell types are therefore named as in that vocabulary: for example
#' the LM22 features \code{Plasma.cells} and \code{T.cells.CD4.memory.activated} are counted as \code{Plasma} and
#' \code{CD4.memory.activated}, as in the names of the deconvolution subgroups.
#'
#' @param deconvolution.subgroupped Deconvolution subgroups as returned by
#'   \code{multideconv::compute.deconvolution.analysis()} (the elements \code{Deconvolution matrix} and
#'   \code{Deconvolution subgroups composition} are used).
#' @param cell.groups The cell groups returned by \code{construct_cell_groups()}: its element \code{Composition} (the
#'   features of each cell group) is used. An unnamed list whose second element holds the composition is also accepted.
#' @param cells_extra Optional character vector of additional cell type names, added to the vocabulary used to
#'   recognise cell types in the feature names.
#'
#' @return A binary matrix (data.frame) where rows are cell groups and columns are cell types (1 = present, 0 = absent).
#'   Each cell group appears twice, as \code{<group>_pos} and \code{<group>_neg}, to match the feature names used by
#'   \code{compute.latent_factors()}; both rows have the same composition (the direction of the relation with the TF
#'   module is not represented). Cell types not present in any cell group are dropped, and cell groups whose cell types
#'   are not recognised (their row is all zeros) are listed, each with a message.
#'
#' @examples
#' data("network.tuto")
#' data("deconv_subgroups.tuto")
#' cell_groups <- construct_cell_groups(network.tuto, deconv_subgroups.tuto)
#' composition <- compute.composition.matrix(deconv_subgroups.tuto, cell_groups)
#' composition[1:4, 1:5]
#'
compute.composition.matrix = function(deconvolution.subgroupped, cell.groups, cells_extra = NULL){

  if(!is.list(deconvolution.subgroupped) || is.null(deconvolution.subgroupped[["Deconvolution matrix"]]) ||
     is.null(deconvolution.subgroupped[["Deconvolution subgroups composition"]])){
    stop("deconvolution.subgroupped must be the output of multideconv::compute.deconvolution.analysis() (elements 'Deconvolution matrix' and 'Deconvolution subgroups composition')")
  }
  # Features of each cell group (element Composition of construct_cell_groups())
  composition = if(is.list(cell.groups) && !is.null(cell.groups[["Composition"]])) cell.groups[["Composition"]] else if(is.list(cell.groups) && is.null(names(cell.groups)) && length(cell.groups) >= 2) cell.groups[[2]] else NULL
  if(!is.list(composition) || length(composition) == 0 || is.null(names(composition)) || !all(sapply(composition, is.character))){
    stop("cell.groups must be the output of construct_cell_groups(), with at least one cell group (element 'Composition': the features of each cell group)")
  }

  # Initialize composition matrix
  cells_types = extract_cells(colnames(deconvolution.subgroupped[["Deconvolution matrix"]]), cells_extra = cells_extra)
  presence_matrix <- data.frame(matrix(data = 0, nrow = length(composition), ncol = length(cells_types)))
  colnames(presence_matrix) <- cells_types

  #Members of every deconvolution subgroup (subgroup name -> features)
  deconv_subgroups = list()
  contador = 1
  subgroups = deconvolution.subgroupped[["Deconvolution subgroups composition"]]
  for (i in seq_along(subgroups)) {
    if(length(subgroups[[i]])!=0){ #Whether a specific cell type does not contains subgroups
      for (j in 1:length(subgroups[[i]])) {
        deconv_subgroups[[contador]] = subgroups[[i]][[j]]
        names(deconv_subgroups)[contador] = names(subgroups[[i]])[j]
        contador = contador + 1
      }
    }
  }

  # Cell types of each cell group
  for (j in seq_along(composition)) {
    features <- composition[[j]] #Features of cell group j
    for (cell_feature in features) {
      cells <- get_all_cells(cell_feature, deconv_subgroups) #Use recursive function to extract all nested cell features from different subgroup levels
      cells_types = extract_cells(cells, cells_extra = cells_extra)
      presence_matrix[j, cells_types] <- 1 #Set 1 if cell feature is present in subgroup
    }
  }

  rownames(presence_matrix) = names(composition)

  not_recognised = rownames(presence_matrix)[rowSums(presence_matrix) == 0]
  if(length(not_recognised) > 0){
    message("Cell groups whose cell types are not recognised in the feature names (no cell type in the composition matrix; see cells_extra): ",
            paste(not_recognised, collapse = ", "))
  }

  idy = which(colSums(presence_matrix)==0)

  if(length(idy)!=0){
    message("Cell types not present in any cell group (left out of the composition matrix): ", paste0(names(idy), collapse = ", "))
    presence_matrix = presence_matrix[, -idy, drop = FALSE]
  }

  # One row per latent factor feature: <group>_pos and <group>_neg (same composition)
  presence_matrix = presence_matrix[rep(1:nrow(presence_matrix), each = 2), , drop = FALSE]

  rownames(presence_matrix) = paste0(rep(names(composition), each = 2),
                                     c("_pos","_neg"))

  return(presence_matrix)

}

#' Construct cell groups based on TF networks and deconvolution
#'
#' A cell group is a set of deconvolution features with a similar biological activity: features related to the
#' same TF module to the same degree (similar correlation with the module score). The features of a cell group
#' do not need to vary together across samples. The cell groups are built in three steps:
#' \enumerate{
#'   \item each TF module is correlated with every deconvolution feature (\code{compute.modules.relationship()});
#'   \item for each module, the significantly correlated features are clustered on their correlation with the
#'     module (\code{identify.cell.groups()});
#'   \item each dendrogram is cut into cell groups, and each cell group gets one score per sample from a CCA with the
#'     TFs of its module (\code{cell.groups.computation()}).
#' }
#' Samples are matched by name: the module scores and TF activities of \code{network} are put in the sample order of
#' the deconvolution matrix, which may contain only part of the samples of the network.
#'
#' @param network A TF module network as returned by \code{compute.WTCNA()} (its module scores, TF activities and
#'   TFs per module are used).
#' @param dt Deconvolution subgroups as returned by \code{multideconv::compute.deconvolution.analysis()}; its element
#'   \code{Deconvolution matrix} (samples x features) is used. Every sample must be in \code{network}.
#' @param batch Optional vector with the cohort of each sample: in the row order of the deconvolution matrix, or named by
#'   sample. Required when \code{network} was built with \code{compute.WTCNA(batch = TRUE)}. With two or more cohorts,
#'   correlations and scores are corrected for cohort.
#' @param pval Numeric. P-value threshold used to select the deconvolution features significantly
#'   correlated with each TF module. Default: 0.05.
#' @param clustering.method Linkage method of the hierarchical clustering of the deconvolution features of each module
#'   on their correlation with it (see \code{identify.cell.groups()}). Default: "ward.D2".
#' @param dendrogram_file Optional character. Suffix of the output file names (see \code{return_dendrogram}
#'   and \code{return}).
#' @param return_dendrogram Logical. If TRUE, saves a PDF of the colored cell-group dendrograms to
#'   \code{Results/Dendrogram_color_clusters_<dendrogram_file>.pdf} (only when \code{dendrogram_file} is set). Default FALSE.
#' @param return Logical. If TRUE, saves the composition and the scores of the cell groups to
#'   \code{Results/Cell.groups.composition_<dendrogram_file>.csv} and
#'   \code{Results/Cell.groups.scores_<dendrogram_file>.csv}. Default FALSE.
#' @param cells_extra Optional character vector with the cell names that are not in the nomenclature of multideconv
#'   (the same as in \code{multideconv::compute.deconvolution.analysis()}), used to list the cell types of each
#'   cell group in \code{Module_relation}.
#' @param min_cluster_size Integer. Minimum number of features of a cell group (minimum cluster size of the
#'   dynamic tree cut of each dendrogram). Default is 3.
#' @param deep_split Logical. If TRUE (default), the dynamic tree cut splits the clusters more finely.
#'
#' @return A named list of 4 elements:
#' \describe{
#'   \item{Cell_groups}{A data frame with the projected cell group scores (samples x groups). A positive score
#'     means that the sample is on the active side of the TF module of that cell group, a negative score that
#'     it is on the inactive side.}
#'   \item{Composition}{A named list where each element is a character vector of the deconvolution features in each group.}
#'   \item{Weights}{A named list of CCA projection parameters (\code{xcoef}, \code{train_means}, \code{train_sds}) for each group.
#'     With \code{batch}, \code{train_means} is \code{NULL}: test data are centred on their own means.}
#'   \item{Module_relation}{A data frame with one row per cell group: \code{cell_group}, \code{module} (its TF
#'     module), \code{relation} (\code{"positive"} or \code{"negative"}), \code{correlation} (mean correlation of
#'     its features with the module score), \code{n_features} and \code{cell_types} (the cell types of its
#'     features, the most frequent first). It translates a score into cell abundance: for a
#'     positively related group, a positive score means more of its features; for a negatively related group, a
#'     positive score means fewer of them.}
#' }
#' Cell groups are named \code{Dendrogram_<module>_group_<n>}, where \code{n} is the cluster number shown in the
#' dendrogram PDF of that module. Features that the tree cut does not assign to any cluster are not part of any
#' cell group (a message reports how many).
#' @examples
#' data("network.tuto")
#' data("deconv_subgroups.tuto")
#' cell_groups <- construct_cell_groups(network.tuto, deconv_subgroups.tuto)
#' head(cell_groups$Module_relation)
#'
#' @export
construct_cell_groups = function(network, dt, batch = NULL, pval = 0.05, clustering.method = "ward.D2", dendrogram_file = NULL, return_dendrogram = FALSE, return = FALSE, cells_extra = NULL, min_cluster_size = 3, deep_split = TRUE){

  if(!is.list(network) || is.data.frame(network) || is.null(network[["TFs module matrix"]]) || is.null(network[["TFs_matrix"]]) || is.null(network[["TFs per module"]])){
    stop("network must be the output of compute.WTCNA()")
  }
  if(!is.list(dt) || is.data.frame(dt) || is.null(dt[["Deconvolution matrix"]])){
    stop("dt must be the output of multideconv::compute.deconvolution.analysis() (element 'Deconvolution matrix')")
  }
  deconv = dt[["Deconvolution matrix"]]
  samples = make.names(rownames(deconv))

  # Multi-cohort network (compute.WTCNA(batch = TRUE)): TF activities are one matrix per cohort
  multi_cohort = is.list(network[["TFs_matrix"]]) && !is.data.frame(network[["TFs_matrix"]])
  if(multi_cohort && is.null(batch)){
    stop("network was built with compute.WTCNA(batch = TRUE): give the cohort of each sample in batch")
  }

  # Cohort of each sample, in the row order of the deconvolution matrix
  if(!is.null(batch)){
    if(!is.null(names(batch))){
      idx_batch = match(samples, make.names(names(batch)))
      if(anyNA(idx_batch)) stop("batch is named but has no cohort for the samples: ", paste(utils::head(rownames(deconv)[is.na(idx_batch)], 5), collapse = ", "))
      batch = unname(batch[idx_batch])
    } else if(length(batch) != nrow(deconv)){
      stop("batch must have one value per sample of the deconvolution matrix (", nrow(deconv), "), in the order of its rows, or be named by sample")
    }
  }

  # Module scores and TF activities in the sample order of the deconvolution matrix (matched by name)
  tfs = if(multi_cohort) do.call(rbind, unname(lapply(network[["TFs_matrix"]], as.matrix))) else network[["TFs_matrix"]] # unname() keeps the sample names
  idx_modules = match(samples, make.names(rownames(network[["TFs module matrix"]])))
  idx_tfs = match(samples, make.names(rownames(tfs)))
  missing_samples = rownames(deconv)[is.na(idx_modules) | is.na(idx_tfs)]
  if(length(missing_samples) > 0){
    stop(length(missing_samples), " samples of the deconvolution matrix are not in the TF network: ", paste(utils::head(missing_samples, 5), collapse = ", "),
         if(length(missing_samples) > 5) ", ...")
  }
  network[["TFs module matrix"]] = network[["TFs module matrix"]][idx_modules, , drop = F]
  network[["TFs_matrix"]] = tfs[idx_tfs, , drop = F]

  corr_modules = compute.modules.relationship(network[["TFs module matrix"]], deconv, batch = batch, plot = F, pval = pval)
  cell_dendrograms = identify.cell.groups(corr_modules, clustering.method = clustering.method, return = F)
  if(is.null(cell_dendrograms)){
    stop("No cell groups found: no TF module is significantly correlated (p-value <= ", pval, ") with at least two deconvolution features")
  }
  cell.groups = cell.groups.computation(deconv, cell_dendrograms, network, batch = batch, return = return, dendrogram_file = dendrogram_file, return_dendrogram = return_dendrogram, cells_extra = cells_extra, min_cluster_size = min_cluster_size, deep_split = deep_split)

  names(cell.groups) = c("Cell_groups", "Composition", "Weights", "Module_relation")
  return(cell.groups)
}

#' Compute composite scores on test set based on previous cell groups
#'
#' Scores the cell groups learned on a training set in an independent (test) set, on the same scale as the training
#' scores. Three steps:
#' \enumerate{
#'   \item the training deconvolution subgroups are rebuilt in the test data with
#'     \code{multideconv::replicate_deconvolution_subgroups()} (each subgroup is the median of its member features);
#'   \item for each cell group, its features are scaled with the training means and standard deviations;
#'   \item they are projected with the training CCA weights of the cell group, giving one score per test sample.
#' }
#' Scoring the training samples themselves reproduces their training scores. Each test sample is scored on its own,
#' so a single sample can be scored (except for cell groups built with \code{batch}, see Details).
#'
#' @param deconv_res Deconvolution subgroups of the training set, as returned by
#'   \code{multideconv::compute.deconvolution.analysis()} (the elements \code{Deconvolution matrix} and
#'   \code{Deconvolution subgroups composition} are used).
#' @param cell_groups The cell groups of the training set, as returned by \code{construct_cell_groups()}: its elements
#'   \code{Composition} (features of each cell group) and \code{Weights} (CCA projection parameters \code{xcoef},
#'   \code{train_means}, \code{train_sds} of each cell group) are used. Older unnamed lists, with these as second and
#'   third elements, are also accepted.
#' @param features Optional character vector with the names of the cell groups to score. If \code{NULL} (default), all
#'   cell groups are scored. Names that are not cell groups are ignored with a message; an error is raised if none
#'   of them is a cell group.
#' @param deconvolution_test A data frame or matrix with the deconvolution features of the test set before
#'   subgrouping (the same kind of matrix as given to \code{multideconv::compute.deconvolution.analysis()} for the
#'   training set, with the same column names), samples as rows.
#'
#' @return A data frame with one row per test sample (sample names as row names) and one column per scored cell group,
#'   on the scale of the training scores. A cell group is scored only if all its features are in the test set without
#'   missing values: a score computed from part of them can be very different from the training score (even
#'   reversed), so the other cell groups are skipped, with a warning listing them. An error is raised if no cell group
#'   can be scored.
#'
#' @details
#' If the cell groups were built with batch correction (\code{batch} in \code{construct_cell_groups()}), the test
#' samples are centred on their own means (as each training cohort was), so \code{deconvolution_test} should contain a
#' single cohort with several samples: the scores then describe each sample relative to the other test samples.
#' A warning is given when such cell groups are scored with fewer than 10 test samples.
#'
#' @examples
#' data("network.tuto")
#' data("deconv_subgroups.tuto")
#' data("deconv.tuto")
#' cell_groups <- construct_cell_groups(network.tuto, deconv_subgroups.tuto)
#' test_scores <- compute.test.set(deconv_subgroups.tuto, cell_groups, deconvolution_test = deconv.tuto[1:20, ])
#' test_scores[1:3, 1:3]
#'
compute.test.set = function(deconv_res, cell_groups, features = NULL, deconvolution_test){

  # Features and CCA parameters of each cell group (by name; older unnamed lists by position)
  unnamed = is.list(cell_groups) && is.null(names(cell_groups)) && length(cell_groups) >= 3
  composition = if(is.list(cell_groups) && !is.null(cell_groups[["Composition"]])) cell_groups[["Composition"]] else if(unnamed) cell_groups[[2]] else NULL
  weights = if(is.list(cell_groups) && !is.null(cell_groups[["Weights"]])) cell_groups[["Weights"]] else if(unnamed) cell_groups[[3]] else NULL
  if(!is.list(composition) || !is.list(weights) || length(composition) == 0 || length(weights) != length(composition) || is.null(names(composition))){
    stop("cell_groups must be the output of construct_cell_groups() (elements 'Composition' and 'Weights')")
  }
  if(is.null(features)) features = names(composition)
  unknown = setdiff(features, names(composition))
  if(length(unknown) > 0 && length(unknown) < length(features)){
    message("Not cell groups of cell_groups, ignored: ", paste(utils::head(unknown, 10), collapse = ", "), if(length(unknown) > 10) paste0(", ... (", length(unknown), " in total)"))
  }

  ################################################################################Replicate cell subgroups (training features, in the same order)
  deconvolution_test = replicate_deconvolution_subgroups(deconv_res, deconvolution_test)

  # Compute composite scores
  idx = which(names(composition) %in% features)
  if(length(idx) == 0){
    stop("None of the features match the cell group names")
  }
  cell_dendrogram = c()
  names = c()
  skipped = c()
  for (i in 1:length(idx)) {
    group = names(composition)[idx[i]]
    cells = composition[[idx[i]]]
    pca_cells = deconvolution_test[,colnames(deconvolution_test) %in% cells, drop = F]
    pca_cells <- pca_cells[, colSums(is.na(pca_cells)) == 0, drop = FALSE] # Drop features with missing values (the training means and SDs are used to scale, so constant features can be scored)
    loadings_cells = if(!is.null(names(weights)) && group %in% names(weights)) weights[[group]] else weights[[idx[i]]]
    if(!all(rownames(loadings_cells$xcoef) %in% colnames(pca_cells))){ # A feature is missing: a partial score would not be comparable with the training score
      skipped = c(skipped, group)
      next
    }
    cell_dendrogram = cbind(cell_dendrogram, compute.test.score(pca_cells, loadings_cells))
    names = c(names, group)
  }

  if(is.null(cell_dendrogram)){
    stop("No cell group could be scored: every cell group has features missing, or with missing values, in the test set (check that the ",
         "column names of deconvolution_test are the deconvolution features used for training)")
  }
  if(length(skipped) > 0){
    warning("Cell groups skipped because some of their features are missing or have missing values in the test set (their score would not ",
            "be comparable with the training scores): ",
            paste(utils::head(skipped, 10), collapse = ", "), if(length(skipped) > 10) paste0(", ... (", length(skipped), " in total)"), call. = FALSE)
  }

  colnames(cell_dendrogram) = names

  # Cell groups built with batch: the test samples are their own reference (centred on their own means)
  batch_trained = any(sapply(names, function(group) is.null((if(!is.null(names(weights)) && group %in% names(weights)) weights[[group]] else weights[[match(group, names(composition))]])$train_means)))
  if(batch_trained && nrow(deconvolution_test) < 10){
    warning("The cell groups were built with batch correction, so the test samples are centred on their own means and their scores depend ",
            "on the other test samples. With fewer than 10 test samples (here ", nrow(deconvolution_test), ") the scores are unreliable",
            if(nrow(deconvolution_test) == 1) " (with one sample, every score is 0)", ": give a whole cohort as test set.", call. = FALSE)
  }

  return(data.frame(cell_dendrogram, check.names = FALSE))

}

#' Extract cell types from deconvolution feature names
#'
#' Reads the cell type of each deconvolution feature from its name, with the cell type vocabulary of multideconv
#' (\code{multideconv::get_cell_type_nomenclature()}, plus \code{cells_extra}). The first vocabulary name found in the
#' feature name is kept, and a \code{_Subgroup.<n>} suffix (or \code{.Iteration.<n>} of multideconv < 0.2.0) is removed.
#' Vocabulary names are matched literally (a "." only matches a ".").
#'
#' @param groups A character vector of deconvolution feature names: raw features (\code{<method>_<signature>_<cell type>},
#'   or \code{<method>_<cell type>}) or subgroups (\code{<cell type>_Subgroup.<n>}).
#' @param cells_extra Optional character vector of additional cell type names, added to the vocabulary.
#'
#' @return A character vector of the distinct cell types found (\code{character(0)} if none). Feature names without
#'   any vocabulary cell type are ignored.
#'
#' @examples
#' extract_cells(c("Epidish_BPRNACan_Monocytes", "Monocytes_Subgroup.2", "Quantiseq_B.cells"))
#' extract_cells("Method_Signature_Mesenchymal.cells", cells_extra = "Mesenchymal.cells")
#'
#' @export
#'
extract_cells = function(groups, cells_extra = NULL){
  # Cell type vocabulary is sourced from multideconv so it stays in sync with that
  # package's deconvolution column naming (see multideconv::get_cell_type_nomenclature)
  names_cells = get_cell_type_nomenclature(cells_extra = cells_extra)
  names_cells = gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", names_cells) # Matched literally ("." or "+" are not regex symbols)

  # Create regex to capture base cell type + cluster
  # Matches known cell type, optionally followed by _Subgroup.x (or _Subgroup.x.Iteration.x from multideconv < 0.2.0)
  regex_pattern <- paste0(
    "(",
    paste(names_cells, collapse = "|"),
    ")",
    "((?:[._](?:Subgroup|Iteration)\\.[0-9]+)*)"
  )


  normalized_names <- vapply(groups, function(x) {
    m <- regexpr(regex_pattern, x, perl = TRUE)
    if (m != -1) {
      matched <- regmatches(x, m)
      # remove Iteration / Subgroup artifacts
      gsub("([._](Subgroup|Iteration)\\.[0-9]+)+", "", matched)
    } else {
      NA_character_
    }
  }, character(1), USE.NAMES = FALSE)

  normalized_names <- unique(normalized_names[!is.na(normalized_names)])

  return(normalized_names)
}

#' Merge highly correlated TF modules
#'
#' Merges TF modules whose scores (eigengenes) are too similar, repeating until no pair of modules has a correlation
#' above \code{corr}: at each step the most correlated pair above \code{corr} is merged, and the correlations are
#' computed again. The merged module keeps the name and color of the first module of the pair; its score is the mean
#' of the scores of all the original modules it contains, scaled (mean 0, standard deviation 1).
#'
#' @param data A data frame of TF module eigengenes (samples x modules), with columns named \code{ME<color>}.
#' @param colors A character vector with the module color label of each TF.
#' @param corr Numeric. Correlation threshold above which two modules are merged.
#' @param method Character. Correlation used: "pearson" or "spearman" (default).
#'
#' @return A list of two elements:
#' \itemize{
#'   \item \code{[[1]]}: Data frame of merged module eigengenes (samples x modules).
#'   \item \code{[[2]]}: Updated character vector of module color labels after merging.
#' }
#'
#' @keywords internal
mergeModules = function(data, colors, corr, method = "spearman"){
  original = data
  members = stats::setNames(as.list(colnames(data)), colnames(data)) # Original modules contained in each module

  repeat{
    if(ncol(data) < 2) break
    cc = stats::cor(data, method = method, use = "pairwise.complete.obs")
    cc[lower.tri(cc, diag = TRUE)] = NA # Each pair once
    if(all(is.na(cc)) || max(cc, na.rm = TRUE) <= corr) break

    best = which(cc == max(cc, na.rm = TRUE), arr.ind = TRUE)[1, ]
    module1 = colnames(data)[best[1]]
    module2 = colnames(data)[best[2]]
    colors[colors %in% c(substring(module1, 3), substring(module2, 3))] = substring(module1, 3)
    members[[module1]] = c(members[[module1]], members[[module2]])
    members[[module2]] = NULL
    data[[module1]] = as.numeric(scale(rowMeans(original[, members[[module1]], drop = FALSE]))) # Mean of the original module scores, scaled
    data[[module2]] = NULL
  }

  return(list(data, colors))
}

#' Calculate dendrogram cluster assignments
#'
#' Cuts each cell-type dendrogram into clusters (the cell groups of a TF module) in two steps:
#' \enumerate{
#'   \item the dendrogram is cut with dynamic tree cutting (\code{dynamicTreeCut::cutreeDynamic()}, \code{method = "tree"});
#'   \item each cluster is made to cover a single range of correlation values with the module. With \code{method = "tree"},
#'     a cluster can be a branch from which smaller branches were split off as other clusters, so its features can lie in
#'     separated ranges of correlation (e.g. 0.14, 0.17-0.19 and 0.26-0.39, with other clusters in between). The assigned
#'     features are sorted by their correlation, and each run of consecutive features with the same cluster becomes a
#'     piece. A piece with fewer than \code{min_cluster_size} features is joined to the neighbouring piece (just below or
#'     just above in correlation) with the smallest gap between their correlations, the smallest gap first, until every
#'     piece has at least \code{min_cluster_size} features. The pieces are the final clusters, numbered from the lowest to
#'     the highest correlation.
#' }
#' Features left unassigned by the dynamic tree cut stay unassigned. Step 2 needs the correlations stored in the
#' dendrograms by \code{identify.cell.groups()}; for dendrograms without them, only step 1 is done (with a message).
#'
#' @param cell.group.dendrogram A list of \code{hclust} objects (or objects convertible to \code{hclust}),
#'   one per TF module, as returned by \code{identify.cell.groups()}.
#' @param deep_split Logical. Passed to \code{dynamicTreeCut::cutreeDynamic()}'s
#'   \code{deepSplit} argument. With \code{method = "tree"} it is a switch, not a level: \code{TRUE}
#'   (default; any non-zero number has the same effect) splits clusters more finely than \code{FALSE}.
#' @param min_cluster_size Integer. Minimum number of features per cluster: passed to
#'   \code{dynamicTreeCut::cutreeDynamic()}'s \code{minClusterSize} argument, and used in step 2. Default is 3.
#'
#' @return A list of integer cluster label vectors, one per dendrogram, in the
#'   same order as \code{cell.group.dendrogram} and named like it, each in the order of the dendrogram labels
#'   (label \code{0} marks unassigned elements).
#'
#' @keywords internal
calculate_dendrogram_cuts = function(cell.group.dendrogram, deep_split = TRUE, min_cluster_size = 3){

  cuts = list()
  no_values = c()
  for (i in 1:length(cell.group.dendrogram)) {
    hc <- if(inherits(cell.group.dendrogram[[i]], "hclust")) cell.group.dendrogram[[i]] else as.hclust(cell.group.dendrogram[[i]])

    labels <- dynamicTreeCut::cutreeDynamic(
      dendro         = hc,
      method         = "tree",
      deepSplit      = deep_split,
      minClusterSize = min_cluster_size,
      verbose        = 0
    )

    # Each cluster must be one range of correlation values: split clusters into their separated pieces,
    # and join pieces smaller than min_cluster_size to their closest neighbour
    values <- hc$correlation
    if(is.null(values) || length(values) != length(labels)){
      no_values <- c(no_values, if(is.null(names(cell.group.dendrogram))) i else names(cell.group.dendrogram)[i])
    } else if(any(labels > 0)){
      assigned <- which(labels > 0)
      assigned <- assigned[order(values[assigned])]
      v <- unname(values[assigned])
      runs <- rle(labels[assigned])
      ends <- cumsum(runs$lengths)
      pieces <- lapply(seq_along(ends), function(j) (ends[j] - runs$lengths[j] + 1):ends[j]) # Positions in the sorted assigned features
      repeat{
        sizes <- lengths(pieces)
        if(length(pieces) < 2 || all(sizes >= min_cluster_size)) break
        # Gap between each pair of neighbouring pieces; a pair is a candidate if one of the two pieces is small
        gaps <- sapply(seq_len(length(pieces) - 1), function(j) v[min(pieces[[j + 1]])] - v[max(pieces[[j]])])
        candidates <- which(sizes[-length(sizes)] < min_cluster_size | sizes[-1] < min_cluster_size)
        j <- candidates[which.min(gaps[candidates])]
        pieces[[j]] <- c(pieces[[j]], pieces[[j + 1]])
        pieces[[j + 1]] <- NULL
      }
      labels[] <- 0
      for(j in seq_along(pieces)) labels[assigned[pieces[[j]]]] <- j # Numbered from the lowest to the highest correlation
    }

    cuts[[i]] <- labels
  }
  names(cuts) <- names(cell.group.dendrogram)

  if(length(no_values) > 0){
    message("Dendrograms without the correlations of identify.cell.groups() (only the dynamic tree cut is used, clusters may cover separated ranges): ",
            paste(no_values, collapse = ", "))
  }

  return(cuts)

}

#' Plot cell-group dendrograms colored by cluster
#'
#' Draws each cell-type dendrogram with its branches colored by the cluster
#' assignments from \code{calculate_dendrogram_cuts()} and saves all of them
#' in a single PDF. The leaves are drawn in the order of their correlation with the module as far as the tree allows
#' (the two sides of a branch are swapped, the tree itself is unchanged), so that each cluster, a range of
#' correlation values, is drawn in one place when possible. A cluster that the tree splits into separate blocks of
#' leaves is drawn in each of them. The number of the cluster is written next to each of its blocks: it is the
#' number in the cell group name (\code{Dendrogram_<module>_group_<number>}).
#' Unassigned features are left uncolored.
#'
#' @param cell.group.dendrogram A named list of \code{hclust} objects (or objects convertible to \code{hclust}),
#'   one per TF module, as returned by \code{identify.cell.groups()} (the correlations they store are used to order the
#'   leaves; without them, the leaves keep the order of the tree).
#' @param cuts_per_dendrogram A list of integer cluster label vectors, one per dendrogram,
#'   as returned by \code{calculate_dendrogram_cuts()}.
#' @param file_name Optional character. Suffix of the output file
#'   (\code{Results/Dendrogram_color_clusters_<file_name>.pdf}). If \code{NULL} (default),
#'   nothing is saved.
#'
#' @return Called for its side effect (saves a PDF); returns \code{NULL} invisibly.
#'
#' @keywords internal
plot_dendrogram_clusters = function(cell.group.dendrogram, cuts_per_dendrogram, file_name = NULL) {
  plots <- list()

  for (i in 1:length(cell.group.dendrogram)) {
    hc     <- if (inherits(cell.group.dendrogram[[i]], "hclust")) cell.group.dendrogram[[i]] else as.hclust(cell.group.dendrogram[[i]])
    dend   <- as.dendrogram(hc)
    if (!is.null(hc$correlation)) dend <- stats::reorder(dend, wts = hc$correlation, agglo.FUN = mean) # Leaves in the order of their correlation, as far as the tree allows
    labels <- cuts_per_dendrogram[[i]]
    names(labels) <- hc$labels

    n_clusters <- max(labels)
    if (n_clusters <= 12) {
      cluster_colors <- RColorBrewer::brewer.pal(max(3, n_clusters), "Paired")[1:n_clusters]
    } else {
      cluster_colors <- colorspace::rainbow_hcl(n_clusters, c = 80, l = 55)
    }

    dend <- dendextend::color_branches(dend, clusters = labels[order.dendrogram(dend)],
                                       col = cluster_colors)
    dend <- dendextend::set(dend, "labels", rep("", length(labels(dend))))  # suppress leaf labels
    dend <- dendextend::set(dend, "branches_lwd", 2)

    # Cluster number next to each block of leaves of a cluster: cell groups are named Dendrogram_<module>_group_<number>
    leaf_clusters <- labels[order.dendrogram(dend)]
    blocks <- rle(unname(leaf_clusters))
    block_end <- cumsum(blocks$lengths)
    block_pos <- block_end - (blocks$lengths - 1) / 2 # Middle of each block
    keep <- blocks$values > 0
    ids <- blocks$values[keep]
    ids_pos <- block_pos[keep]

    p <- dendextend::as.ggdend(dend) |>
      ggplot2::ggplot() +
      ggplot2::annotate("text", x = ids_pos, y = -0.02 * attr(dend, "height"), label = ids,
                        colour = cluster_colors[match(ids, unique(leaf_clusters[leaf_clusters > 0]))], # branch colors follow the order of the clusters along the leaves
                        hjust = 0, size = 4.5, fontface = "bold") +
      ggplot2::coord_flip() +
      ggplot2::scale_y_reverse() +
      ggplot2::ggtitle(paste("Dendrogram:", names(cell.group.dendrogram)[i]),
                       subtitle = paste0("Numbers: clusters (cell groups Dendrogram_", names(cell.group.dendrogram)[i], "_group_<number>)")) +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text   = ggplot2::element_blank(),
                     axis.ticks  = ggplot2::element_blank(),
                     plot.title  = ggplot2::element_text(hjust = 0.5, face = "bold"),
                     plot.subtitle = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::labs(x = NULL, y = "Height")

    plots[[i]] <- p
  }

  if (!is.null(file_name)) {
    out_path <- file.path("Results", paste0("Dendrogram_color_clusters_", file_name, ".pdf"))
    pdf(out_path, width = 10, height = 12)
    for (p in plots) print(p)
    dev.off()
    message("Dendrogram plots saved to: ", out_path)
  }

}

#' Run Reactome pathway enrichment for a single TF module
#'
#' Given the hub TFs of a module, extracts their target genes from the TF-gene
#' universe, keeps the 20% most variable targets, and runs Reactome over-representation
#' analysis (\code{ReactomePA::enrichPathway()}, human) using all the genes of the expression matrix as the
#' background universe. Gene symbols are converted to Entrez IDs with \code{org.Hs.eg.db}; symbols without an Entrez ID
#' are left out (a message gives how many). Pathways with a Benjamini-Hochberg adjusted p-value < 0.05 are kept.
#' TF names are matched after \code{make.names()} on both sides, so the syntactic names returned by
#' \code{compute.TFs.activity()} (e.g. \code{NKX2.1}) match the names of the network (\code{NKX2-1}).
#'
#' The result describes the genes regulated by the hub TFs of the module: any set of TFs has hundreds of targets that
#' carry the pathways those TFs regulate, so random sets of TFs also give enriched pathways. It should be read as what
#' the hub TFs of the module regulate, not as evidence that the module itself is special.
#'
#' @param tpm.counts A numeric matrix of normalized expression values (genes x samples), with gene symbols as row names.
#' @param module_color Character. Name of the TF module (color label) to analyze.
#' @param hub_genes A list as returned by \code{identify_hub_TFs()}, where the first
#'   element maps module names to vectors of hub TF gene symbols.
#' @param tfs_universe A data frame of TF-target interactions with at minimum columns
#'   \code{source} (TF) and \code{target} (gene).
#' @param min_genes Integer. Minimum number of genes to test (20% most variable targets); with fewer, no enrichment is
#'   computed and a message is given. Default is 5.
#'
#' @return A \code{ReactomePA} enrichResult object restricted to pathways with
#'   adjusted p-value < 0.05, or \code{NULL} if no enrichment is found, if none of the
#'   targets is in the row names of \code{tpm.counts}, or if fewer than \code{min_genes} genes would be tested.
#'
#' @keywords internal
module_enrich = function(tpm.counts, module_color, hub_genes, tfs_universe, min_genes = 5){
  targets = tfs_universe[make.names(tfs_universe$source) %in% make.names(hub_genes[[1]][[module_color]]),] #Extract targets from TFs (syntactic names too, e.g. NKX2.1 for NKX2-1)
  targets = unique(targets$target) #Keep only unique targets from TFs

  targets_genes = tpm.counts[rownames(tpm.counts)%in%targets, , drop = F] #Extract gene expression from targets
  if(nrow(targets_genes) == 0){
    return(NULL)
  }
  n_tested = round(0.2*nrow(targets_genes)) #Keep only highly variable targets (20%)
  if(n_tested < min_genes){
    message("Module ", module_color, ": only ", n_tested, " gene(s) would be tested (20% most variable of its ", nrow(targets_genes),
            " targets in the data); at least ", min_genes, " are needed, no enrichment computed")
    return(NULL)
  }
  targets_genes = targets_genes[order(matrixStats::rowVars(targets_genes), decreasing = T), , drop = F][1:n_tested, , drop = F]

  # Messages of the annotation and enrichment packages are not shown
  entrz <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = rownames(targets_genes), columns = "ENTREZID", keytype = "SYMBOL")) #Change to EntrezID
  universe = suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = rownames(tpm.counts), columns = "ENTREZID", keytype = "SYMBOL")) #Change to EntrezID
  unmapped_tested = length(unique(entrz$SYMBOL[is.na(entrz$ENTREZID)]))
  unmapped_universe = length(unique(universe$SYMBOL[is.na(universe$ENTREZID)]))
  if(unmapped_tested > 0 || unmapped_universe > 0){
    message("Module ", module_color, ": ", unmapped_tested, " of the ", n_tested, " tested genes and ", unmapped_universe, " of the ", nrow(tpm.counts),
            " background genes have no Entrez ID and are left out")
  }

  reac <- suppressMessages(ReactomePA::enrichPathway(gene    = unique(stats::na.omit(entrz$ENTREZID)),
                                                     organism     = 'human',
                                                     universe = unique(stats::na.omit(universe$ENTREZID)),
                                                     pvalueCutoff = 0.05))
  if(is.null(reac)){
    return(NULL)
  }

  reac@result = reac@result[reac@result$p.adjust<0.05,]

  if(nrow(reac@result) == 0){
    return(NULL)
  }
  return(reac)
}

#' Compute composite score for cell groups
#'
#' Computes a composite score by performing Canonical Correlation Analysis (CCA)
#' between the deconvolution features of a cell group and the TFs of its module.
#' The score is the weighted sum of the (scaled) features of the cell group, with the weights of the first
#' canonical component. As the sign of CCA weights is arbitrary, the score is oriented so that it is aligned
#' with the TF module: a positive score means that the sample is on the active side of the module for that cell
#' group. The score therefore increases with the features of the cell group when they are positively related
#' to the module, and decreases with them when they are negatively related. The relation of the cell group to
#' its module (mean correlation of its features with the module score) is returned to translate the score
#' into cell abundance.
#'
#' @param cell_group A numeric matrix of cell deconvolution features for a cell group (samples x features), with sample
#'   names as row names. The module score and TF activities of \code{tfs.module.network} are matched to these samples by
#'   name; every sample must be in the network.
#' @param module_group Character. Name (color) of the TF module the cell group was built from; must match a
#'   column name of the module matrix in \code{tfs.module.network} exactly.
#' @param tfs.module.network Output of compute.WTCNA() (single- or multi-cohort).
#' @param batch Optional vector with the cohort of each sample, in the row order of \code{cell_group}. It is treated as
#'   categorical: per-batch means are regressed out of the cell group features, the module score and the module TFs
#'   before the CCA.
#'
#' @return An unnamed list of three elements:
#' \itemize{
#'   \item \code{[[1]]}: Numeric matrix (samples x 1) with the composite score, i.e. the scaled cell group
#'     features projected onto the first canonical component, oriented to be aligned with the TF module.
#'   \item \code{[[2]]}: Projection parameters used to score new samples: \code{xcoef} (canonical weights of
#'     the first component), \code{train_means} and \code{train_sds} (column means/SDs used for scaling).
#'     The sign of \code{xcoef} follows the orientation of the score.
#'     \code{train_means} is \code{NULL} when batch correction was applied: new samples are then centred on
#'     their own means, as each training cohort was.
#'   \item \code{[[3]]}: Relation of the cell group to its TF module: the mean correlation of its features with
#'     the module score (positive or negative).
#' }
#'
#' @details
#' The CCA is reliable only when the number of variables (features of the cell group plus TFs of the module) is well
#' below the number of samples. When it reaches the number of samples, the first canonical correlation is 1 whatever
#' the data, the weights are not unique and the score is decided by numerical noise (it can change completely with
#' the sample order or a tiny change in the data). Large modules with few samples (e.g. cross-validation folds) are the
#' most affected.
#'
#' @examples
#' data("network.tuto")
#' data("deconv_subgroups.tuto")
#' cell_groups <- construct_cell_groups(network.tuto, deconv_subgroups.tuto)
#' group <- cell_groups$Module_relation[1, ]
#' features <- deconv_subgroups.tuto[[1]][, cell_groups$Composition[[group$cell_group]], drop = FALSE]
#' score <- compute_composite_score(features, group$module, network.tuto)
#' head(score[[1]])
#'
#' @export
#'
compute_composite_score = function(cell_group, module_group, tfs.module.network, batch = NULL){

  modules = tfs.module.network[["TFs module matrix"]]
  tfs_all = tfs.module.network[["TFs_matrix"]]
  if(is.null(modules) || is.null(tfs_all) || is.null(tfs.module.network[["TFs per module"]])){
    stop("tfs.module.network must be the output of compute.WTCNA()")
  }
  if(!module_group %in% colnames(modules) || is.null(tfs.module.network[["TFs per module"]][[module_group]])){
    stop("module_group '", module_group, "' is not a TF module of tfs.module.network")
  }
  if(!is.null(batch) && length(batch) != nrow(cell_group)){
    stop("batch must have one value per sample of cell_group (", nrow(cell_group), ")")
  }

  # Module score and TF activities in the sample order of cell_group (matched by name)
  if(is.list(tfs_all) && !is.data.frame(tfs_all)) tfs_all <- do.call(rbind, unname(lapply(tfs_all, as.matrix))) # Multi-cohort network: stacked cohort by cohort; unname() keeps the sample names
  if(is.null(rownames(cell_group))) stop("cell_group must have the sample names as row names")
  idx_tfs = match(make.names(rownames(cell_group)), make.names(rownames(tfs_all)))
  idx_modules = match(make.names(rownames(cell_group)), make.names(rownames(modules)))
  if(anyNA(idx_tfs) || anyNA(idx_modules)){
    missing_samples = rownames(cell_group)[is.na(idx_tfs) | is.na(idx_modules)]
    stop(length(missing_samples), " samples of cell_group are not in tfs.module.network: ", paste(utils::head(missing_samples, 5), collapse = ", "))
  }
  tfs_all <- tfs_all[idx_tfs, , drop = F]
  modules <- modules[idx_modules, , drop = F]

  tf_per_module_matrix = tfs_all[, colnames(tfs_all) %in% tfs.module.network[["TFs per module"]][[module_group]]]
  tf_module_matrix = modules[, colnames(modules) == module_group, drop = F]

  # ---- Regress out batch if provided ----
  if(!is.null(batch) && length(unique(batch)) > 1){
    batch <- as.factor(batch)  # Categorical: regress out per-batch means
    cell_group <- apply(cell_group, 2, function(x) residuals(lm(x ~ batch)))
    tf_module_matrix <- apply(tf_module_matrix, 2, function(x) residuals(lm(x ~ batch)))
    tf_per_module_matrix <- apply(tf_per_module_matrix, 2, function(x) residuals(lm(x ~ batch)))
  }

  # Relation of the cell group to its TF module: mean correlation of its features with the module score
  relation = mean(stats::cor(cell_group, tf_module_matrix))

  cell_group_matrix = scale(cell_group)

  # CCA analysis: Note 'cell_group_matrix' is already scale
  cca_result <- stats::cancor(cell_group_matrix, tf_per_module_matrix) ## Now this is done per tfs_per_module (tfs are already scale check compute.TFs.activity())

  # Apply canonical weights to each matrix
  cell_group <- cell_group[, rownames(cca_result$xcoef), drop = F] #Ensure no features have been discard on the way and if yes, subset the matrix (CCA can discard collinear or zero variance features)

  scaled_obj <- scale(as.matrix(cell_group)) ### Why do you scale again? Because the CCA weights are based on the scale matrix
  train_means <- colMeans(cell_group)
  train_sds   <- apply(cell_group, 2, sd)
  if(!is.null(batch) && length(unique(batch)) > 1) train_means <- NULL # Batch-corrected: new cohorts are centred on their own means

  ### xcoef[,1] corresponds to the most correlated linear component (we used the scale cell_group_matrix because the coef came from the scale matrix)
  xcoef = cca_result$xcoef[, 1, drop = F]
  weighted_cell_group_matrix <- scaled_obj %*% xcoef

  # The CCA weights are only defined up to their sign (w and -w are equally good solutions, and which one is returned
  # depends e.g. on the sample order): orient the score so that it is aligned with the TF module. It increases with the
  # features of the cell group when they are positively related to the module, and decreases with them when negatively related
  # Positive CellTFusion score = higher TF module activity
  # ---------------------------------------------------------
  module_score <- as.numeric(tf_module_matrix[, 1])

  orientation_cor <- cor(
    as.numeric(weighted_cell_group_matrix),
    module_score,
    use = "complete.obs"
  )

  if(isTRUE(orientation_cor < 0)){
    xcoef <- -xcoef
    weighted_cell_group_matrix <- -weighted_cell_group_matrix
  }

  projection_params = list(
    xcoef       = xcoef,
    train_means = train_means,
    train_sds   = train_sds
  )

  return(list(weighted_cell_group_matrix, projection_params, relation))

}

#' Recursively retrieve the original deconvolution features of a subgroup
#'
#' Replaces a deconvolution subgroup by its member features, going down through nested subgroups, and returns the
#' original deconvolution feature names (e.g. \code{B.cells_Subgroup.1} gives \code{CBSX_BPRNACan_B.cells},
#' \code{CBSX_BPRNACan3DProMet_B.cells}, ...). A name that is not a subgroup is returned unchanged. The cell types are
#' read from the returned names afterwards with \code{extract_cells()}.
#'
#' @param subgroup_name Character. A single deconvolution feature name: a subgroup or an original feature.
#' @param cell_subgroups A named list mapping subgroup names to vectors of member
#'   names (which can themselves be subgroup names), e.g. the flattened \code{Deconvolution subgroups composition}
#'   of \code{multideconv::compute.deconvolution.analysis()}.
#'
#' @return A character vector of the unique original deconvolution feature names contained in
#'   \code{subgroup_name} (the name itself if it is not a subgroup).
#'
#' @keywords internal
get_all_cells <- function(subgroup_name, cell_subgroups) {
  if (subgroup_name %in% names(cell_subgroups)) { #Check if subgroup_name is key in cell_subgroups
    # If the subgroup contains further subgroups, retrieve their base elements
    nested_cells <- unlist(lapply(cell_subgroups[[subgroup_name]], get_all_cells, cell_subgroups = cell_subgroups))
    return(unique(nested_cells)) # Only return base elements
  } else {
    # If subgroup_name is not a key in cell_subgroups, it is considered a base composition
    return(subgroup_name)
  }
}

#' Project test-set cell group scores using training CCA parameters
#'
#' Scales a test-set cell group matrix using the mean and standard deviation
#' stored from training, then projects it onto the first canonical component
#' learned during training. All the features of the cell group are needed: a score computed from part of them
#' (the others set to their training mean) can be very different from the training score, even reversed, because
#' the CCA weights of correlated features can have opposite signs.
#'
#' @param cell_group A numeric matrix of cell deconvolution features for test
#'   samples (samples x features). It must contain every feature of the cell group (the row names of \code{xcoef});
#'   other columns are ignored.
#' @param projection_params A list containing:
#'   \describe{
#'     \item{xcoef}{Named numeric matrix of CCA canonical weights (features x 1).}
#'     \item{train_means}{Named numeric vector of training column means, or \code{NULL} if the training used
#'       batch correction, in which case the test data are centred on their own column means.}
#'     \item{train_sds}{Named numeric vector of training column standard deviations.}
#'   }
#'
#' @return A numeric matrix (samples x 1) of projected canonical scores for the
#'   test cohort, or \code{NULL} (with a warning) if some features of the cell group are missing.
#'
#' @keywords internal
compute.test.score = function(cell_group, projection_params){

  xcoef       = projection_params$xcoef
  train_means = projection_params$train_means
  train_sds   = projection_params$train_sds

  # Every feature of the cell group is needed (a partial score is not comparable with the training score)
  features = rownames(xcoef)
  missing_features = setdiff(features, colnames(cell_group))
  if(length(missing_features) > 0){
    warning("Features of the cell group missing in the test data (", paste(utils::head(missing_features, 5), collapse = ", "), "): no score computed")
    return(NULL)
  }

  cell_group  = cell_group[, features, drop = F]
  train_means = train_means[features]
  train_sds   = train_sds[features]

  # Scale using the training means/SDs, as in training
  if(is.null(train_means)) train_means = colMeans(cell_group) # Batch-corrected training: centre the test cohort on its own means
  cell_group_scaled = scale(cell_group, center = train_means, scale = train_sds)

  return(as.matrix(cell_group_scaled) %*% xcoef)
}

#' Prepare the inputs of a score-trait test
#'
#' Internal helper of the \code{scores.*} tests. Takes the score matrix (first element of \code{scores}), matches the
#' rows of \code{coldata} to its samples by name, leaves out samples with a missing trait (with a message) and checks
#' that the trait is a grouping variable.
#'
#' @param scores A list whose first element is a samples x features score matrix, the matrix itself, or the output of
#'   \code{compute.latent_factors()} (its element \code{Z} is used).
#' @param coldata A data frame (or tibble) of sample annotations.
#' @param trait Character. Name of the column of \code{coldata} to test.
#' @param binary Logical. If TRUE, the trait must have exactly two groups.
#' @param test Character. Name of the test, used in the messages.
#'
#' @return A list with \code{scores} (the input as a list), \code{X} (score matrix of the samples with a trait value)
#'   and \code{trait} (factor of their trait values).
#'
#' @keywords internal
prepare_trait_test = function(scores, coldata, trait, binary = FALSE, test = "test"){
  if(is.matrix(scores) || is.data.frame(scores)) scores = list(scores)
  if(is.list(scores) && !is.null(scores$Z)) scores = list(scores$Z) # Output of compute.latent_factors(): the factor scores
  if(!is.list(scores) || !(is.matrix(scores[[1]]) || is.data.frame(scores[[1]]))){
    stop("scores must be a samples x features score matrix, or a list whose first element is that matrix")
  }
  X = scores[[1]]
  coldata = as.data.frame(coldata) # Tibbles are accepted
  if(!trait %in% colnames(coldata)) stop("trait '", trait, "' is not a column of coldata")

  # Samples matched by name (coldata without sample names: same order as the scores)
  idx = if(!is.null(rownames(X))) match(make.names(rownames(X)), make.names(rownames(coldata))) else NA
  if(!anyNA(idx)){
    coldata = coldata[idx, , drop = FALSE]
  } else if(.row_names_info(coldata) < 0 && nrow(coldata) == nrow(X)){
    message("coldata has no sample names: its rows are assumed to be in the order of the scores")
  } else {
    stop("The samples of the scores are not all in the row names of coldata (",
         paste(utils::head(rownames(X)[is.na(idx)], 5), collapse = ", "), ")")
  }

  values = coldata[[trait]]
  missing_trait = is.na(values)
  if(any(missing_trait)){
    message(sum(missing_trait), " samples without a value for '", trait, "' are left out")
  }
  if(is.numeric(values) && length(unique(values[!missing_trait])) > 10){
    stop("'", trait, "' looks continuous (", length(unique(values[!missing_trait])), " distinct values): the ", test,
         " compares groups of samples; give a categorical trait")
  }
  values = factor(values[!missing_trait])
  n_levels = nlevels(values)
  if(binary && n_levels != 2) stop("The ", test, " requires a binary trait (exactly two groups); '", trait, "' has ", n_levels)
  if(n_levels < 2) stop("'", trait, "' has a single group: nothing to compare")

  list(scores = scores, X = as.data.frame(X)[!missing_trait, , drop = FALSE], trait = values)
}

#' Build the result of a score-trait test
#'
#' Internal helper of the \code{scores.*} tests: adjusts the p-values, selects the significant features and returns
#' the significant part of \code{scores} with a table of the tests of all features.
#'
#' @param input The output of \code{prepare_trait_test()}.
#' @param statistic Numeric vector with the test statistic of each feature.
#' @param p Numeric vector with the p-value of each feature.
#' @param test Character. Name of the test.
#' @param statistic_name Character. Name of the statistic (e.g. "t").
#' @param pval Numeric. Significance threshold for the adjusted p-values.
#' @param p_adjust Character. Method of \code{stats::p.adjust()}.
#'
#' @return A list with \code{result} (see the \code{scores.*} tests), \code{sig} (indices of the significant
#'   features) and \code{adjusted_text} (adjusted p-value of each feature for the plot subtitles, empty without
#'   adjustment).
#'
#' @keywords internal
finish_trait_test = function(input, statistic, p, test, statistic_name, pval, p_adjust){
  scores = input$scores
  features = colnames(input$X)
  p_adjusted = stats::p.adjust(p, method = p_adjust)
  sig = which(!is.na(p_adjusted) & p_adjusted <= pval)

  statistics = data.frame(feature = features, test = test, statistic_name = statistic_name, statistic = statistic,
                          p_value = p, p_adjusted = p_adjusted, significant = seq_along(features) %in% sig,
                          stringsAsFactors = FALSE)

  result = list()
  result[[1]] = scores[[1]][, sig, drop = FALSE]
  if (length(scores) >= 2) result[[2]] = scores[[2]][sig]
  if (length(scores) >= 3) result[[3]] = scores[[3]][sig]
  if (length(scores) >= 4 && is.data.frame(scores[[4]])) result[[4]] = scores[[4]][sig, , drop = FALSE]
  if (!is.null(names(scores))) names(result) = names(scores)[seq_along(result)]
  result[["Statistics"]] = statistics

  threshold = paste0(if (p_adjust == "none") "p-value" else paste0("adjusted p-value (", p_adjust, ")"), " <= ", pval)
  if (length(sig) == 0) {
    message("No significant feature after ", test, " (", threshold, ")")
  } else {
    message(length(sig), " significant feature(s) after ", test, " (", threshold, "): ",
            paste(utils::head(features[sig], 10), collapse = ", "), if (length(sig) > 10) ", ...")
  }

  # Adjusted p-value for the plot subtitles (empty without adjustment)
  adjusted_text = if (p_adjust == "none") rep("", length(p)) else
    paste0(" | adjusted p (", p_adjust, ") ", ifelse(p_adjusted < 0.001, "< 0.001", paste0("= ", round(p_adjusted, 4))))

  list(result = result, sig = sig, adjusted_text = adjusted_text)
}

#' T-test for cell group comparisons
#'
#' Performs a two-sample t-test (Welch's, unequal variances) comparing cell group scores between two groups
#' of a binary trait.
#' Significant features are plotted as boxplots and saved as PDF files in the "Results/" directory
#' (\code{Results/Ttest_<trait>_<feature>.pdf}).
#'
#' @param scores A samples x features score matrix, a list whose first element is that matrix (e.g. the output of
#'   \code{construct_cell_groups()}), or the output of \code{compute.latent_factors()} (its factor scores \code{Z} are
#'   tested).
#' @param coldata A data frame (or tibble) containing sample-level annotations including the trait to test. Its rows
#'   are matched to the samples of the scores by row name; without sample names, they must be in the order of the
#'   scores. Samples with a missing trait are left out (with a message).
#' @param trait Character. Name of the column in `coldata` used as the grouping variable (categorical, two groups).
#' @param pval Numeric. Threshold for significance: a feature is significant if its (adjusted) p-value is
#'   \code{<= pval} (default = 0.05).
#' @param p_adjust Character. Method to correct the p-values for multiple testing across the features (see
#'   \code{stats::p.adjust()}, e.g. \code{"BH"}). Default \code{"none"}.
#'
#' @return A list with the significant features: its first element is the score matrix restricted to them (no
#'   column if none is significant), followed by the matching parts of the other elements of \code{scores} (e.g.
#'   \code{Composition}, \code{Weights} and \code{Module_relation} of \code{construct_cell_groups()}), with their names.
#'   The last element, \code{Statistics}, is a data frame with the test of every feature (\code{feature}, \code{test},
#'   \code{statistic_name}, \code{statistic}, \code{p_value}, \code{p_adjusted}, \code{significant}).
#'
#' @examples
#' data("network.tuto")
#' data("deconv_subgroups.tuto")
#' data("traitdata.tuto")
#' cell_groups <- construct_cell_groups(network.tuto, deconv_subgroups.tuto)
#' res <- scores.ttest(cell_groups$Cell_groups[, 1:5], traitdata.tuto, "binaryResponse")
#' res$Statistics
#'
#' @export
#'
scores.ttest <- function(scores, coldata, trait, pval = 0.05, p_adjust = "none") {
  input = prepare_trait_test(scores, coldata, trait, binary = TRUE, test = "Welch's t-test")
  X = input$X

  tests = lapply(seq_len(ncol(X)), function(j) tryCatch(stats::t.test(X[, j] ~ input$trait, var.equal = FALSE), error = function(e) NULL))
  res = finish_trait_test(input, sapply(tests, function(t) if (is.null(t)) NA else unname(t$statistic)),
                          sapply(tests, function(t) if (is.null(t)) NA else t$p.value), "Welch's t-test", "t", pval, p_adjust)

  for (j in res$sig) {
    data = data.frame(Value = X[, j], Trait = input$trait)

    n_labels  <- data %>% dplyr::count(Trait) %>%
      dplyr::mutate(label = paste0(as.character(Trait), "\n(n=", n, ")"))
    label_map <- setNames(n_labels$label, as.character(n_labels$Trait))

    stat_test <- data %>%
      rstatix::t_test(Value ~ Trait, var.equal = FALSE) %>%
      rstatix::add_xy_position(x = "Trait") %>%
      dplyr::mutate(p.format = ifelse(p < 0.001, "< 0.001", as.character(round(p, 4))))

    pdf(paste0("Results/Ttest_", trait, "_", colnames(X)[j], ".pdf"),
        width = 7, height = 6)
    print(
      ggplot2::ggplot(data, ggplot2::aes(x = Trait, y = Value, fill = Trait)) +
        ggplot2::geom_violin(alpha = 0.35, trim = FALSE, color = NA) +
        ggplot2::geom_boxplot(width = 0.14, outlier.shape = NA, alpha = 0.85,
                              color = "gray25", linewidth = 0.5) +
        ggplot2::geom_jitter(width = 0.07, size = 1.6, alpha = 0.45, color = "gray20") +
        ggpubr::stat_pvalue_manual(stat_test, label = "p = {p.format}",
                                   tip.length = 0.02, size = 4.5, bracket.size = 0.5,
                                   inherit.aes = FALSE) +
        ggplot2::scale_x_discrete(labels = label_map) +
        ggplot2::scale_fill_brewer(palette = "Set2") +
        ggplot2::labs(
          title    = colnames(X)[j],
          subtitle = paste0("Welch's t-test | ", trait, res$adjusted_text[j]),
          x        = NULL,
          y        = "Score"
        ) +
        ggplot2::theme_classic(base_size = 13) +
        ggplot2::theme(
          legend.position  = "none",
          axis.text.x      = ggplot2::element_text(size = 12, color = "black"),
          axis.text.y      = ggplot2::element_text(size = 11),
          axis.title.y     = ggplot2::element_text(size = 13),
          plot.title       = ggplot2::element_text(size = 14, face = "bold"),
          plot.subtitle    = ggplot2::element_text(size = 10, color = "gray45"),
          panel.grid.major.y = ggplot2::element_line(color = "gray92", linewidth = 0.4)
        )
    )
    dev.off()
  }

  return(res$result)
}

#' Kruskal-Wallis test for multi-group comparisons
#'
#' Performs a Kruskal-Wallis test to compare scores across multiple trait levels.
#' Significant results are visualized as annotated boxplots with Dunn post-hoc tests (BH-adjusted), saved as
#' \code{Results/Kruskal_<trait>_<feature>.pdf}.
#'
#' @inheritParams scores.ttest
#' @param trait Character. Name of the column in `coldata` used as the grouping variable (categorical, two or more groups).
#'
#' @return See \code{scores.ttest()}.
#'
#' @examples
#' data("network.tuto")
#' data("deconv_subgroups.tuto")
#' data("traitdata.tuto")
#' cell_groups <- construct_cell_groups(network.tuto, deconv_subgroups.tuto)
#' res <- scores.kruskal.test(cell_groups$Cell_groups[, 1:5], traitdata.tuto, "Enrollment.IC")
#' res$Statistics
#'
#' @export
#'
scores.kruskal.test <- function(scores, coldata, trait, pval = 0.05, p_adjust = "none") {
  input = prepare_trait_test(scores, coldata, trait, test = "Kruskal-Wallis test")
  X = input$X

  tests = lapply(seq_len(ncol(X)), function(j) tryCatch(stats::kruskal.test(X[, j], input$trait), error = function(e) NULL))
  res = finish_trait_test(input, sapply(tests, function(t) if (is.null(t)) NA else unname(t$statistic)),
                          sapply(tests, function(t) if (is.null(t)) NA else t$p.value), "Kruskal-Wallis test", "H", pval, p_adjust)

  for (j in res$sig) {
    data = data.frame(Value = X[, j], Trait = input$trait)

    # Dunn post-hoc test + positions for significance annotation
    pwc <- data %>%
      rstatix::dunn_test(Value ~ Trait, p.adjust.method = "BH") %>%
      rstatix::add_xy_position(x = "Trait")

    n_labels  <- data %>% dplyr::count(Trait) %>%
      dplyr::mutate(label = paste0(as.character(Trait), "\n(n=", n, ")"))
    label_map <- setNames(n_labels$label, as.character(n_labels$Trait))

    pdf(paste0("Results/Kruskal_", trait, "_", colnames(X)[j], ".pdf"),
        width = 8, height = 6)
    print(
      ggplot2::ggplot(data, ggplot2::aes(x = Trait, y = Value, fill = Trait)) +
        ggplot2::geom_violin(alpha = 0.35, trim = FALSE, color = NA) +
        ggplot2::geom_boxplot(width = 0.14, outlier.shape = NA, alpha = 0.85,
                              color = "gray25", linewidth = 0.5) +
        ggplot2::geom_jitter(width = 0.07, size = 1.6, alpha = 0.45, color = "gray20") +
        ggpubr::stat_pvalue_manual(pwc, hide.ns = TRUE, label = "p.adj.signif",
                                   tip.length = 0.02, size = 5, bracket.size = 0.5,
                                   inherit.aes = FALSE) +
        ggplot2::scale_x_discrete(labels = label_map) +
        ggplot2::scale_fill_brewer(palette = "Set2") +
        ggplot2::labs(
          title    = colnames(X)[j],
          subtitle = paste0("Kruskal-Wallis | Dunn post-hoc (BH) | ", trait, res$adjusted_text[j]),
          caption  = rstatix::get_pwc_label(pwc),
          x        = NULL,
          y        = "Score"
        ) +
        ggplot2::theme_classic(base_size = 13) +
        ggplot2::theme(
          legend.position  = "none",
          axis.text.x      = ggplot2::element_text(size = 12, color = "black"),
          axis.text.y      = ggplot2::element_text(size = 11),
          axis.title.y     = ggplot2::element_text(size = 13),
          plot.title       = ggplot2::element_text(size = 14, face = "bold"),
          plot.subtitle    = ggplot2::element_text(size = 10, color = "gray45"),
          plot.caption     = ggplot2::element_text(size = 9, color = "gray55"),
          panel.grid.major.y = ggplot2::element_line(color = "gray92", linewidth = 0.4)
        )
    )
    dev.off()
  }

  return(res$result)
}

#' Wilcoxon rank-sum test for binary traits
#'
#' Performs a Wilcoxon rank-sum (Mann-Whitney U) test comparing scores between
#' two levels of a binary clinical trait.
#' Significant features are plotted as boxplots and saved to the "Results/" folder
#' (\code{Results/Wilcoxon_<trait>_<feature>.pdf}).
#'
#' @inheritParams scores.ttest
#'
#' @return See \code{scores.ttest()}.
#'
#' @examples
#' data("network.tuto")
#' data("deconv_subgroups.tuto")
#' data("traitdata.tuto")
#' cell_groups <- construct_cell_groups(network.tuto, deconv_subgroups.tuto)
#' res <- scores.wilcox.test(cell_groups$Cell_groups[, 1:5], traitdata.tuto, "binaryResponse")
#' res$Statistics
#'
#' @export
#'
scores.wilcox.test <- function(scores, coldata, trait, pval = 0.05, p_adjust = "none") {
  input = prepare_trait_test(scores, coldata, trait, binary = TRUE, test = "Wilcoxon test")
  X = input$X

  tests = lapply(seq_len(ncol(X)), function(j) tryCatch(stats::wilcox.test(X[, j] ~ input$trait, exact = FALSE), error = function(e) NULL))
  res = finish_trait_test(input, sapply(tests, function(t) if (is.null(t)) NA else unname(t$statistic)),
                          sapply(tests, function(t) if (is.null(t)) NA else t$p.value), "Wilcoxon rank-sum test", "W", pval, p_adjust)

  for (j in res$sig) {
    data = data.frame(Value = X[, j], Trait = input$trait)

    n_labels  <- data %>% dplyr::count(Trait) %>%
      dplyr::mutate(label = paste0(as.character(Trait), "\n(n=", n, ")"))
    label_map <- setNames(n_labels$label, as.character(n_labels$Trait))

    stat_test <- data %>%
      rstatix::wilcox_test(Value ~ Trait, exact = FALSE) %>%
      rstatix::add_xy_position(x = "Trait") %>%
      dplyr::mutate(p.format = ifelse(p < 0.001, "< 0.001", as.character(round(p, 4))))

    pdf(paste0("Results/Wilcoxon_", trait, "_", colnames(X)[j], ".pdf"),
        width = 7, height = 6)
    print(
      ggplot2::ggplot(data, ggplot2::aes(x = Trait, y = Value, fill = Trait)) +
        ggplot2::geom_violin(alpha = 0.35, trim = FALSE, color = NA) +
        ggplot2::geom_boxplot(width = 0.14, outlier.shape = NA, alpha = 0.85,
                              color = "gray25", linewidth = 0.5) +
        ggplot2::geom_jitter(width = 0.07, size = 1.6, alpha = 0.45, color = "gray20") +
        ggpubr::stat_pvalue_manual(stat_test, label = "p = {p.format}",
                                   tip.length = 0.02, size = 4.5, bracket.size = 0.5,
                                   inherit.aes = FALSE) +
        ggplot2::scale_x_discrete(labels = label_map) +
        ggplot2::scale_fill_brewer(palette = "Set2") +
        ggplot2::labs(
          title    = colnames(X)[j],
          subtitle = paste0("Wilcoxon rank-sum test | ", trait, res$adjusted_text[j]),
          x        = NULL,
          y        = "Score"
        ) +
        ggplot2::theme_classic(base_size = 13) +
        ggplot2::theme(
          legend.position  = "none",
          axis.text.x      = ggplot2::element_text(size = 12, color = "black"),
          axis.text.y      = ggplot2::element_text(size = 11),
          axis.title.y     = ggplot2::element_text(size = 13),
          plot.title       = ggplot2::element_text(size = 14, face = "bold"),
          plot.subtitle    = ggplot2::element_text(size = 10, color = "gray45"),
          panel.grid.major.y = ggplot2::element_line(color = "gray92", linewidth = 0.4)
        )
    )
    dev.off()
  }

  return(res$result)
}

#' One-way ANOVA test for multi-group comparisons
#'
#' Performs one-way ANOVA to test for differences in scores across multiple levels of a trait.
#' Tukey post-hoc tests are used to identify pairwise differences and significance is visualized as annotated
#' boxplots (\code{Results/ANOVA_<trait>_<feature>.pdf}).
#'
#' @inheritParams scores.ttest
#' @param trait Character. Name of the column in `coldata` used as the grouping variable (categorical, two or more groups).
#'
#' @return See \code{scores.ttest()}.
#'
#' @examples
#' data("network.tuto")
#' data("deconv_subgroups.tuto")
#' data("traitdata.tuto")
#' cell_groups <- construct_cell_groups(network.tuto, deconv_subgroups.tuto)
#' res <- scores.anova.test(cell_groups$Cell_groups[, 1:5], traitdata.tuto, "Enrollment.IC")
#' res$Statistics
#'
#' @export
#'
scores.anova.test = function(scores, coldata, trait, pval = 0.05, p_adjust = "none"){
  input = prepare_trait_test(scores, coldata, trait, test = "ANOVA")
  X = input$X

  tests = lapply(seq_len(ncol(X)), function(j) tryCatch(summary(stats::aov(X[, j] ~ input$trait))[[1]][1, c("F value", "Pr(>F)")], error = function(e) NULL))
  res = finish_trait_test(input, sapply(tests, function(t) if (is.null(t)) NA else t[["F value"]]),
                          sapply(tests, function(t) if (is.null(t)) NA else t[["Pr(>F)"]]), "one-way ANOVA", "F", pval, p_adjust)

  for (j in res$sig) {
    data = data.frame(Value = X[, j], Trait = input$trait)

    # Tukey post-hoc test + position for significance annotation
    pwc <- data %>%
      rstatix::tukey_hsd(Value ~ Trait) %>%
      rstatix::add_xy_position(x = "Trait")

    n_labels  <- data %>% dplyr::count(Trait) %>%
      dplyr::mutate(label = paste0(as.character(Trait), "\n(n=", n, ")"))
    label_map <- setNames(n_labels$label, as.character(n_labels$Trait))

    pdf(paste0("Results/ANOVA_", trait, "_", colnames(X)[j], ".pdf"),
        width = 8, height = 6)
    print(
      ggplot2::ggplot(data, ggplot2::aes(x = Trait, y = Value, fill = Trait)) +
        ggplot2::geom_violin(alpha = 0.35, trim = FALSE, color = NA) +
        ggplot2::geom_boxplot(width = 0.14, outlier.shape = NA, alpha = 0.85,
                              color = "gray25", linewidth = 0.5) +
        ggplot2::geom_jitter(width = 0.07, size = 1.6, alpha = 0.45, color = "gray20") +
        ggpubr::stat_pvalue_manual(pwc, hide.ns = TRUE, label = "p.adj.signif",
                                   tip.length = 0.02, size = 5, bracket.size = 0.5,
                                   inherit.aes = FALSE) +
        ggplot2::scale_x_discrete(labels = label_map) +
        ggplot2::scale_fill_brewer(palette = "Set2") +
        ggplot2::labs(
          title    = colnames(X)[j],
          subtitle = paste0("One-way ANOVA | Tukey HSD | ", trait, res$adjusted_text[j]),
          caption  = rstatix::get_pwc_label(pwc),
          x        = NULL,
          y        = "Score"
        ) +
        ggplot2::theme_classic(base_size = 13) +
        ggplot2::theme(
          legend.position  = "none",
          axis.text.x      = ggplot2::element_text(size = 12, color = "black"),
          axis.text.y      = ggplot2::element_text(size = 11),
          axis.title.y     = ggplot2::element_text(size = 13),
          plot.title       = ggplot2::element_text(size = 14, face = "bold"),
          plot.subtitle    = ggplot2::element_text(size = 10, color = "gray45"),
          plot.caption     = ggplot2::element_text(size = 9, color = "gray55"),
          panel.grid.major.y = ggplot2::element_line(color = "gray92", linewidth = 0.4)
        )
    )
    dev.off()
  }

  return(res$result)
}

#' Fisher's exact test for score-trait association
#'
#' Performs a Fisher's exact test between each score, binarised at its median into High/Low groups
#' (a score equal to the median is "Low"), and a categorical trait. Significant features are plotted as barplots and
#' saved to the "Results/" folder (\code{Results/Fisher_<trait>_<feature>.pdf}). The statistic reported is the odds
#' ratio for a two-group trait (\code{NA} otherwise).
#'
#' @inheritParams scores.ttest
#' @param trait Character. Name of the column in `coldata` to test (categorical, two or more groups).
#'
#' @return See \code{scores.ttest()}.
#'
#' @examples
#' data("network.tuto")
#' data("deconv_subgroups.tuto")
#' data("traitdata.tuto")
#' cell_groups <- construct_cell_groups(network.tuto, deconv_subgroups.tuto)
#' res <- scores.fisher.test(cell_groups$Cell_groups[, 1:5], traitdata.tuto, "binaryResponse")
#' res$Statistics
#'
#' @export
#'
scores.fisher.test = function(scores, coldata, trait, pval = 0.05, p_adjust = "none"){
  input = prepare_trait_test(scores, coldata, trait, test = "Fisher's exact test")
  X = input$X

  # Each score split at its median into High/Low
  levels_list = lapply(seq_len(ncol(X)), function(j) factor(ifelse(X[, j] > stats::median(X[, j]), "High", "Low"), levels = c("High", "Low")))
  tests = lapply(levels_list, function(l) tryCatch(stats::fisher.test(table(l, input$trait), workspace = 2e7), error = function(e) NULL))
  res = finish_trait_test(input, sapply(tests, function(t) if (is.null(t) || is.null(t$estimate)) NA else unname(t$estimate)),
                          sapply(tests, function(t) if (is.null(t)) NA else t$p.value), "Fisher's exact test", "odds ratio", pval, p_adjust)

  for (j in res$sig) {
    p_value = res$result$Statistics$p_value[j]
    df = data.frame("score_level" = levels_list[[j]], "Trait" = input$trait)
    p_label <- ifelse(p_value < 0.001, "< 0.001", paste0("= ", round(p_value, 4)))
    pdf(paste0("Results/Fisher_", trait, "_", colnames(X)[j], ".pdf"), width = 8, height = 6)
    print(
      ggstatsplot::ggbarstats(df, score_level, Trait, results.subtitle = FALSE,
                              title        = colnames(X)[j],
                              subtitle     = paste0("Fisher's exact test | p ", p_label, " | ", trait, res$adjusted_text[j]),
                              xlab         = trait,
                              legend.title = "Score level") +
        ggplot2::theme(
          plot.title    = ggplot2::element_text(size = 14, face = "bold"),
          plot.subtitle = ggplot2::element_text(size = 10, color = "gray45"),
          axis.text     = ggplot2::element_text(size = 12),
          axis.title    = ggplot2::element_text(size = 13),
          legend.title  = ggplot2::element_text(size = 12),
          legend.text   = ggplot2::element_text(size = 11)
        )
    )
    dev.off()
  }

  return(res$result)

}

#' Perform statistical analysis on scores using a specified test
#'
#' Unified interface for testing associations between any score matrix and a
#' clinical/experimental trait. Accepts cell group scores, NMF latent factors
#' (output of \code{compute.latent_factors()}), or any samples x features matrix.
#'
#' @param scores A list, NMF output from \code{compute.latent_factors()}, or a
#'   score matrix (samples x features). When a list, the first element must be
#'   the score matrix; optional second to fourth elements (e.g. the output of
#'   \code{construct_cell_groups()}) are subsetted to the significant features and returned with their names.
#' @param coldata A data frame (or tibble) containing clinical or experimental metadata for samples.
#'   Must include the column specified in `trait`. Its rows are matched to the samples of the scores by row name;
#'   without sample names, they must be in the order of the scores. Samples with a missing trait are left out.
#' @param trait Character. The name of the column in `coldata` representing the clinical
#'   or experimental trait to test against (e.g., response, subtype, etc.). It must be categorical.
#' @param method Character. Statistical test to perform. One of:
#'   \itemize{
#'     \item `"fisher"` - Fisher's exact test (scores binarised at median)
#'     \item `"wilcox"` - Wilcoxon rank-sum test (non-parametric, binary traits)
#'     \item `"anova"` - One-way ANOVA (parametric, two or more groups)
#'     \item `"kruskal"` - Kruskal-Wallis test (non-parametric, two or more groups)
#'     \item `"ttest"` - Welch's t-test (parametric, binary traits)
#'   }
#'   Only one can be used per call; defaults to the first one (`"fisher"`).
#' @param pval Numeric. Threshold for significance: a feature is significant if its (adjusted) p-value is
#'   \code{<= pval} (default: 0.05).
#' @param p_adjust Character. Method to correct the p-values for multiple testing across the features (see
#'   \code{stats::p.adjust()}, e.g. \code{"BH"}). Default \code{"none"}.
#'
#' @details
#' The function automatically calls the corresponding statistical test function
#' based on the `method` argument:
#' \itemize{
#'   \item \code{scores.fisher.test()}
#'   \item \code{scores.wilcox.test()}
#'   \item \code{scores.anova.test()}
#'   \item \code{scores.kruskal.test()}
#'   \item \code{scores.ttest()}
#' }
#'
#' Each test produces both a statistical result and visual outputs (PDF plots)
#' stored in the `"Results/"` folder for the significant features. These visualizations include the relevant
#' test results (p-values) annotated on the plots.
#'
#' @return A list whose first element is the score matrix restricted to the significant features (no column if none
#'   is significant). Optional second to fourth elements (from the input list, e.g. the composition, weights and
#'   \code{Module_relation} of the cell groups) are subsetted accordingly and keep their names. The last element,
#'   \code{Statistics}, is a data frame with the test of every feature (\code{feature}, \code{test},
#'   \code{statistic_name}, \code{statistic}, \code{p_value}, \code{p_adjusted}, \code{significant}).
#'
#' @examples
#' data("network.tuto")
#' data("deconv_subgroups.tuto")
#' data("traitdata.tuto")
#' cell_groups <- construct_cell_groups(network.tuto, deconv_subgroups.tuto)
#' sig <- scores.stat.analysis(cell_groups, traitdata.tuto, trait = "Enrollment.IC",
#'                             method = "kruskal", pval = 0.05, p_adjust = "BH")
#' head(sig$Statistics)
#' names(sig$Composition)
#'
#' @export
#'
scores.stat.analysis <- function(scores, coldata, trait,
                                 method = c("fisher", "wilcox", "anova", "kruskal", "ttest"),
                                 pval = 0.05, p_adjust = "none") {
  method <- match.arg(method)

  result <- switch(method,
                   fisher   = scores.fisher.test(scores, coldata, trait, pval, p_adjust),
                   wilcox   = scores.wilcox.test(scores, coldata, trait, pval, p_adjust),
                   anova    = scores.anova.test(scores, coldata, trait, pval, p_adjust),
                   kruskal  = scores.kruskal.test(scores, coldata, trait, pval, p_adjust),
                   ttest    = scores.ttest(scores, coldata, trait, pval, p_adjust))

  return(result)
}

#' Compute latent factors from cell group scores using NMF
#'
#' Decomposes signed CCA composite scores into non-negative latent factors
#' by splitting each score into its positive and negative direction components,
#' then applying Non-negative Matrix Factorization via RcppML (fast C++ backend).
#'
#' @param X Numeric matrix of size samples x cell groups (signed CCA composite scores).
#' @param rank Integer; number of NMF factors. If NULL, estimated automatically
#'   at the elbow of the reconstruction MSE across ranks 2:8 (the rank with the largest
#'   second difference of the MSE curve). Ranks are limited to the number of cell groups (no more factors than cell
#'   groups); with fewer than 4 cell groups there are too few candidate ranks for an elbow, and rank = 2 is used
#'   (rank = 1 with a single cell group).
#' @param seed Random seed used for the NMF fits. Default 123. The caller's random number generator
#'   state is restored when the function returns.
#'
#' @param file_name Optional character suffix for the saved patient-mixture plot.
#' @param return Logical. If TRUE (default), saves the patient-mixture barplot to
#'   \code{Results/NMF_patient_mixture_<file_name>.pdf} (\code{Results/NMF_patient_mixture.pdf} if
#'   \code{file_name} is \code{NULL}).
#'
#' @return A named list with:
#' \describe{
#'   \item{Z}{Sample-level NMF factor scores (samples x rank). Non-negative. Each factor (column) sums to 1 over the
#'     samples: the size of each factor is in \code{nmf_model$d} (the NMF is \code{W \%*\% diag(d) \%*\% t(Z)}). Tests
#'     on one factor at a time do not depend on this scale; to compare factors within a sample, multiply each column
#'     by its \code{d}, as done for \code{patient_mixture}.}
#'   \item{W}{Feature weights per factor ((2 x n_CGs) x rank). Non-negative.}
#'   \item{nmf_input}{The positive-negative split matrix fed to NMF (samples x (2 x n_CGs)).}
#'   \item{nmf_model}{The \code{RcppML::nmf()} model object (includes the scaling vector \code{d}).}
#'   \item{patient_mixture}{Long-format data frame of per-sample factor proportions used for the mixture plot: the
#'     factor scores multiplied by the size of each factor (\code{nmf_model$d}), divided by their sum in each sample.}
#' }
#'
#' @details
#' Signed CCA scores are decomposed as:
#' \itemize{
#'   \item \code{score_pos = max(score, 0)}: how far the sample is on the active side of the TF module of that
#'     cell group (aligned with the TF program)
#'   \item \code{score_neg = max(-score, 0)}: how far the sample is on the inactive side (anti-aligned)
#' }
#' Both are concatenated column-wise before NMF. Column names are suffixed
#' with "_pos" and "_neg" to track direction. The scores must not contain missing values.
#'
#' @examples
#' data("network.tuto")
#' data("deconv_subgroups.tuto")
#' cell_groups <- construct_cell_groups(network.tuto, deconv_subgroups.tuto)
#' latent <- compute.latent_factors(cell_groups$Cell_groups, return = FALSE)
#' head(latent$Z)
#' latent$nmf_model$d # size of each factor
#'
compute.latent_factors <- function(X, rank = NULL, seed = 123, file_name = NULL, return = TRUE) {

  if (!requireNamespace("RcppML", quietly = TRUE))
    stop("Install RcppML: install.packages('RcppML')")

  withr::local_preserve_seed()
  set.seed(seed)
  X <- as.matrix(X)
  if (anyNA(X)) stop("X has missing values: NMF needs complete cell group scores")

  # -- Positive-negative split -----------------------------------------------
  X_pos <- pmax(X, 0)
  X_neg <- pmax(-X, 0)
  colnames(X_pos) <- paste0(colnames(X), "_pos")
  colnames(X_neg) <- paste0(colnames(X), "_neg")
  nmf_input <- cbind(X_pos, X_neg)

  A <- t(nmf_input)   # features x samples (RcppML convention)

  # -- Rank estimation -------------------------------------------------------
  if (is.null(rank) && ncol(X) < 4) { # Fewer than 3 candidate ranks (2 to the number of cell groups): no elbow
    rank <- min(2L, ncol(X))
    message("Too few cell groups (", ncol(X), ") to estimate the rank: using rank = ", rank)
  }
  if (is.null(rank)) {
    ranks <- 2:min(8, ncol(X)) # No more factors than cell groups
    message("Estimating optimal rank via elbow on reconstruction MSE (ranks ", min(ranks), ":", max(ranks), ")...")
    .mse <- function(A, fit) {
      recon <- fit$w %*% diag(fit$d, nrow = length(fit$d)) %*% fit$h
      mean((A - recon)^2)
    }
    mse_vals <- vapply(ranks, function(k) {
      fit <- RcppML::nmf(A, k = k, seed = seed, verbose = FALSE)
      .mse(A, fit)
    }, numeric(1))
    names(mse_vals) <- as.character(ranks)
    d2   <- diff(diff(mse_vals))
    rank <- as.integer(names(mse_vals)[which.max(d2) + 1L])
    message("MSE by rank: ", paste(sprintf("k%s=%.3g", names(mse_vals), mse_vals), collapse = "  "))
    message("Selected rank (elbow): ", rank)
  }

  # -- NMF fit ---------------------------------------------------------------
  message("Running RcppML NMF with rank = ", rank, "...")
  nmf_res <- RcppML::nmf(A, k = rank, seed = seed, verbose = FALSE)
  W <- as.matrix(nmf_res$w)    # features x rank
  Z <- t(as.matrix(nmf_res$h)) # samples x rank

  rownames(Z) <- rownames(X)
  colnames(Z) <- paste0("Factor", seq_len(rank))
  colnames(W) <- paste0("Factor", seq_len(rank))
  rownames(W) <- colnames(nmf_input)

  # -- Patient mixture plot ---------------------------------------------------

  # Factor proportions per sample: factor scores with the size of each factor (d), normalised to sum to 1
  Z_sized <- sweep(Z, 2, nmf_res$d, "*")
  Z_norm <- sweep(Z_sized, 1, rowSums(Z_sized), "/")
  Z_norm[is.nan(Z_norm)] <- 0

  factor_levels <- paste0("Factor", seq_len(rank))
  fac_colors    <- setNames(scales::hue_pal()(rank), factor_levels)

  patient_mixture <- as.data.frame(Z_norm) %>%
    dplyr::mutate(sample = rownames(Z_norm)) %>%
    tidyr::pivot_longer(-sample,
                         names_to  = "factor",
                         values_to = "proportion") %>%
    dplyr::mutate(factor = factor(factor, levels = factor_levels))
  # order samples by dominant factor
  sample_order <- patient_mixture %>%
    dplyr::group_by(sample) %>%
    dplyr::slice_max(proportion, n = 1, with_ties = FALSE) %>%
    dplyr::arrange(factor, desc(proportion)) %>%
    dplyr::pull(sample)

  p_mixture <- ggplot2::ggplot(
    patient_mixture %>%
      dplyr::mutate(sample = factor(sample, levels = sample_order)),
    ggplot2::aes(x = sample, y = proportion, fill = factor)
  ) +
    ggplot2::geom_col(width = 1, position = "stack") +
    ggplot2::scale_fill_manual(values = fac_colors, name = "Factor") +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::labs(
      title    = "Patient TME composition - NMF factor mixture",
      subtitle = "Each bar = one patient | Height = fraction of each factor",
      x        = "Patients",
      y        = "Factor proportion"
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      axis.text.x     = ggplot2::element_blank(),
      axis.ticks.x    = ggplot2::element_blank(),
      panel.grid      = ggplot2::element_blank(),
      legend.position = "right"
    )

  # save plot
  if(return){
    if (!dir.exists("Results")) dir.create("Results", recursive = TRUE)
      fname <- if (!is.null(file_name)) {
        paste0("Results/NMF_patient_mixture_", file_name, ".pdf")
      } else {
        "Results/NMF_patient_mixture.pdf"
      }
    
      ggplot2::ggsave(fname, p_mixture, width = 14, height = 4)    
  }

  return(list(
    Z         = Z,
    W         = W,
    nmf_input = nmf_input,
    nmf_model = nmf_res,
    patient_mixture = patient_mixture
  ))
}

#' Extract top-contributing features from NMF latent factors
#'
#' For each NMF factor (column of the W basis matrix), selects the features (cell group parts
#' \code{<group>_pos} and \code{<group>_neg}) whose weight is at or above the specified quantile of the factor's
#' positive weights (zero weights are left out before computing the quantile).
#'
#' @param latent_factors A list returned by \code{compute.latent_factors()},
#'   containing at minimum a matrix \code{W} (features x factors) of NMF basis weights.
#' @param quantile_cutoff Numeric between 0 and 1. Features with weight at or above
#'   this quantile of the factor's positive weights are retained. Default 0.7.
#'
#' @return A named list (one element per factor) of named numeric vectors, where
#'   names are feature names (with their \code{_pos} / \code{_neg} suffix) and values are their NMF basis weights,
#'   sorted descending.
#'
#' @keywords internal
extract_contributing_features <- function(latent_factors,
  quantile_cutoff = 0.7) {

  W <- latent_factors$W

  top_contributors <- list()

  for (factor_name in colnames(W)) {

    w <- W[, factor_name]
    w <- w[!is.na(w) & w > 0]

    if (length(w) == 0) {
      message("No features for ", factor_name)
      next
    }

    # threshold by quantile on all features regardless of _pos/_neg
    thr <- quantile(w, quantile_cutoff, na.rm = TRUE)
    sel <- w[w >= thr]

    # sort by weight descending
    sel <- sort(sel, decreasing = TRUE)

    top_contributors[[factor_name]] <- sel
  }

  return(top_contributors)
}

#' Identify cell-type niches from NMF latent factors
#'
#' For each NMF factor, extracts the top-contributing cell group parts (\code{<group>_pos} / \code{<group>_neg}, by
#' basis weight, see \code{extract_contributing_features()}), maps them to their constituent cell types
#' (\code{compute.composition.matrix()}), and keeps only cell types enriched relative to the background composition.
#' Saves a star-network PDF per factor and returns the weighted cell-type associations.
#'
#' For each factor and cell type:
#' \itemize{
#'   \item the weight is the sum of the NMF weights of the selected parts whose cell group contains the cell type;
#'   \item the enrichment is the fraction of selected parts containing the cell type divided by the fraction of all
#'     cell groups containing it. It is a ratio, not a statistical test: with few selected parts, a cell type
#'     present in only one of them can pass the threshold.
#' }
#' The direction is not taken into account: a selected part can mean more or fewer cells of its group (depending on
#' \code{_pos} / \code{_neg} and on the relation of the group to its TF module, see \code{Module_relation} in
#' \code{construct_cell_groups()}), and both are added together. A cell type of a niche is therefore associated with
#' the factor, without saying whether the factor has more or fewer of it.
#'
#' @param latent_factors A list returned by \code{compute.latent_factors()},
#'   containing \code{W} (features x factors) and the NMF model object.
#' @param dt A named list of deconvolution subgroup results, used to build
#'   the cell-group composition matrix via \code{compute.composition.matrix()}.
#' @param cell.groups A list of cell group definitions as returned by
#'   \code{construct_cell_groups()}.
#' @param enrich_thresh Numeric. Minimum enrichment ratio (see above) for a cell type to be retained per factor.
#'   Default 1.5.
#' @param quantile_cutoff Numeric between 0 and 1. Quantile threshold for
#'   selecting top-contributing cell group parts per factor (see \code{extract_contributing_features()}). Default 0.7.
#' @param cells_extra Optional character vector of additional cell type names, added to the vocabulary used to recognise
#'   cell types in the feature names (see \code{compute.composition.matrix()}).
#' @param return Logical. If \code{TRUE}, saves network PDF plots to
#'   \code{Results/Network_<factor>_<file_name>.pdf} (\code{Results/Network_<factor>.pdf} without \code{file_name}).
#'   Default \code{TRUE}.
#' @param file_name Character. Suffix appended to output file names.
#'
#' @return A named list (one element per factor with at least one enriched cell type) of named numeric vectors giving
#'   the enriched cell types and their weights, sorted descending.
#'
#' @examples
#' data("network.tuto")
#' data("deconv_subgroups.tuto")
#' cell_groups <- construct_cell_groups(network.tuto, deconv_subgroups.tuto)
#' latent <- compute.latent_factors(cell_groups$Cell_groups, return = FALSE)
#' niches <- compute_cells_niches(latent, deconv_subgroups.tuto, cell_groups, return = FALSE)
#' niches
#'
#' @export
compute_cells_niches <- function(latent_factors, dt, cell.groups,
                                 enrich_thresh   = 1.5,
                                 quantile_cutoff = 0.7,
                                 cells_extra     = NULL,
                                 return = TRUE, file_name = NULL) {

  if (!dir.exists("Results")) dir.create("Results", recursive = TRUE)

  # -- extract top contributing cell groups per factor -----------------------
  # returns named list: factor -> named numeric vector
  # names = cell group parts with their _pos/_neg suffix (as the rows of the composition matrix), values = NMF weights (all > 0)
  top_features_list <- extract_contributing_features(
    latent_factors,
    quantile_cutoff = quantile_cutoff
  )

  # -- binary composition matrix ---------------------------------------------
  # rows = cell groups, columns = cell types
  comp_matrix     <- compute.composition.matrix(dt, cell.groups,
                                                cells_extra = cells_extra)
  background_freq <- colMeans(comp_matrix)

  factor_celltype_weights <- list()

  for (factor_name in names(top_features_list)) {

    top_features <- top_features_list[[factor_name]]
    if (length(top_features) == 0) next

    # restrict composition to contributing cell groups
    comp_sub <- comp_matrix[
      rownames(comp_matrix) %in% names(top_features), , drop = FALSE
    ]
    if (nrow(comp_sub) == 0) next

    # -- edge weights: factor -> cell types --------------------------------
    # for each cell type: sum of NMF weights of cell groups containing it
    # all weights are positive - no direction distinction

    # Interpretation of edge
    #   Magnitude = cumulative strength of CGs containing it

    edge_weights <- colSums(
      comp_sub * top_features[rownames(comp_sub)]
    )
    edge_weights <- edge_weights[edge_weights > 0]

    if (length(edge_weights) == 0) next

    # -- enrichment relative to background --------------------------------
    # is this cell type more common in the top cell groups than globally?
    # enrichment >= enrich_thresh -> cell type enriched across cell groups relative to baseline ("Is this cell type more common in my selected subset than expected based on its overall frequency?")

    fg_freq <- colMeans(comp_sub)
    enrich  <- fg_freq / background_freq
    enrich[!is.finite(enrich)] <- NA

    # Interpretation
    # Suppose CD4.cells appears in 70% of the selected top features (fg_freq = 0.7)
    # - Its background frequency = 0.34408602
    # - Enrichment = 0.7 / 0.344 ~= 2.03 -> more than 2x enriched relative to baseline

    keep         <- names(which(enrich >= enrich_thresh))
    edge_weights <- edge_weights[names(edge_weights) %in% keep]

    if (length(edge_weights) == 0) next

    factor_celltype_weights[[factor_name]] <- sort(edge_weights, decreasing = TRUE)

    # -- network plot ------------------------------------------------------
    edges <- data.frame(
      from = factor_name,
      to   = names(edge_weights),
      stringsAsFactors = FALSE
    )

    g <- igraph::graph_from_data_frame(edges, directed = FALSE)

    igraph::V(g)$type <- ifelse(igraph::V(g)$name == factor_name, "factor", "celltype")

    # node size - factor fixed, cell types scaled by edge weight
    ct_names  <- igraph::V(g)$name[igraph::V(g)$type == "celltype"]
    ct_sizes  <- scales::rescale(edge_weights[ct_names], to = c(6, 14))

    igraph::V(g)$size <- ifelse(
      igraph::V(g)$name == factor_name,
      18,
      ct_sizes[igraph::V(g)$name]
    )

    # node color - factor gold, cell types by edge weight magnitude
    max_ew <- max(edge_weights)
    igraph::V(g)$color <- sapply(igraph::V(g)$name, function(v) {
      if (v == factor_name) return("#E69F00")
      # shade from light to dark orange proportional to weight
      w <- edge_weights[v] / max_ew
      colorRampPalette(c("#F5CBA7", "#D55E00"))(100)[round(w * 99) + 1]
    })

    igraph::V(g)$frame.color <- "grey20"
    igraph::V(g)$label.color <- "black"
    igraph::V(g)$label.cex   <- 0.9

    # edge width scaled by weight
    ew_ordered  <- edge_weights[edges$to]
    igraph::E(g)$weight <- ew_ordered
    igraph::E(g)$width  <- scales::rescale(ew_ordered, to = c(1, 5))
    igraph::E(g)$color  <- "#D55E00"

    # Node size is proportional to the absolute value of edge_weights.
    #
    # edge_weights for a cell type represent how strongly the factor is associated with that cell type through the cell groups:
    #
    # Magnitude = cumulative contribution across all cell groups for this factor.

    lay <- igraph::layout_as_star(g,
                           center = which(igraph::V(g)$name == factor_name))

    if(return){
      pdf(paste0("Results/Network_", factor_name, if(!is.null(file_name)) paste0("_", file_name), ".pdf"))
      plot(
        g,
        layout              = lay,
        vertex.label.family = "sans",
        vertex.label.dist   = 2,
        vertex.label.degree = -pi / 2,
        main                = paste("Cell-type niche structure of", factor_name),
        margin              = 0.2
      )
      dev.off()
    }
  }

  # Single list of cell types per factor (no more/fewer split), ordered by weight descending
  return(factor_celltype_weights)
}


#' Project cell group scores onto trained NMF latent factors
#'
#' Applies the positive-negative split to new samples and projects them
#' onto the latent space learned during training using the Moore-Penrose
#' pseudoinverse of the scaled basis matrix W. Negative values of the projection are set to 0 (on the tutorial data
#' the result is the same as an exact non-negative projection). Cell groups missing in \code{scores_test} are
#' tolerated with a warning naming them: the projection uses the other cell groups (on the tutorial data, with 20 of
#' 144 cell groups missing, the factors still correlate 0.999 with the full projection).
#'
#' @param latent_spaces A list returned by \code{compute.latent_factors()},
#'   containing at minimum \code{W} (features x rank) and the RcppML model
#'   object with the \code{nmf_model$d} scaling vector.
#' @param scores_test A samples x cell groups matrix of signed CCA composite
#'   scores for the test cohort. Column names must match those used in training
#'   (before the _pos/_neg suffix was added).
#'
#' @return A non-negative matrix of size samples x rank containing the
#'   projected NMF factor scores for the test cohort.
#'
#' @details
#' The RcppML NMF decomposition is A = W %*% diag(d) %*% H, not A = W %*% H.
#' The diagonal scaling vector d must be absorbed into W before computing
#' the pseudoinverse, otherwise the projected scores are on the wrong scale.
#'
#' @examples
#' data("network.tuto")
#' data("deconv_subgroups.tuto")
#' cell_groups <- construct_cell_groups(network.tuto, deconv_subgroups.tuto)
#' latent <- compute.latent_factors(cell_groups$Cell_groups, return = FALSE)
#' projected <- project_factors(latent, cell_groups$Cell_groups[1:10, ])
#' head(projected)
#'
#' @export
project_factors <- function(latent_spaces, scores_test) {

  scores_test <- as.matrix(scores_test)

  # -- Step 1: apply same pos-neg split as training ---------------------------
  test_pos <- pmax(scores_test, 0)
  test_neg <- pmax(-scores_test, 0)
  colnames(test_pos) <- paste0(colnames(scores_test), "_pos")
  colnames(test_neg) <- paste0(colnames(scores_test), "_neg")
  nmf_input_test <- cbind(test_pos, test_neg)   # samples x (2 * n_CGs)

  # -- Step 2: extract W and absorb d scaling ---------------------------------
  # RcppML stores: A ~= W %*% diag(d) %*% H
  # Must absorb d into W before pseudoinverse: W_scaled = W %*% diag(d)
  W <- latent_spaces$W                          # features x rank
  d <- latent_spaces$nmf_model$d                # length = rank

  W_scaled <- W %*% diag(d, nrow = length(d))  # features x rank, d absorbed

  # -- Step 3: restrict to features shared between training and test ----------
  common_features <- intersect(rownames(W_scaled), colnames(nmf_input_test))

  if (length(common_features) == 0)
    stop("No overlapping features between training NMF input and test scores. ",
         "Check that column names of scores_test match training cell group names.")

  missing_groups <- unique(sub("_(pos|neg)$", "", setdiff(rownames(W_scaled), common_features)))
  if (length(missing_groups) > 0)
    warning(length(missing_groups), " cell group(s) of the training missing in the test scores, the projection uses the others: ",
            paste(utils::head(missing_groups, 10), collapse = ", "), if (length(missing_groups) > 10) ", ...", call. = FALSE)

  W_sub <- W_scaled[common_features, , drop = FALSE]   # common x rank
  Y_sub <- t(nmf_input_test[, common_features,
                              drop = FALSE])            # common x samples

  # -- Step 4: project via pseudoinverse -------------------------------------
  # H_test = pinv(W_sub) %*% Y_sub  ->  Z_test = t(H_test)
  Z_test <- t(MASS::ginv(W_sub) %*% Y_sub)             # samples x rank

  # Clamp small negatives from numerical noise in pseudoinverse
  Z_test <- pmax(Z_test, 0)

  rownames(Z_test) <- rownames(scores_test)
  colnames(Z_test) <- paste0("Factor", seq_len(ncol(Z_test)))

  return(Z_test)
}

#' Project test-set samples onto training NMF factors
#'
#' Convenience wrapper that combines \code{compute.test.set()} and
#' \code{project_factors()} into a single call. Given a trained
#' \code{CellTFusion()} result and a test deconvolution matrix, it
#' (1) reconstructs cell group composite scores for the test samples using
#' the training CCA loadings, and (2) projects those scores onto the trained
#' NMF latent space.
#'
#' @param train_processed A list returned by \code{CellTFusion()}, containing
#'   at minimum \code{Processed_deconvolution}, \code{Cell_groups}, and
#'   \code{Latent_spaces}.
#' @param test_deconv A numeric matrix or data frame of the deconvolution features of the test samples before
#'   subgrouping (samples x features; the same kind of matrix as given to \code{CellTFusion()} as \code{deconv}),
#'   with the column names used during training. If the model was trained with batch correction, the test samples are
#'   centred on their own means (as each training cohort was), so they should come from a single cohort.
#'
#' @return A numeric matrix (test samples x NMF factors) of non-negative
#'   projected factor scores. Cell groups with features missing in \code{test_deconv} are skipped by
#'   \code{compute.test.set()} (with a warning), and the projection uses the other cell groups.
#'
#' @examples
#' \dontrun{
#' res <- CellTFusion(raw.counts_train, deconv = deconv_train)
#' test_factors <- project_test_factors(res, deconv_test)
#' }
#'
#' @export
project_test_factors <- function(train_processed, test_deconv) {

  missing_elements <- setdiff(c("Processed_deconvolution", "Cell_groups", "Latent_spaces"), names(train_processed))
  if (length(missing_elements) > 0)
    stop("train_processed must be the output of CellTFusion(); missing element(s): ", paste(missing_elements, collapse = ", "))

  cell_groups_scores <- compute.test.set(
    deconv_res        = train_processed$Processed_deconvolution,
    cell_groups       = train_processed$Cell_groups,
    features          = names(train_processed$Cell_groups[[2]]),
    deconvolution_test = test_deconv
  )

  project_factors(train_processed$Latent_spaces, cell_groups_scores)
}

#' Prepare CellTFusion cross-validation folds for pipeML
#'
#' Fold construction function for \code{pipeML::compute_features.training.ML()} (argument
#' \code{fold_construction_fun}) that computes \code{CellTFusion} latent factors within each
#' cross-validation fold, without information leakage: cell groups and latent factors are learned with
#' \code{CellTFusion()} on the training samples of each fold only, and the test samples of the fold are
#' projected onto them with \code{project_test_factors()}. It works for classification and survival tasks.
#'
#' @param data Data frame provided by \code{pipeML}: samples as rows and genes as columns (the
#'   \code{features_train} given to \code{compute_features.training.ML()}, i.e. the transposed count matrix),
#'   plus the outcome columns (\code{target} for classification, \code{time} and \code{event} for survival).
#' @param folds Named list with the training rows of each fold, provided by \code{pipeML}.
#' @param bestune Provided by \code{pipeML} when building the features of the final model; \code{NULL} while
#'   running the folds.
#' @param deconv Deconvolution matrix of all the samples (samples x features, e.g. from
#'   \code{multideconv::compute.deconvolution()}), with rows in the same order as \code{data}. Deconvolution is
#'   computed independently for each sample, so computing it once for all samples does not leak information
#'   between folds.
#' @param raw.counts Optional count matrix (genes x samples; TPM values if \code{normalized = FALSE} is passed
#'   to \code{CellTFusion()}) with the samples of \code{data}. Its columns are matched to the samples of
#'   \code{deconv} by name (an error is raised if they are not the same samples); without column names, they must be
#'   in the same order as \code{data}.
#'   \code{pipeML} converts the column names of \code{features_train} with \code{make.names()} (e.g.
#'   \code{"HLA-A"} becomes \code{"HLA.A"}), so gene symbols in \code{data} may no longer match TF regulons and
#'   gene sets; passing \code{raw.counts} keeps the original gene symbols. If \code{NULL} (default),
#'   \code{t(data)} is used.
#' @param batch Optional vector with the cohort of each sample, in the same order as \code{data} (e.g.
#'   \code{coldata$Cohort}). If given, \code{CellTFusion()} runs in multi-cohort mode (see \code{CellTFusion()})
#'   and the test samples of each cohort are projected separately, centered on their own means. Default
#'   \code{NULL} (single cohort).
#' @param ncores Integer. Number of folds computed in parallel. Default 1 (sequential).
#' @param return Logical; in final mode, whether \code{CellTFusion()} saves its outputs and plots in
#'   \code{Results/}. Outputs are never saved for the folds. Default FALSE.
#' @param verbose Logical; whether \code{CellTFusion()} prints its progress messages. Default FALSE.
#' @param ... Other arguments passed to \code{CellTFusion()} (e.g. \code{normalized}, \code{TF.collection},
#'   \code{min_targets_size}, \code{minMod}, \code{corr_mod}, \code{corr}, \code{pval}, \code{cancer_type}).
#'   \code{dt}, \code{tfs} and \code{pathways} are not allowed, as they must be computed within each fold.
#'
#' @return
#' \itemize{
#'   \item In fold mode (\code{bestune = NULL}): each fold is saved to \code{Results/fold_<fold name>.rds}, which
#'     \code{pipeML} reads back, as a list with \code{train_data} (latent factors of the training samples plus
#'     the outcome columns), \code{test_data} (latent factors of the test samples, plus \code{time} and
#'     \code{event} for survival), \code{obs_test} (observed outcome of the test samples), \code{rowIndex} (rows
#'     of the test samples) and \code{fold_name}. The list of folds is returned invisibly. If \code{CellTFusion()}
#'     fails in a fold (e.g. no cell groups are found on its training samples), that fold is skipped with a
#'     warning and is not saved; an error is raised only if all the folds fail.
#'   \item In final mode (\code{bestune} provided): a list with (1) the latent factors of all the samples plus
#'     the outcome columns, (2) the \code{CellTFusion()} result on all the samples, available as
#'     \code{Custom_output} in the \code{pipeML} result and used to compute the features of new samples with
#'     \code{project_test_factors()}, and (3) \code{bestune}.
#' }
#'
#' @details
#' Pass the fixed arguments (\code{deconv}, and optionally \code{raw.counts}, \code{batch}
#' and any \code{CellTFusion()} argument) through \code{fold_construction_args_fixed}. Tunable
#' arguments (\code{fold_construction_args_tunable}) are not supported. With \code{ncores > 1}, the first fold is
#' computed before the others start in parallel, so that the TF and pathway collections cached in
#' \code{Results/} are downloaded only once; the package must be installed for the parallel workers.
#'
#' @examples
#' \dontrun{
#' raw.counts <- CellTFusion::raw.counts.tuto
#' traitdata  <- CellTFusion::traitdata.tuto
#' deconv     <- CellTFusion::deconv.tuto  # deconvolution of all the samples
#'
#' # Classification
#' ml_res <- pipeML::compute_features.training.ML(
#'   features_train = t(raw.counts),
#'   task_type      = "classification",
#'   target_var     = traitdata$Best.Confirmed.Overall.Response,
#'   trait.positive = "PD",
#'   metric         = "AUROC",
#'   k_folds        = 5,
#'   n_rep          = 1,
#'   fold_construction_fun        = prepare_celltfusion_folds,
#'   fold_construction_args_fixed = list(deconv = deconv, raw.counts = raw.counts)
#' )
#'
#' # Features of new samples: projected onto the model trained on all the training samples
#' features_test <- project_test_factors(ml_res$Custom_output, deconv_test)
#'
#' # Survival (time and event of each sample)
#' surv_res <- pipeML::compute_features.training.ML(
#'   features_train = t(raw.counts),
#'   task_type      = "survival",
#'   time_var       = time,
#'   event_var      = event,
#'   k_folds        = 5,
#'   n_rep          = 1,
#'   fold_construction_fun        = prepare_celltfusion_folds,
#'   fold_construction_args_fixed = list(deconv = deconv, raw.counts = raw.counts)
#' )
#' }
#'
#' @export
prepare_celltfusion_folds <- function(data, folds = NULL, bestune = NULL, deconv, raw.counts = NULL,
                                      batch = NULL, ncores = 1,
                                      return = FALSE, verbose = FALSE, ...) {

  if (any(c("dt", "tfs", "pathways") %in% names(list(...)))) {
    stop("'dt', 'tfs' and 'pathways' cannot be given: they must be computed within each fold")
  }

  # Outcome columns given by pipeML in data: "target" (classification) or "time" + "event" (survival)
  if ("target" %in% colnames(data)) {
    outcome_cols <- "target"
  } else if (all(c("time", "event") %in% colnames(data))) {
    outcome_cols <- c("time", "event")
  } else {
    stop("data must contain a 'target' column (classification) or 'time' and 'event' columns (survival)")
  }
  survival <- !identical(outcome_cols, "target")
  outcome <- data[, outcome_cols, drop = FALSE]

  # Counts (genes x samples), deconvolution and cohort of all the samples, in the same order as data
  counts <- if (is.null(raw.counts)) t(data[, setdiff(colnames(data), outcome_cols), drop = FALSE]) else raw.counts
  if (is.data.frame(batch) && ncol(batch) == 1) batch <- batch[[1]]  # one-column data frame or tibble
  if (ncol(counts) != nrow(data) || nrow(deconv) != nrow(data) || (!is.null(batch) && (is.data.frame(batch) || length(batch) != nrow(data)))) {
    stop("raw.counts (columns), deconv (rows) and batch must have the same samples, in the same order, as data")
  }
  if (!is.null(batch)) batch <- as.character(batch)
  if (is.null(rownames(deconv))) rownames(deconv) <- rownames(data)
  if (!is.null(raw.counts) && !is.null(colnames(counts))) {
    # Counts matched to the deconvolution samples by name (not by position)
    idx <- match(make.names(rownames(deconv)), make.names(colnames(counts)))
    if (anyNA(idx)) {
      stop("raw.counts and deconv do not have the same samples: ", sum(is.na(idx)), " samples of deconv are not columns of raw.counts (e.g. ",
           paste(utils::head(rownames(deconv)[is.na(idx)], 3), collapse = ", "), ")")
    }
    counts <- counts[, idx, drop = FALSE]
  }
  colnames(counts) <- rownames(deconv)

  # CellTFusion on a subset of samples: cell groups and latent factors are learned on these samples only
  fit_celltfusion <- function(idx, save) {
    CellTFusion(raw.counts = counts[, idx, drop = FALSE], deconv = deconv[idx, , drop = FALSE],
                batch = if (!is.null(batch)) batch[idx], return = save, verbose = verbose, ...)
  }

  # Latent factors of test samples projected onto a trained CellTFusion result
  # (with batch, each cohort is projected separately, centered on its own means)
  project_samples <- function(structure, idx) {
    cohorts <- if (!is.null(batch)) batch[idx] else rep("all", length(idx))
    projected <- lapply(split(idx, cohorts), function(i) project_test_factors(structure, deconv[i, , drop = FALSE]))
    do.call(rbind, unname(projected))[rownames(deconv)[idx], , drop = FALSE]
  }

  # -----------------------------
  # CASE 1: bestune provided - CellTFusion on all the training samples
  # -----------------------------
  if (!is.null(bestune)) {
    structure <- fit_celltfusion(seq_len(nrow(data)), save = return)
    train_data_final <- cbind(data.frame(structure$Latent_spaces$Z), outcome)
    return(list(train_data_final, structure, bestune))
  }

  # -----------------------------
  # CASE 2: bestune NOT provided - compute folds
  # -----------------------------
  if (is.null(folds)) stop("Provide 'folds' (a list of training row indices) or 'bestune'.")
  fold_names <- if (is.null(names(folds))) paste0("Fold", seq_along(folds)) else names(folds)

  run_fold <- function(i) {
    train_idx <- folds[[i]]
    test_idx  <- setdiff(seq_len(nrow(data)), train_idx)

    # TRAIN: cell groups and latent factors learned on the training samples of the fold
    structure <- fit_celltfusion(train_idx, save = FALSE)
    train_data <- cbind(data.frame(structure$Latent_spaces$Z), outcome[train_idx, , drop = FALSE])

    # TEST: test samples projected onto the latent factors learned on the training samples
    test_data <- data.frame(project_samples(structure, test_idx))
    if (survival) test_data <- cbind(test_data, outcome[test_idx, , drop = FALSE])  # pipeML evaluates each fold with them

    list(
      train_data = train_data,
      test_data  = test_data,
      obs_test   = if (survival) outcome[test_idx, , drop = FALSE] else outcome$target[test_idx],
      rowIndex   = test_idx,
      fold_name  = fold_names[i]
    )
  }

  # A fold that fails (e.g. no cell groups found) returns its error message, so the other folds are not lost
  try_fold <- function(i) tryCatch(run_fold(i), error = function(e) conditionMessage(e))

  if (ncores > 1 && length(folds) > 1) {
    # The first fold runs first so the TF and pathway collections cached in Results/ are downloaded only once
    processed_folds <- list(try_fold(1))
    cl <- parallel::makeCluster(min(ncores, length(folds) - 1))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterCall(cl, setwd, getwd())
    processed_folds <- c(processed_folds, parallel::parLapply(cl, seq_along(folds)[-1], try_fold))
  } else {
    processed_folds <- lapply(seq_along(folds), try_fold)
  }
  names(processed_folds) <- fold_names

  # Failed folds are skipped with a warning (pipeML then evaluates the remaining folds)
  failed <- vapply(processed_folds, is.character, logical(1))
  if (all(failed)) stop("CellTFusion() failed in all the folds. Error in the first fold: ", processed_folds[[1]])
  for (i in which(failed)) warning("Fold '", fold_names[i], "' failed and is skipped: ", processed_folds[[i]], call. = FALSE)
  processed_folds <- processed_folds[!failed]

  # Save each fold (pipeML reads back every Results/fold_*.rds, so drop stale files from earlier runs first)
  dir.create("Results", showWarnings = FALSE)
  file.remove(list.files("Results", pattern = "^fold_.*\\.rds$", full.names = TRUE))
  for (i in seq_along(processed_folds)) {
    saveRDS(processed_folds[[i]], file = file.path("Results", paste0("fold_", names(processed_folds)[i], ".rds")))
  }

  invisible(processed_folds)
}

#' Save a grid of scatter plots for significant module-feature pairs
#'
#' For each pair of columns from \code{matA} and \code{matB} whose p-value in
#' \code{p_mat} is at or below \code{pval}, draws a scatter plot with a linear regression
#' line and annotates it with the correlation coefficient and p-value taken from
#' \code{cor_mat} and \code{p_mat}. The full grid is saved as an SVG file in \code{Results/}.
#' The line is a linear fit of the values given in \code{matA} and \code{matB}.
#' Each panel is 3.2 x 3.2 inches, so the size of the file follows from the number of panels. The column
#' names are the axis titles; names longer than 30 characters are split over several lines at \code{"_"} or \code{"."},
#' and the margins are set from the number of lines.
#'
#' @param matA A numeric matrix or data frame (samples x modules) for the y-axis of each panel.
#'   \code{compute.modules.relationship()} passes the ranks for Spearman correlations, and the values
#'   corrected for batch when \code{batch} is used.
#' @param matB A numeric matrix or data frame (samples x features) for the x-axis of each panel.
#' @param cor_mat A numeric matrix of correlation coefficients between columns of
#'   \code{matA} (rows) and \code{matB} (columns).
#' @param p_mat A numeric matrix of p-values corresponding to \code{cor_mat}.
#' @param file_name Character. Base name for the output SVG file (saved to
#'   \code{Results/<file_name>_scatter_grid.svg}).
#' @param pval Numeric. P-value cutoff for displaying a pair. Default 0.05.
#' @param only_sig Logical. If \code{TRUE} (default), only pairs with
#'   \code{p <= pval} are plotted.
#' @param ncol Integer or \code{NULL}. Number of columns in the plot grid.
#'   If \code{NULL}, set to \code{ceiling(sqrt(n_pairs))}.
#' @param max_pairs Integer. Maximum number of panels: if there are more pairs, only the strongest
#'   correlations are drawn (with a message). Default 100.
#' @param axis_suffix Character added to the axis titles after the column names (e.g. \code{" (rank)"}). Default "".
#'
#' @return Called for its side effect (saves SVG); returns \code{NULL} invisibly.
#'
#' @keywords internal
plot_module_scatter_grid <- function(matA, matB, cor_mat, p_mat,
                                     file_name,
                                     pval = 0.05,
                                     only_sig = TRUE,
                                     ncol = NULL,     # allow NULL for auto
                                     max_pairs = 100,
                                     axis_suffix = "") {

  pairs <- expand.grid(A = colnames(matA),
                       B = colnames(matB),
                       stringsAsFactors = FALSE)

  pairs$keep <- mapply(function(a, b) {
    if (only_sig) isTRUE(p_mat[a, b] <= pval) else TRUE
  }, pairs$A, pairs$B)

  pairs <- pairs[pairs$keep, ]

  # Too many panels cannot be read: keep the strongest correlations
  if (nrow(pairs) > max_pairs) {
    message(nrow(pairs), if (only_sig) " significant", " pairs: only the ", max_pairs, " strongest correlations are drawn in the scatter grid")
    strength <- mapply(function(a, b) abs(cor_mat[a, b]), pairs$A, pairs$B)
    pairs <- pairs[order(-strength)[seq_len(max_pairs)], ]
  }

  if (nrow(pairs) == 0) {
    message("No scatter plots to draw.")
    return(invisible(NULL))
  }

  nplots <- nrow(pairs)

  # ---- automatic ncol if not specified ----
  if (is.null(ncol)) {
    ncol <- ceiling(sqrt(nplots))
  }

  nrow_grid <- ceiling(nplots / ncol)
  svg(paste0("Results/", file_name, "_scatter_grid.svg"),
      width = 3.2 * ncol, height = 3.2 * nrow_grid) # each panel is 3.2 x 3.2 inches
  on.exit(dev.off(), add = TRUE) # the device is closed even if a panel fails

  # Long names are split over several lines (at "_" or ".") so that they fit in the panel
  wrap_label <- function(x, width = 30) {
    parts <- strsplit(x, "(?<=[_.])", perl = TRUE)[[1]]
    lines <- character(0)
    current <- ""
    for (part in parts) {
      if (nchar(current) > 0 && nchar(current) + nchar(part) > width) {
        lines <- c(lines, current)
        current <- part
      } else {
        current <- paste0(current, part)
      }
    }
    c(lines, current)
  }
  x_lines <- max(sapply(unique(pairs$B), function(x) length(wrap_label(paste0(x, axis_suffix)))))
  y_lines <- max(sapply(unique(pairs$A), function(x) length(wrap_label(paste0(x, axis_suffix)))))

  par(
    mfrow = c(nrow_grid, ncol),
    mar = c(2.9 + 1.1 * x_lines, 3.1 + 1.1 * y_lines, 2.5, 1)   # bottom, left, top, right: room for every line of the axis titles
  )
  par(cex = 1) # same text size whatever the number of panels (mfrow reduces it otherwise)

  for (k in seq_len(nplots)) {

    i <- pairs$A[k]
    j <- pairs$B[k]

    df <- data.frame(
      x = matB[, j],
      y = matA[, i]
    )

    r  <- cor_mat[i, j]
    pv <- p_mat[i, j]

    plot(df$x, df$y,
         pch = 16, cex = 1,
         xlab = "",
         ylab = "",
         cex.axis = 0.9)

    # The names of the two columns are the axis titles, centred on the panel so that long names have its whole width
    xlab_lines <- wrap_label(paste0(j, axis_suffix))
    ylab_lines <- wrap_label(paste0(i, axis_suffix))
    mtext(xlab_lines, side = 1, line = 2.6 + 1.1 * (seq_along(xlab_lines) - 1),
          at = graphics::grconvertX(0.5, "nfc", "user"))
    mtext(ylab_lines, side = 2, line = 2.8 + 1.1 * rev(seq_along(ylab_lines) - 1),
          at = graphics::grconvertY(0.5, "nfc", "user"))

    abline(lm(y ~ x, df), col = "steelblue", lwd = 2.5)  # thicker regression line

    mtext(
      paste0("R = ", signif(r, 3),
             " | p = ", format.pval(pv, digits = 2)),
      side = 3,
      line = 0.5,
      cex = 1
    )
  }
}


#' Summarize TF module-trait associations as ANOVA boxplot grids
#'
#' For each categorical trait in \code{coldata}, runs a one-way ANOVA for each
#' TF module, performs Tukey HSD post-hoc tests for significant modules, and
#' saves a multi-panel SVG boxplot grid to \code{Results/}.
#'
#' @param tfs.modules A numeric matrix or data frame of TF module scores
#'   (samples x modules), typically the first element (\code{TFs module matrix}) of the
#'   output from \code{compute.WTCNA()}.
#' @param coldata A data frame of sample metadata. Only character, factor and logical
#'   columns are used as traits. Traits with a single level, or with as many levels as
#'   samples, are skipped with a message; samples with a missing value are left out for that trait.
#' @param pval Numeric. ANOVA (adjusted) p-value threshold for significance (\code{<= pval}). Default 0.05.
#' @param file.name Character. Base name appended to output SVG file names
#'   (\code{Results/ANOVA_boxplot_summary_<file.name>_<trait>.svg}, or \code{Results/ANOVA_boxplot_summary_<trait>.svg}
#'   without \code{file.name}). In \code{<trait>}, characters other than letters, digits, ".", "_", "-" and spaces
#'   (e.g. "/") are replaced by "_".
#' @param ncol Integer. Number of columns in the boxplot facet grid. Default 5. The size of the figure is
#'   computed from the number of panels, and the y-axis of each panel from its data.
#' @param p_adjust Character. Method of \code{stats::p.adjust()} to correct the ANOVA p-values across the modules,
#'   for each trait. Default \code{"none"}.
#'
#' @return Invisibly, a data frame with the ANOVA of every module for every trait tested (\code{trait},
#'   \code{module}, \code{F}, \code{p}, \code{p_adjusted}), or \code{NULL} when no trait could be tested.
#'   Called mainly for its side effect (saves one SVG file per categorical trait with at
#'   least one significant module).
#'
#' @keywords internal
compute.metadata.association.boxplot_summary <- function(
    tfs.modules, coldata, pval = 0.05, file.name,
    ncol = 5, p_adjust = "none"
){

  coldata_cat <- coldata %>%
    dplyr::select(dplyr::where(is.character) | dplyr::where(is.factor) | dplyr::where(is.logical))

  if (ncol(coldata_cat) == 0) {
    message("No categorical traits found.")
    return(invisible(NULL))
  }

  anova_all <- NULL

  # ---- LOOP OVER EACH TRAIT ----
  for(tr in colnames(coldata_cat)){

    # Samples with a value for this trait
    trait_values <- droplevels(as.factor(coldata_cat[[tr]]))
    keep <- !is.na(trait_values)
    n_levels <- nlevels(droplevels(trait_values[keep]))
    if(n_levels < 2 || sum(keep) <= n_levels){
      message("Trait '", tr, "' skipped: ", if(n_levels < 2) "it has a single level" else "it has as many levels as samples")
      next
    }

    df_long <- cbind(as.data.frame(tfs.modules), group = trait_values)[keep, , drop = FALSE] %>%
      as.data.frame() %>%
      tidyr::pivot_longer(
        cols = colnames(tfs.modules),
        names_to = "module",
        values_to = "value"
      ) %>%
      dplyr::mutate(
        trait = tr,
        group = droplevels(group)
      )

    # ---- ANOVA ----
    anova_df <- df_long %>%
      dplyr::group_by(module) %>%
      dplyr::do({
        res <- rstatix::anova_test(data = ., dv = value, between = group)
        tibble::tibble(F = res$F, p = res$p)
      }) %>%
      dplyr::ungroup()
    anova_df$p_adjusted <- stats::p.adjust(anova_df$p, method = p_adjust) # Across the modules, for this trait
    anova_all <- rbind(anova_all, data.frame(trait = tr, as.data.frame(anova_df), stringsAsFactors = FALSE))

    sig_pairs <- anova_df[anova_df$p_adjusted <= pval, ]
    if(nrow(sig_pairs) == 0){
      message("No module significantly associated with trait '", tr, "'")
      next
    }

    feature_labels <- sig_pairs %>%
      dplyr::mutate(
        feature_label = paste0(module, "\nF=", signif(F,3), ", ", if(p_adjust == "none") "p=" else "p.adj=", signif(.data[["p_adjusted"]],3))
      ) %>%
      dplyr::distinct(module, feature_label) %>%
      tibble::deframe()

    df_long <- df_long %>%
      dplyr::semi_join(sig_pairs, by = "module")


    # ---- Tukey ----
    tukey_df <- df_long %>%
      dplyr::group_by(module) %>%
      dplyr::filter(dplyr::n_distinct(group) > 1) %>%
      dplyr::group_modify(~ {
        tuk <- rstatix::tukey_hsd(.x, value ~ group)
        if(nrow(tuk) == 0) return(NULL)
        rstatix::add_xy_position(tuk, x = "group")
      }) %>%
      dplyr::ungroup()

    # ---- PLOT ----
    group_n <- table(trait_values[keep])
    group_labels <- setNames(paste0(names(group_n), "\n(n=", group_n, ")"), names(group_n))

    p <- ggplot2::ggplot(df_long, ggplot2::aes(x = group, y = value, fill = group)) +
      ggplot2::geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85) +
      ggplot2::geom_jitter(width = 0.15, size = 0.8, alpha = 0.5, color = "black") +
      ggplot2::facet_wrap(~ module, ncol = ncol, scales = "free_y",
                 labeller = ggplot2::labeller(module = feature_labels)) +
      ggplot2::scale_x_discrete(labels = group_labels) +
      ggplot2::scale_fill_brewer(palette = "Set2") +
      ggplot2::labs(
        x = NULL,
        y = "Module score",
        title = paste0("TF modules associated with ", tr),
        subtitle = paste0("One-way ANOVA, ", if(p_adjust == "none") "p-value" else paste0(p_adjust, " adjusted p-value"), " <= ", pval,
                          "  |  Brackets: Tukey HSD between groups")
      ) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        strip.text = ggplot2::element_text(size = 11),
        axis.text.x = ggplot2::element_text(size = 10, color = "black"),
        axis.text.y = ggplot2::element_text(size = 10),
        axis.title = ggplot2::element_text(size = 12),
        plot.title = ggplot2::element_text(size = 15, face = "bold"),
        plot.subtitle = ggplot2::element_text(size = 11, color = "gray35"),
        panel.grid.major.x = ggplot2::element_blank(),
        legend.position = "none"
      )

    if(nrow(tukey_df) > 0){
      p <- p + ggpubr::stat_pvalue_manual(
        tukey_df,
        hide.ns = TRUE,
        label = "p.adj.signif",
        size = 4,
        tip.length = 0.01,
        bracket.size = 0.4,
        inherit.aes = FALSE
      )
    }

    # ---- SAVE SVG PER TRAIT (size from the number of panels) ----
    n_panels <- nrow(sig_pairs)
    ggplot2::ggsave(
      paste0("Results/ANOVA_boxplot_summary", if(!is.null(file.name)) paste0("_", file.name), "_",
             gsub("[^A-Za-z0-9._ -]", "_", tr), ".svg"), # Characters not allowed in file names (e.g. "/") replaced by "_"
      plot = p,
      device = grDevices::svg,
      width = max(5, 3.2 * min(ncol, n_panels) + 0.5),
      height = 3.4 * ceiling(n_panels / ncol) + 1.4,
      limitsize = FALSE
    )
  }

  invisible(anova_all)
}

#' Run multivariate feature-based GSEA using limma and Hallmark gene sets
#'
#' This function fits a multivariate linear model for each gene using all features
#' in \code{features_df} as continuous covariates (each feature standardized, and the cohort added when \code{batch}
#' is given). For each feature, it extracts the moderated t-statistics and p-values, ranks genes, and performs GSEA
#' (\code{fgsea}) using the human Hallmark gene sets from MSigDB. Optional dotplots for the top enriched
#' pathways can be saved as PDFs.
#'
#' @param RNA.tpm A numeric matrix or data frame of log-normalized gene expression values (e.g. log2(TPM + 1);
#'   genes as rows, gene symbols as row names, samples as columns).
#' @param features_df A data frame of continuous features (samples in rows,
#'   features in columns) to be modeled as covariates. Rows must be the same samples,
#'   in the same order, as the columns of \code{RNA.tpm} (checked by name; an error is raised otherwise). Features
#'   without variation are left out with a message.
#' @param plot_dot Logical; if TRUE, generates and saves dotplots of top
#'   enriched Hallmark pathways for each feature. Default is TRUE.
#' @param top_n Integer; number of top pathways to display in the dotplot. Default is 10.
#' @param file_name Character; optional suffix for saved PDF files
#'   (\code{Results/GSEA_<feature>_<file_name>.pdf}, or \code{Results/GSEA_<feature>.pdf} without it). Default is NULL.
#' @param width Numeric; width of the PDF plot in inches. Default is 8.
#' @param height Numeric; height of the PDF plot in inches. Default is 10.
#' @param batch Optional vector with the cohort of each sample (in the row order of \code{features_df}, or named by
#'   sample). With two or more cohorts, the cohort is added to the model as a categorical covariate, so that expression
#'   differences between cohorts are not attributed to the features. Default \code{NULL}.
#' @param seed Integer. Random seed used by \code{fgsea::fgsea()} (its p-values come from random permutations), so that
#'   the results are reproducible. Default 123. The caller's random number generator state is not changed.
#'
#' @return A list containing:
#' \describe{
#'   \item{DE_results}{A named list of \code{topTable} results for each feature,
#'     including logFC (change in expression per standard deviation of the feature), moderated t-statistics,
#'     p-values, and adjusted p-values.}
#'   \item{GSEA_results}{A named list of \code{fgsea::fgsea()} result tables
#'     (columns include \code{pathway}, \code{pval}, \code{padj}, \code{NES}) for each feature,
#'     ordered by p-value. fgsea may not be able to estimate the p-value of some pathways: their \code{pval},
#'     \code{padj} and \code{NES} are then \code{NA} (their \code{ES} is kept).}
#' }
#'
#' @details
#' The function works as follows:
#' \enumerate{
#'   \item Hallmark gene sets are retrieved from MSigDB using \code{msigdbr}.
#'   \item A multivariate linear model is fitted for each gene using \code{limma::lmFit}: all the features are in
#'     the same model, so the effect of each feature is adjusted for the others (and for the cohort with \code{batch}).
#'   \item Empirical Bayes moderation is applied via \code{limma::eBayes}.
#'   \item For each feature:
#'     \enumerate{
#'       \item Differential expression results are extracted using \code{topTable} for
#'         the coefficient of that feature.
#'       \item Genes are ranked by moderated t-statistics.
#'       \item Hallmark GSEA is performed on the ranked gene list with \code{fgsea::fgsea()}.
#'       \item Optionally, a dotplot of the top enriched pathways is generated.
#'     }
#' }
#'
#' @examples
#' \dontrun{
#'   results <- compute_factor_gsea(RNA.tpm = expression_matrix,
#'                               features_df = feature_table,
#'                               plot_dot = TRUE,
#'                               top_n = 10,
#'                               file_name = "Feature1_vs_all")
#' }
#'
#' @export
compute_factor_gsea <- function(RNA.tpm,
                             features_df,
                             plot_dot = TRUE,
                             top_n = 10,
                             file_name = NULL,
                             width = 8,
                             height = 10,
                             batch = NULL,
                             seed = 123) {


  features_df = data.frame(features_df)
  if(!identical(make.names(rownames(features_df)), make.names(colnames(RNA.tpm)))){
    stop("Samples must be the same and in the same order in features_df rows and RNA.tpm columns")
  }

  # Features standardized: logFC per standard deviation of the feature (t-statistics and GSEA do not depend on it)
  constant <- vapply(features_df, function(x) stats::sd(x, na.rm = TRUE) == 0 || all(is.na(x)), logical(1))
  if(any(constant)){
    message("Features without variation left out: ", paste(colnames(features_df)[constant], collapse = ", "))
    features_df <- features_df[, !constant, drop = FALSE]
  }
  if(ncol(features_df) == 0) stop("No feature with variation in features_df")
  features_df[] <- lapply(features_df, function(x) as.numeric(scale(x)))

  # Cohort of each sample (by name if named), added to the model as a categorical covariate
  cohort <- NULL
  if(!is.null(batch)){
    if(!is.null(names(batch))){
      idx <- match(make.names(rownames(features_df)), make.names(names(batch)))
      if(anyNA(idx)) stop("batch is named but has no cohort for some samples of features_df")
      batch <- batch[idx]
    } else if(length(batch) != nrow(features_df)){
      stop("batch must have one value per sample of features_df (", nrow(features_df), "), in its row order, or be named by sample")
    }
    if(anyNA(batch)) stop("batch has missing values")
    if(length(unique(batch)) > 1) cohort <- factor(batch)
  }

  # -----------------------------------------
  # Retrieve Hallmark gene sets (argument renamed from category to collection in recent msigdbr versions)
  # -----------------------------------------
  hallmark_df <- if("collection" %in% names(formals(msigdbr::msigdbr))) msigdbr::msigdbr(species = "Homo sapiens", collection = "H") else
    msigdbr::msigdbr(species = "Homo sapiens", category = "H")
  gene_sets <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)

  # -----------------------------------------
  # Fit multivariate limma model (ALL features, and the cohort if given)
  # -----------------------------------------
  design <- if(is.null(cohort)) stats::model.matrix(~ ., data = features_df) else
    stats::model.matrix(~ ., data = cbind(features_df, .cohort = cohort))

  fit <- limma::lmFit(RNA.tpm, design)
  fit <- limma::eBayes(fit)

  deg_list <- list()
  gsea_results_list <- list()

  # -----------------------------------------
  # Loop over features (extract coefficients)
  # -----------------------------------------
  for (feature_name in colnames(features_df)) {

    message("GSEA for ", feature_name)

    # Extract DE results for this feature
    deg <- limma::topTable(fit,
                    coef = feature_name,
                    number = Inf,
                    adjust.method = "BH",
                    sort.by = "P")

    deg_list[[feature_name]] <- deg

    # Rank genes by moderated t-statistic (genes with a name)
    tstats <- fit$t[, feature_name]
    names(tstats) <- rownames(RNA.tpm)
    tstats <- tstats[nchar(names(tstats)) > 0 & !is.na(names(tstats))]
    tstats <- sort(tstats, decreasing = TRUE)

    # Run Hallmark GSEA via fgsea
    gsea_res <- withr::with_seed(seed, fgsea::fgsea( # Fixed seed: reproducible permutation p-values
      pathways    = gene_sets,
      stats       = tstats,
      eps         = 0,
      nPermSimple = 10000
    ))
    gsea_res <- gsea_res[order(gsea_res$pval), ]

    gsea_results_list[[feature_name]] <- gsea_res

    # Optional plotting
    if (plot_dot && !is.null(gsea_res) && nrow(gsea_res) > 0) {

      if (!dir.exists("Results")) dir.create("Results")

      plot_df <- utils::head(gsea_res[order(gsea_res$padj), ], top_n)
      plot_df$pathway <- sub("^HALLMARK_", "", plot_df$pathway)
      plot_df$pathway <- factor(plot_df$pathway, levels = rev(plot_df$pathway))

      p <- ggplot2::ggplot(plot_df,
               ggplot2::aes(x = NES, y = pathway,
                            size = -log10(padj + 1e-10),
                            colour = NES)) +
        ggplot2::geom_point() +
        ggplot2::scale_colour_gradient2(low = "blue", mid = "white", high = "red") +
        ggplot2::labs(title = paste("Hallmark GSEA -", feature_name),
                      x = "NES", y = NULL,
                      size = "-log10(padj)", colour = "NES") +
        ggplot2::theme_bw()

      pdf(paste0("Results/GSEA_", feature_name, if(!is.null(file_name)) paste0("_", file_name), ".pdf"),
          width = width,
          height = height)
      print(p)
      dev.off()
    }
  }

  return(list(
    DE_results = deg_list,
    GSEA_results = gsea_results_list
  ))
}


#' Derive TME meta-programs by clustering Hallmarks across NMF factors
#'
#' Hierarchically clusters (Ward.D2 on Euclidean distance) Hallmark gene sets by
#' their NES profile across NMF factors to identify recurrent transcriptional
#' programs in the TME. Use \code{annotate_metaprograms_TME()} afterwards to label
#' each meta-program with a Bagaev et al. (2021) MFP subtype. The meta-program numbers (\code{MP1}, \code{MP2}, ...)
#' follow the order of the clusters and have no meaning of their own.
#'
#' @param gsea_results Output of \code{compute_factor_gsea()} (a list with a
#'   \code{GSEA_results} element, one table per NMF factor). To derive meta-programs from several cohorts, put the
#'   \code{GSEA_results} lists of the runs together first, with distinct factor names, e.g.
#'   \code{list(GSEA_results = c(setNames(gsea_A$GSEA_results, paste0("A_", names(gsea_A$GSEA_results))),
#'   setNames(gsea_B$GSEA_results, paste0("B_", names(gsea_B$GSEA_results)))))}.
#' @param k Integer. Number of meta-programs (clusters) to extract. If
#'   \code{NULL} (default), chosen at the elbow of the within-cluster sum of squares
#'   (the k with the largest second difference) for k = 2 to min(10, n_hallmarks - 1).
#' @param file_name Optional character suffix for the saved heatmap
#'   (\code{Results/TCGA_meta_programs_<file_name>.pdf}, or \code{Results/TCGA_meta_programs.pdf} without it).
#' @param plot Logical. If \code{TRUE} (default), saves a clustering heatmap of the NES (Hallmarks x factors); NES
#'   beyond -3 and 3 are shown with the extreme colours. With a single factor, the factor columns are not clustered.
#'
#' @return A data frame with one row per meta-program and columns
#'   \code{meta_program} (\code{"MP1"}, \code{"MP2"}, ...) and \code{hallmarks}
#'   (comma-separated Hallmark gene set names).
#'
#' @examples
#' \dontrun{
#' gsea <- compute_factor_gsea(RNA.tpm = counts.norm, features_df = res$Latent_spaces$Z, plot_dot = FALSE)
#' meta_programs <- derive_meta_programs(gsea, file_name = "my_cohort")
#' }
#'
#' @export
derive_meta_programs <- function(gsea_results,
                                 k          = NULL,
                                 file_name = NULL,
                                 plot       = TRUE) {

  ### Build NES matrix from GSEA results
  nes_mat <- build_nes_matrix(gsea_results)

  # -- cluster Hallmarks by NES profile across factors -----------------------
  # distance between Hallmarks = dissimilarity of their NES profiles
  d  <- dist(nes_mat, method = "euclidean")
  hc <- hclust(d, method = "ward.D2")

  # -- determine k -----------------------------------------------------------
  if (is.null(k)) {
    # use elbow on within-cluster sum of squares
    max_k <- min(10, nrow(nes_mat) - 1)
    wss   <- sapply(2:max_k, function(ki) {
      cuts <- cutree(hc, k = ki)
      sum(sapply(unique(cuts), function(cl) {
        members <- nes_mat[cuts == cl, , drop = FALSE]
        if (nrow(members) == 1) return(0)
        sum(apply(members, 2, var)) * (nrow(members) - 1)
      }))
    })
    # second derivative elbow
    d2 <- diff(diff(wss))
    k  <- which.max(d2) + 2L
    message("Selected k = ", k, " meta-programs")
  }

  cluster_assignments <- cutree(hc, k = k)

  # -- build meta-program list -----------------------------------------------
  meta_programs <- lapply(seq_len(k), function(cl) {
    names(cluster_assignments)[cluster_assignments == cl]
  })

  # -- name each meta-program by its dominant Hallmarks ---------------------
  meta_names <- paste0("MP", seq_len(k))
  names(meta_programs) <- meta_names

  # -- heatmap ---------------------------------------------------------------
  if (plot) {

    # annotation: which cluster each Hallmark belongs to
    ann_row <- data.frame(
      Meta_program = meta_names[cluster_assignments],
      row.names    = names(cluster_assignments)
    )

    # color per meta-program
    mp_colors <- setNames(
      scales::hue_pal()(k),
      meta_names
    )

    pdf(paste0("Results/TCGA_meta_programs", if (!is.null(file_name)) paste0("_", file_name), ".pdf"),
        width = 8, height = 10)
    on.exit(dev.off(), add = TRUE) # The file is closed even if the heatmap fails
    pheatmap::pheatmap(
      nes_mat,
      cluster_rows     = hc,
      cluster_cols     = ncol(nes_mat) > 1, # A single factor cannot be clustered
      cutree_rows      = k,
      annotation_row   = ann_row,
      annotation_colors = list(Meta_program = mp_colors),
      color            = colorRampPalette(
                           c("#2471A3", "white", "#C0392B"))(100),
      breaks           = seq(-3, 3, length.out = 101),
      fontsize_row     = 8,
      fontsize_col     = 9,
      main             = paste0("Hallmark meta-programs across NMF factors\n",
                                "(k = ", k, ")")
    )

  }

  df <- data.frame(
    meta_program = names(meta_programs),
    hallmarks    = vapply(meta_programs, paste, character(1L), collapse = ", "),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  return(df)
}

#' Build a Hallmarks x factors NES matrix from GSEA results
#'
#' Combines the per-factor GSEA outputs from \code{compute_factor_gsea()} into
#' a single matrix. All tested Hallmarks are kept regardless of significance;
#' Hallmarks absent from a factor's results, or whose NES could not be computed (fgsea may fail to estimate some
#' pathways, which then have an \code{NA} NES), are filled with 0, i.e. treated as not enriched.
#'
#' @param gsea_results Output list from \code{compute_factor_gsea()}, containing
#'   a \code{GSEA_results} element (named list of \code{fgsea} result tables,
#'   one per factor), or that \code{GSEA_results} list itself. The factors must be named.
#'
#' @return A numeric matrix of NES values with Hallmarks as rows and NMF factors
#'   as columns. Missing Hallmark-factor combinations are set to 0.
#'
#' @examples
#' \dontrun{
#' gsea <- compute_factor_gsea(RNA.tpm = counts.norm, features_df = res$Latent_spaces$Z, plot_dot = FALSE)
#' nes <- build_nes_matrix(gsea)
#' }
#'
#' @export
build_nes_matrix <- function(gsea_results) {

  # Output of compute_factor_gsea(), or its GSEA_results list given directly
  gsea_list <- if (is.list(gsea_results) && !is.null(gsea_results$GSEA_results)) gsea_results$GSEA_results else gsea_results
  if (!is.list(gsea_list) || is.data.frame(gsea_list) || length(gsea_list) == 0) {
    stop("gsea_results must be the output of compute_factor_gsea() (or its GSEA_results list), with at least one factor")
  }
  if (is.null(names(gsea_list)) || any(names(gsea_list) == "") || anyDuplicated(names(gsea_list)) > 0) {
    stop("The factors of gsea_results must have distinct names")
  }

  # get all Hallmarks seen across any factor
  all_hallmarks <- unique(unlist(
    lapply(gsea_list, function(res) {
      if (is.null(res) || nrow(res) == 0) return(character(0))
      res$pathway
    })
  ))

  # build NES matrix: hallmarks x factors
  nes_mat <- sapply(names(gsea_list), function(fac) {

    res <- gsea_list[[fac]]

    # empty factor
    if (is.null(res) || nrow(res) == 0) {
      return(setNames(rep(0, length(all_hallmarks)), all_hallmarks))
    }

    # NES per hallmark - fill missing with 0
    nes <- setNames(res$NES, res$pathway)
    nes[is.na(nes)] <- 0
    scores <- setNames(rep(0, length(all_hallmarks)), all_hallmarks)
    scores[names(nes)] <- nes
    scores
  })

  # nes_mat is now hallmarks x factors
  nes_mat
}

#' Map study factors to TCGA meta-programs
#'
#' For each study NMF factor, scores it against each TCGA meta-program
#' by computing the mean NES of the meta-program's Hallmarks in that factor (NES read with
#' \code{build_nes_matrix()}, as for the references: Hallmarks without a NES count as 0).
#' The meta-program with the highest positive mean NES is the best match (\code{NA} if no
#' meta-program has a positive score). If the
#' meta-program reference has a \code{TME_subtype} column (see
#' \code{annotate_metaprograms_TME()}), it is appended to the output.
#'
#' @param gsea_study Output from compute_factor_gsea() on study cohort.
#' @param cancer_type Character. TCGA cancer type abbreviation identifying which pre-built
#'   TCGA meta-program reference shipped with the package to use. Available: \code{"blca"},
#'   \code{"luad"}, \code{"skcm"}. Ignored if \code{mp_file} is given.
#' @param mp_file Optional. Either a meta-program data frame (as returned by
#'   \code{derive_meta_programs()} / \code{annotate_metaprograms_TME()}, with columns \code{meta_program} and
#'   \code{hallmarks}) or the path to an RData file containing such an object, which must be named
#'   \code{meta_programs}. If NULL, the pre-built reference for \code{cancer_type} is used.
#' @param plot Logical. If TRUE (default), saves a barplot of factor-to-meta-program scores
#'   to \code{Results/Factor_MP_mapping_<file_name>.pdf} (\code{Results/Factor_MP_mapping.pdf} if
#'   \code{file_name} is \code{NULL}).
#' @param file_name Optional character suffix for saving output plots.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{factor_mapping}: Data frame with one row per study factor, in the order of the factors, and columns
#'     \code{factor}, \code{best_MP} (\code{NA} if no meta-program has a positive score), \code{best_score},
#'     \code{all_scores} (all meta-program scores as \code{"MP:score"} pairs) and, if available, \code{TME_subtype}.
#'   \item \code{reference}: The meta-program reference data frame used.
#' }
#'
#' @examples
#' \dontrun{
#' gsea <- compute_factor_gsea(RNA.tpm = counts.norm, features_df = res$Latent_spaces$Z, plot_dot = FALSE)
#' mapping <- map_factors_to_metaprograms(gsea, cancer_type = "blca")
#' mapping$factor_mapping
#' }
#'
#' @export
map_factors_to_metaprograms <- function(gsea_study,
                                cancer_type,
                                mp_file = NULL,
                                plot = TRUE, file_name = NULL) {

  if(is.null(mp_file)){
    cancer_type <- tolower(cancer_type)

    mp_file <- system.file("extdata",
                           paste0("TCGA_meta_programs_", cancer_type, ".RData"),
                           package = "CellTFusion")
    if (mp_file == "") {
      stop("No meta-program file found for cancer type '", cancer_type,
           "'. Available types: blca, luad, skcm.")
    }
  }

  if(is.data.frame(mp_file)){
    meta_programs <- mp_file
  }else{
    reference_env <- new.env()
    load(mp_file, envir = reference_env)
    if(!exists("meta_programs", envir = reference_env, inherits = FALSE)){
      stop("The file ", mp_file, " has no object named 'meta_programs' (objects: ", paste(ls(reference_env), collapse = ", "), ")")
    }
    meta_programs <- reference_env$meta_programs
  }
  if(!is.data.frame(meta_programs) || !all(c("meta_program", "hallmarks") %in% colnames(meta_programs))){
    stop("The meta-program reference must be a data frame with the columns 'meta_program' and 'hallmarks' (see derive_meta_programs())")
  }

  # -- study NES matrix: Hallmarks x study factors, read as for the references
  nes_study <- build_nes_matrix(gsea_study)

  # -- score each study factor against each meta-program ---------------------
  # score = mean NES of the meta-program Hallmarks in this study factor
  results <- lapply(colnames(nes_study), function(fac) {
    message("Processing: ", fac)
    fac_nes <- nes_study[, fac]

    mp_scores <- sapply(meta_programs$meta_program, function(mp) {
      hallmarks <- strsplit(
        meta_programs$hallmarks[meta_programs$meta_program == mp], ", ")[[1]]
      h <- intersect(hallmarks, names(fac_nes))
      if (length(h) == 0) return(NA_real_)
      mean(fac_nes[h], na.rm = TRUE)
    })

    idx        <- which.max(mp_scores)
    best_mp    <- if (length(idx) > 0 && mp_scores[idx] > 0) names(idx)     else NA_character_
    best_score <- if (length(idx) > 0 && mp_scores[idx] > 0) mp_scores[idx] else NA_real_

    data.frame(
      factor           = fac,
      best_MP          = best_mp,
      best_score       = round(best_score, 3),
      all_scores       = paste(
        sprintf("%s:%.3f", names(mp_scores), mp_scores),
        collapse = ", "
      ),
      stringsAsFactors = FALSE
    )
  })

  result_df <- do.call(rbind, results)

  # attach TME subtype from reference if available (rows kept in the order of the factors)
  if ("TME_subtype" %in% colnames(meta_programs)) {
    result_df$TME_subtype <- meta_programs$TME_subtype[match(result_df$best_MP, meta_programs$meta_program)]
  }
  rownames(result_df) <- NULL

  # -- plot: mapping confidence - all MP scores per factor -------------------
  if (plot) {

    if (!dir.exists("Results")) dir.create("Results")

    mfp_colors <- c(
      "IE"              = "#2166AC",
      "IE/F"            = "#92C5DE",
      "F"               = "#F4A582",
      "D"               = "#B2182B",
      "uncharacterized" = "#AAAAAA"
    )

    # parse all_scores back into long format
    score_long <- do.call(rbind, lapply(seq_len(nrow(result_df)), function(i) {
      pairs  <- strsplit(result_df$all_scores[i], ", ")[[1]]
      kv     <- strsplit(pairs, ":")
      data.frame(
        factor = result_df$factor[i],
        MP     = vapply(kv, `[[`, character(1L), 1L),
        score  = as.numeric(vapply(kv, `[[`, character(1L), 2L)),
        stringsAsFactors = FALSE
      )
    }))

    # mark which MP is the winner for each factor
    score_long$is_best <- paste(score_long$factor, score_long$MP) %in%
      paste(result_df$factor, result_df$best_MP)

    # attach TME_subtype label to winner for color
    tme_lookup <- if ("TME_subtype" %in% colnames(result_df)) {
      setNames(result_df$TME_subtype, result_df$factor)
    } else {
      setNames(rep("uncharacterized", nrow(result_df)), result_df$factor)
    }
    score_long$TME_subtype <- ifelse(
      score_long$is_best,
      tme_lookup[score_long$factor],
      NA_character_
    )

    factor_order <- colnames(nes_study)
    score_long$factor <- factor(score_long$factor, levels = factor_order)
    score_long$MP     <- factor(score_long$MP,
                                levels = unique(meta_programs$meta_program))

    p <- ggplot2::ggplot(score_long,
                         ggplot2::aes(x = MP, y = score)) +
      ggplot2::geom_col(fill = "#DDDDDD", width = 0.6) +
      ggplot2::geom_col(
        data = score_long[score_long$is_best, ],
        ggplot2::aes(fill = TME_subtype),
        width = 0.6
      ) +
      ggplot2::scale_fill_manual(values = mfp_colors, name = "TME subtype",
                                 na.translate = FALSE) +
      ggplot2::geom_hline(yintercept = 0, linewidth = 0.3, colour = "black") +
      ggplot2::facet_wrap(~ factor, scales = "free_y") +
      ggplot2::labs(
        title    = "Study factor to meta-program mapping confidence",
        subtitle = "Grey bars = all MP scores. Colored bar = winning MP (color = TME subtype).",
        x        = "Meta-program",
        y        = "Mean NES score"
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(
        axis.text.x     = ggplot2::element_text(angle = 40, hjust = 1, size = 7),
        strip.text      = ggplot2::element_text(size = 8),
        legend.position = "right"
      )

    fname <- paste0("Results/Factor_MP_mapping",
                    if (!is.null(file_name)) paste0("_", file_name) else "",
                    ".pdf")

    n_fac  <- length(unique(score_long$factor))
    n_cols <- min(3L, n_fac)
    n_rows <- ceiling(n_fac / n_cols)

    ggplot2::ggsave(fname, plot = p,
                    width     = n_cols * 3.5,
                    height    = n_rows * 3,
                    device    = grDevices::cairo_pdf,
                    limitsize = FALSE)

    message("Saved factor-to-MP mapping plot to: ", fname)
  }

  list(factor_mapping = result_df, reference = meta_programs)
}

#' Annotate meta-programs with Bagaev TME subtypes
#'
#' Assigns each meta-program a TME subtype (IE, IE/F, F, D) by weighted majority vote
#' across the TCGA factors with a positive score for it: the scores of the factors are
#' summed per MFP subtype and the subtype with the highest total is assigned. Factors without a subtype in
#' \code{factor_tme_df} do not vote (a message lists them).
#'
#' @param meta_programs_df Data frame output of \code{derive_meta_programs()},
#'   with columns \code{meta_program} and \code{hallmarks}.
#' @param factor_tme_df Data frame output of \code{map_factors_to_TME()},
#'   with columns \code{factor} and \code{best_MFP}.
#' @param factors_mp_df The \code{factor_mapping} data frame returned by
#'   \code{map_factors_to_metaprograms()} run on the TCGA factors themselves, with
#'   columns \code{factor} and \code{all_scores}.
#'
#' @return \code{meta_programs_df}, rows in the same order, with an additional (or replaced) \code{TME_subtype}
#'   column (\code{"IE"}, \code{"IE/F"}, \code{"F"}, \code{"D"}, or \code{"uncharacterized"} when no
#'   factor with a subtype has a positive score for the meta-program, the vote is tied, or the winning label
#'   is \code{"uncharacterized"}).
#'
#' @examples
#' \dontrun{
#' # On a TCGA cohort: factors, their GSEA, meta-programs, and the MFP subtype of each factor
#' gsea <- compute_factor_gsea(RNA.tpm = counts.norm, features_df = res$Latent_spaces$Z, plot_dot = FALSE)
#' meta_programs <- derive_meta_programs(gsea)
#' factors_mp <- map_factors_to_metaprograms(gsea, mp_file = meta_programs, plot = FALSE)$factor_mapping
#' tme <- map_factors_to_TME(cancer_name = "BLCA", Z = res$Latent_spaces$Z)
#' meta_programs <- annotate_metaprograms_TME(meta_programs, tme, factors_mp)
#' }
#'
#' @export
annotate_metaprograms_TME <- function(meta_programs_df, factor_tme_df,
                                      factors_mp_df) {

  if (!is.data.frame(meta_programs_df) || !"meta_program" %in% colnames(meta_programs_df)) {
    stop("meta_programs_df must be the output of derive_meta_programs() (column 'meta_program')")
  }
  if (!is.data.frame(factor_tme_df) || !all(c("factor", "best_MFP") %in% colnames(factor_tme_df))) {
    stop("factor_tme_df must be the output of map_factors_to_TME() (columns 'factor' and 'best_MFP')")
  }
  if (!is.data.frame(factors_mp_df) || !all(c("factor", "all_scores") %in% colnames(factors_mp_df))) {
    stop("factors_mp_df must be the factor_mapping of map_factors_to_metaprograms() (columns 'factor' and 'all_scores')")
  }
  labelled <- factor_tme_df$factor[!is.na(factor_tme_df$best_MFP)]
  if (!any(factors_mp_df$factor %in% labelled)) {
    stop("No factor of factors_mp_df has a subtype in factor_tme_df: the factor names of the two tables must match")
  }
  unlabelled <- setdiff(factors_mp_df$factor, labelled)
  if (length(unlabelled) > 0) {
    message("Factors without a TME subtype in factor_tme_df (they do not vote): ", paste(unlabelled, collapse = ", "))
  }

  all_mps <- meta_programs_df$meta_program

  # parse all_scores into a factor x MP score matrix
  score_mat <- do.call(rbind, lapply(seq_len(nrow(factors_mp_df)), function(i) {
    pairs  <- strsplit(factors_mp_df$all_scores[i], ", ")[[1]]
    kv     <- strsplit(pairs, ":")
    scores <- setNames(
      as.numeric(vapply(kv, `[[`, character(1L), 2L)),
      vapply(kv, `[[`, character(1L), 1L)
    )
    scores[all_mps]
  }))
  rownames(score_mat) <- factors_mp_df$factor

  mfp_lookup <- setNames(factor_tme_df$best_MFP, factor_tme_df$factor)

  # for each MP, weighted majority vote across factors with positive score
  assign_mfp <- function(mp) {
    w      <- score_mat[, mp]
    w[is.na(w)] <- 0
    pos    <- w > 0
    if (!any(pos)) return("uncharacterized")
    labels  <- mfp_lookup[rownames(score_mat)[pos]]
    weights <- w[pos][!is.na(labels)] # Factors without a subtype do not vote
    labels  <- labels[!is.na(labels)]
    if (length(labels) == 0) return("uncharacterized")
    agg    <- tapply(weights, labels, sum)
    winner <- names(which(agg == max(agg)))
    if (length(winner) > 1 || winner == "uncharacterized") return("uncharacterized")
    winner
  }

  mp_mfp_df <- data.frame(
    meta_program = all_mps,
    TME_subtype  = vapply(all_mps, assign_mfp, character(1L)),
    stringsAsFactors = FALSE
  )

  # Rows kept in the order of meta_programs_df
  meta_programs_df$TME_subtype <- mp_mfp_df$TME_subtype[match(meta_programs_df$meta_program, mp_mfp_df$meta_program)]
  meta_programs_df

}

# -----------------------------------------------------------------------------
# Bagaev MFP annotation helpers
# -----------------------------------------------------------------------------

#' Annotate NMF factors with Bagaev et al. (2021) MFP subtypes
#'
#' For a single cancer type, matches TCGA patients present in the NMF factor
#' score matrix \code{Z} to the Bagaev et al. (2021) MFP annotation shipped with the
#' package, then tests each factor across the four MFP subtypes (IE, IE/F, F, D)
#' with a Kruskal-Wallis test. Factors with p < 0.05 (fixed threshold, not corrected for the number of factors) are
#' labelled with the subtype that has the highest median factor score; the others are labelled
#' \code{"uncharacterized"}.
#'
#' Samples are matched to patients by the first three parts of their barcode (\code{TCGA-XX-XXXX}; barcodes written
#' with "." instead of "-" are accepted). When a patient has several samples, a tumour sample (sample type code 01 to 09,
#' e.g. \code{-01A}) is kept rather than a normal one (e.g. \code{-11A}), and a message gives how many samples were
#' left out.
#'
#' @param cancer_name Character. Cancer type abbreviation matching the
#'   \code{TCGA_project} column of the internal annotation (case-insensitive,
#'   e.g. \code{"skcm"}).
#' @param Z Numeric matrix. Samples x factors NMF score matrix (row names =
#'   TCGA barcodes).
#' @param plot Logical. If TRUE (default), saves violin/boxplots of factor scores by MFP group
#'   to \code{Results/TME_factors_MFP_<cancer_name>_<file_name>.pdf} (\code{cancer_name} in lower case).
#' @param file_name Optional character suffix for saving output plots.
#'
#' @return A data frame with columns \code{factor}, \code{best_MFP},
#'   \code{kw_pval} (3 significant digits), \code{median_IE}, \code{median_IEF}, \code{median_F},
#'   \code{median_D}, \code{n_samples}, or \code{NULL} (with a warning) if the
#'   cancer type is not annotated or fewer than 10 patients are matched.
#'
#' @examples
#' \dontrun{
#' # Z: latent factors of a TCGA cohort, samples named by their TCGA barcodes
#' tme <- map_factors_to_TME(cancer_name = "BLCA", Z = res$Latent_spaces$Z)
#' }
#'
#' @export
map_factors_to_TME <- function(cancer_name, Z, plot = TRUE, file_name = NULL) {

  annot_file <- system.file("extdata", "annotation.tsv", package = "CellTFusion")
  if (annot_file == "") stop("annotation.tsv not found in inst/extdata.")
  annot <- read.delim(annot_file, row.names = 1, check.names = FALSE,
                      stringsAsFactors = FALSE)
  
  # -- subset annotation to this cancer type --------------------------------
  annot_ct <- annot[toupper(annot$TCGA_project) == toupper(cancer_name), ,
                    drop = FALSE]
  if (nrow(annot_ct) == 0) {
    warning("No annotation rows found for cancer type '", cancer_name, "'.")
    return(NULL)
  }

  # -- harmonize barcodes: keep only first 3 TCGA parts (TCGA-XX-XXXX) -----
  trunc_barcode <- function(x, sep = "-") {
    parts <- strsplit(x, sep, fixed = TRUE)
    vapply(parts, function(p) paste(p[seq_len(min(3L, length(p)))], collapse = "-"),
           character(1L))
  }

  # Z rownames may use "." instead of "-" (R coerces "-" in colnames/rownames)
  z_sep   <- if (any(grepl(".", rownames(Z), fixed = TRUE)) &&
                 !any(grepl("-", rownames(Z), fixed = TRUE))) "." else "-"

  # Several samples of one patient: tumour samples (sample type 01-09) first, so that they are the ones kept
  sample_code <- vapply(strsplit(rownames(Z), z_sep, fixed = TRUE), function(p)
    if (length(p) >= 4) suppressWarnings(as.integer(substr(p[4], 1, 2))) else NA_integer_, integer(1L))
  tumour <- !is.na(sample_code) & sample_code >= 1 & sample_code <= 9
  Z <- Z[order(!tumour), , drop = FALSE]
  rownames(Z) = trunc_barcode(rownames(Z), sep = z_sep)
  n_duplicated <- sum(duplicated(rownames(Z)))
  if (n_duplicated > 0) {
    message(n_duplicated, " samples of patients with several samples were left out (one sample per patient, tumour samples first)")
  }
  rownames(annot_ct) = trunc_barcode(rownames(annot_ct), sep = "-")

  # remove duplicates introduced by truncation - keep first occurrence
  Z        <- Z[!duplicated(rownames(Z)), , drop = FALSE]
  annot_ct <- annot_ct[!duplicated(rownames(annot_ct)), , drop = FALSE]

  common_short <- intersect(rownames(Z), rownames(annot_ct))

  if (length(common_short) < 10) {
    warning("Only ", length(common_short), " matched patients for '",
            cancer_name, "' - skipping Bagaev annotation.")
    return(NULL)
  }

  Z_matched    <- Z[rownames(Z) %in% common_short, , drop = FALSE]
  annot_matched <- annot_ct[rownames(annot_ct) %in% common_short, , drop = FALSE]

  # align rows using truncated barcodes
  Z_matched    <- Z_matched[rownames(annot_matched), , drop = FALSE]

  if(!isTRUE(all.equal(rownames(Z_matched), rownames(annot_matched)))) {
    stop("NMF matrix and TME annotation are not aligned")
  } 

  # -- drop rows with NA/empty MFP ------------------------------------------
  mfp_vals  <- annot_matched$MFP
  Z_matched <- Z_matched[!is.na(mfp_vals), , drop = FALSE]
  mfp_vals  <- factor(mfp_vals[!is.na(mfp_vals)], levels = c("IE", "IE/F", "F", "D"))

  n_samples <- nrow(Z_matched)

  # -- Kruskal-Wallis + median-based labeling --------------------------------
  result_rows <- lapply(colnames(Z_matched), function(fac) {
    fac_scores <- Z_matched[, fac]

    # Kruskal-Wallis across the 4 MFP groups
    kw       <- kruskal.test(fac_scores ~ mfp_vals)
    kw_pval  <- kw$p.value

    # median factor score per group - label = group with highest median
    group_medians <- tapply(fac_scores, mfp_vals, median, na.rm = TRUE)

    best_mfp <- if (kw_pval < 0.05) {
      names(which.max(group_medians))
    } else {
      "uncharacterized"
    }

    data.frame(
      factor           = fac,
      best_MFP         = best_mfp,
      kw_pval          = signif(kw_pval, 3),
      median_IE        = round(group_medians["IE"],   4),
      median_IEF       = round(group_medians["IE/F"], 4),
      median_F         = round(group_medians["F"],    4),
      median_D         = round(group_medians["D"],    4),
      n_samples        = n_samples,
      stringsAsFactors = FALSE,
      row.names        = NULL
    )
  })

  result_df = do.call(rbind, result_rows)

  # -- plot: factor scores by MFP group -------------------------------------
  if (plot) {

    if (!dir.exists("Results")) dir.create("Results")

    # MFP colors consistent with Bagaev palette
    mfp_colors <- c("IE" = "#2166AC", "IE/F" = "#92C5DE", "F" = "#F4A582", "D" = "#B2182B")

    # build long-format data frame for ggplot
    plot_df <- data.frame(
      sample  = rep(rownames(Z_matched), ncol(Z_matched)),
      factor  = rep(colnames(Z_matched), each = nrow(Z_matched)),
      score   = as.vector(Z_matched),
      MFP     = rep(as.character(mfp_vals), ncol(Z_matched)),
      stringsAsFactors = FALSE
    )
    plot_df$MFP    <- factor(plot_df$MFP,    levels = c("IE", "IE/F", "F", "D"))
    plot_df$factor <- factor(plot_df$factor, levels = colnames(Z_matched))

    # label each facet with KW p-value and best_MFP
    facet_labels <- setNames(
      paste0(result_df$factor,
             "\np = ", result_df$kw_pval,
             "  [", result_df$best_MFP, "]"),
      result_df$factor
    )

    p <- ggplot2::ggplot(plot_df,
                         ggplot2::aes(x = MFP, y = score, fill = MFP)) +
      ggplot2::geom_violin(alpha = 0.6, trim = TRUE, scale = "width") +
      ggplot2::geom_boxplot(width = 0.15, outlier.size = 0.5,
                            fill = "white", alpha = 0.8) +
      ggplot2::scale_fill_manual(values = mfp_colors) +
      ggplot2::facet_wrap(~ factor, scales = "free_y",
                          labeller = ggplot2::labeller(factor = facet_labels)) +
      ggplot2::labs(
        title    = paste0("NMF factor scores by Bagaev MFP subtype - ", toupper(cancer_name)),
        subtitle = paste0(n_samples, " matched samples"),
        x        = "MFP subtype",
        y        = "NMF factor score"
      ) +
      ggplot2::theme_bw(base_size = 11) +
      ggplot2::theme(
        legend.position  = "none",
        strip.text       = ggplot2::element_text(size = 8),
        axis.text.x      = ggplot2::element_text(angle = 30, hjust = 1)
      )

    fname <- paste0("Results/TME_factors_MFP_",
                    tolower(cancer_name),
                    if (!is.null(file_name)) paste0("_", file_name) else "",
                    ".pdf")

    # dynamic sizing: ~2 columns, height scales with number of factors
    n_fac  <- ncol(Z_matched)
    n_cols <- 3L
    n_rows <- ceiling(n_fac / n_cols)

    ggplot2::ggsave(fname, plot = p,
                    width  = n_cols * 3,
                    height = n_rows * 3,
                    limitsize = FALSE)

    message("Saved TME factor plot to: ", fname)
  }

  result_df
}

#' Kaplan-Meier survival analysis on clinical groups or CellTFusion features
#'
#' Performs Kaplan-Meier survival analysis and log-rank testing, either on a predefined
#' grouping variable (e.g. a clinical risk group or a supervised cell group split), or
#' automatically for each column of a feature matrix (e.g. latent factor scores from
#' \code{compute.latent_factors()}, TF module scores from \code{compute.WTCNA()}, or
#' cell group scores from \code{construct_cell_groups()}). In the latter case, samples
#' are split into "High"/"Low" groups per feature using a quantile cutoff, and only
#' features with a significant log-rank test are returned.
#'
#' @param survival.data A data frame (or tibble) of clinical/survival metadata, with samples as rows.
#' @param PFS Character. Column name in \code{survival.data} with the survival/follow-up time.
#' @param PFS_event Character. Column name in \code{survival.data} with the event indicator
#'   (1 = event occurred, 0 = censored).
#' @param file_name Optional character. Suffix used when saving Kaplan-Meier plots to \code{Results/}.
#' @param features Optional. A samples x features numeric matrix or data frame (e.g.
#'   \code{latent_spaces$Z}). If provided (and \code{group_column} is
#'   \code{NULL}), each feature is tested individually. Rows must be the same samples, in the same order, as the rows of
#'   \code{survival.data} (checked by name). Mutually exclusive with \code{group_column}.
#' @param p.value Numeric. Log-rank test threshold used to keep a feature as significant when \code{features} is used:
#'   (adjusted) p-value \code{<= p.value}. Default is 0.05.
#' @param thres Numeric between 0 and 1. Quantile cutoff used to split each feature into
#'   High/Low groups when \code{features} is used (High = at or above the cutoff). Default is 0.5 (median split).
#' @param group_column Optional character. Column name in \code{survival.data} defining a
#'   predefined categorical grouping (e.g. cell group membership, cluster, or trait class). Samples without a value
#'   are left out (with a message). Mutually exclusive with \code{features}.
#' @param p_adjust Character. Method to correct the log-rank p-values for multiple testing across the features (see
#'   \code{stats::p.adjust()}, e.g. \code{"BH"}), when \code{features} is used. Default \code{"none"}.
#'
#' @return
#' If \code{group_column} is provided, a list with:
#' \describe{
#'   \item{km_plot}{The \code{survminer::ggsurvplot} object.}
#'   \item{median_PFS}{Median survival time per group (\code{survminer::surv_median()}).}
#'   \item{p_value}{Log-rank test p-value across groups.}
#'   \item{median_follow_up}{Median follow-up time (reverse Kaplan-Meier).}
#' }
#' If \code{features} is provided, a named list with one \code{Surv() ~ feature} formula per significant feature,
#' and a last element \code{Statistics}: a data frame with the test of every feature (\code{feature}, \code{cutoff},
#' \code{n_high}, \code{n_low}, \code{p_value}, \code{p_adjusted}, \code{significant}). When no feature is significant,
#' the list only contains \code{Statistics} (with a message).
#'
#' A Kaplan-Meier plot (with risk table) is saved as an SVG file to
#' \code{Results/SurvPlot_<group-or-feature>_<file_name>.svg} (\code{Results/SurvPlot_<group-or-feature>.svg} without
#' \code{file_name}; characters not allowed in file names, e.g. "/", replaced by "_"): always for
#' \code{group_column}, and one per significant feature for \code{features}.
#'
#' @details
#' Requires the \code{survival}, \code{survminer}, and \code{gridExtra} packages (see \code{Suggests}).
#'
#' @examples
#' \dontrun{
#' # traitdata: clinical data with the survival columns (here "PFS" and "PFS_event");
#' # the tutorial data (traitdata.tuto) does not include them
#'
#' # Predefined clinical/cell groups
#' compute.survival.analysis(
#'   survival.data = traitdata,
#'   PFS           = "PFS",
#'   PFS_event     = "PFS_event",
#'   group_column  = "Best.Confirmed.Overall.Response",
#'   file_name     = "Tutorial"
#' )
#'
#' # Automatic screening of latent factor scores
#' compute.survival.analysis(
#'   survival.data = traitdata,
#'   PFS           = "PFS",
#'   PFS_event     = "PFS_event",
#'   features      = res$Latent_spaces$Z,
#'   p.value       = 0.05,
#'   thres         = 0.5,
#'   file_name     = "Tutorial"
#' )
#' }
compute.survival.analysis = function(survival.data, PFS, PFS_event, file_name = NULL,
                                     features = NULL, p.value = 0.05, thres = 0.5, group_column = NULL,
                                     p_adjust = "none") {

  if (!requireNamespace("survival", quietly = TRUE) ||
      !requireNamespace("survminer", quietly = TRUE) ||
      !requireNamespace("gridExtra", quietly = TRUE)) {
    stop("Packages 'survival', 'survminer', and 'gridExtra' are required for compute.survival.analysis().")
  }

  survival.data <- as.data.frame(survival.data) # Tibbles are accepted
  if (!all(c(PFS, PFS_event) %in% colnames(survival.data))) {
    stop("survival.data must have the columns '", PFS, "' and '", PFS_event, "'")
  }
  file_suffix <- if (!is.null(file_name)) paste0("_", file_name) else ""

  ###### CASE 1: Use predefined groups (like patient clusters or risk groups)
  if (!is.null(group_column)) {

    if (!group_column %in% colnames(survival.data)) stop("group_column '", group_column, "' is not a column of survival.data")
    groups <- survival.data[[group_column]]
    keep <- !is.na(groups) # Samples without a group are left out (not a group of their own)
    if (any(!keep)) message(sum(!keep), " samples without a value for '", group_column, "' are left out")

    data_for_model <- data.frame("time" = survival.data[[PFS]],
                                 "status" = survival.data[[PFS_event]],
                                 "coxHL" = paste0("Group_", groups))[keep, , drop = FALSE]

    # Fit Kaplan-Meier survival curves across groups
    km_fit <- survival::survfit(survival::Surv(time, status) ~ coxHL, data = data_for_model)

    # Median PFS per group
    medians <- survminer::surv_median(km_fit)

    # Median follow-up (reverse KM)
    rev_km <- survival::survfit(survival::Surv(time, 1 - status) ~ 1, data = data_for_model)
    median_follow_up <- summary(rev_km)$table["median"]

    # Log-rank test for survival difference among groups
    pval <- survminer::surv_pvalue(km_fit, data = data_for_model)$pval

    # Legend labels from strata names
    strata_names <- gsub("coxHL=", "", names(km_fit$strata))

    p <- survminer::ggsurvplot(km_fit,
                               data = data_for_model,
                               size = 1,
                               conf.int.style = "step",
                               pval = TRUE,
                               risk.table = TRUE,
                               risk.table.col = "strata",
                               legend.labs = strata_names,
                               risk.table.height = 0.3,
                               ggtheme = ggplot2::theme_grey(),
                               title = paste0("Kaplan-Meier: Surv(time, status) ~ ", group_column),
                               xlab = "Time to death/recurrence/progression")

    p$table <- p$table + ggplot2::theme(legend.position = "none")

    combined <- gridExtra::arrangeGrob(p$plot, p$table, ncol = 1, heights = c(2/3, 1/3))
    ggplot2::ggsave(
      filename = paste0("Results/SurvPlot_", gsub("[^A-Za-z0-9._ -]", "_", group_column), file_suffix, ".svg"), # e.g. "/" replaced by "_"
      plot = combined,
      width = 10,
      height = 5,
      units = "in",
      device = "svg"
    )

    return(list(
      km_plot = p,
      median_PFS = medians,
      p_value = pval,
      median_follow_up = median_follow_up
    ))

  ###### CASE 2: Use features to define high/low groups automatically
  } else if (!is.null(features)) {
    if(!identical(make.names(rownames(features)), make.names(rownames(survival.data)))){
      stop("Samples must be the same and in the same order in features and survival.data rows")
    }
    features <- as.data.frame(features, check.names = FALSE)

    # Log-rank test of each feature split into High/Low at its quantile
    tests <- lapply(colnames(features), function(feature) {
      values <- features[[feature]]
      cutoff <- stats::quantile(values, thres, na.rm = TRUE)
      data_for_model <- data.frame("time" = survival.data[[PFS]],
                                   "status" = survival.data[[PFS_event]],
                                   "coxHL" = factor(ifelse(values >= cutoff, paste0("High_", feature), paste0("Low_", feature))))
      km_fit <- tryCatch(survival::survfit(survival::Surv(time, status) ~ coxHL, data = data_for_model), error = function(e) NULL)
      pval <- if (is.null(km_fit) || length(km_fit$strata) < 2) NA_real_ else survminer::surv_pvalue(km_fit, data = data_for_model)$pval
      list(data = data_for_model, fit = km_fit, cutoff = unname(cutoff), pval = pval,
           n_high = sum(values >= cutoff, na.rm = TRUE), n_low = sum(values < cutoff, na.rm = TRUE))
    })
    p_values <- vapply(tests, function(x) x$pval, numeric(1))
    p_adjusted <- stats::p.adjust(p_values, method = p_adjust)
    sig <- which(!is.na(p_adjusted) & p_adjusted <= p.value)

    statistics <- data.frame(feature = colnames(features), cutoff = vapply(tests, function(x) x$cutoff, numeric(1)),
                             n_high = vapply(tests, function(x) x$n_high, numeric(1)), n_low = vapply(tests, function(x) x$n_low, numeric(1)),
                             p_value = p_values, p_adjusted = p_adjusted, significant = seq_along(p_values) %in% sig,
                             stringsAsFactors = FALSE)

    significant_combinations <- list()
    for (n in sig) {
      feature <- colnames(features)[n]
      km_fit <- tests[[n]]$fit
      data_for_model <- tests[[n]]$data
      significant_combinations[[feature]] <- stats::as.formula(paste0("survival::Surv(time, status) ~ `", feature, "`"))

      strata_names <- gsub("coxHL=", "", names(km_fit$strata))

      p <- survminer::ggsurvplot(km_fit,
                                 data = data_for_model,
                                 size = 1,
                                 palette = c("#E7B800", "#2E9FDF"),
                                 conf.int.style = "step",
                                 pval = TRUE,
                                 risk.table = TRUE,
                                 risk.table.col = "strata",
                                 legend.labs = strata_names,
                                 risk.table.height = 0.3,
                                 ggtheme = ggplot2::theme_grey(),
                                 title = paste0("Kaplan-Meier: Surv(time, status) ~ ", feature),
                                 xlab = "Time to death/recurrence/progression")

      p$table <- p$table + ggplot2::theme(legend.position = "none")
      combined <- gridExtra::arrangeGrob(p$plot, p$table, ncol = 1, heights = c(2/3, 1/3))
      ggplot2::ggsave(
        filename = paste0("Results/SurvPlot_", gsub("[^A-Za-z0-9._ -]", "_", feature), file_suffix, ".svg"), # e.g. "/" replaced by "_"
        plot = combined,
        width = 10,
        height = 5,
        units = "in",
        device = "svg"
      )
    }

    threshold <- paste0(if (p_adjust == "none") "p-value" else paste0("adjusted p-value (", p_adjust, ")"), " <= ", p.value)
    if (length(sig) == 0) {
      message("No feature with a significant log-rank test (", threshold, ")")
    } else {
      message(length(sig), " feature(s) with a significant log-rank test (", threshold, "): ", paste(utils::head(colnames(features)[sig], 10), collapse = ", "))
    }
    significant_combinations[["Statistics"]] <- statistics
    return(significant_combinations)

  } else {
    stop("Please provide either 'group_column' for predefined groups or 'features' for feature-based analysis.")
  }
}
