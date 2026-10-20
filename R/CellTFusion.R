
utils::globalVariables(c("Trait", "Value", ".", "tf", "Module", "values", ".data", "value", "p", "r", "i",
                         "Regulator", "Target", "as.hclust", "condition", "feature_label", "group", "lm", "module",
                         "mor", "mtext", "order.dendrogram", "padj", "pathway", "proportion", "score", "score_level",
                         "setNames", "svg", "target", "n", "n_distinct", "yg", "x", "y", "coefficient", "feature",
                         "weight", "top", "y0",
                         "cell_group", "stable_distinct", "threshold"))

#' @importFrom grDevices colorRampPalette svg
#' @importFrom graphics mtext
#' @importFrom stats as.hclust lm order.dendrogram setNames
#' @importFrom RColorBrewer brewer.pal
#' @importFrom colorspace rainbow_hcl
#' @importFrom ppcor pcor.test
#' @importFrom multideconv compute.deconvolution compute.deconvolution.analysis replicate_deconvolution_subgroups
#' @rawNamespace export(compute.TFs.activity)
#' @rawNamespace export(compute.WTCNA)
#' @rawNamespace export(compute.pathway.activity)
#' @rawNamespace export(compute.metadata.association)
#' @rawNamespace export(compute.latent_factors)
#' @rawNamespace export(compute.modules.relationship)
#' @rawNamespace export(compute.modules.enrichment)
#' @rawNamespace export(compute.test.set)
#' @rawNamespace export(compute.survival.analysis)
#' @rawNamespace export(identify.cell.groups)
NULL


#' Compute one-step CellTFusion
#'
#' Runs the full CellTFusion pipeline in a single call: cell-type deconvolution, TF activity inference,
#' TF module construction, pathway activity, cell group construction, NMF latent factors and cell niches.
#'
#' @param raw_counts A matrix of raw gene expression counts (genes as rows, samples as columns). If
#'   \code{normalized = FALSE}, a matrix of TPM-normalized values instead (not log-transformed).
#'   Always required, also when \code{dt}, \code{tfs} and \code{pathways} are all supplied precomputed (its samples
#'   are checked against them).
#' @param deconv (Optional) A data frame with deconvolution features (cell-type proportions as columns x samples as rows).
#'   If \code{NULL} (default), it is computed with \code{multideconv::compute.deconvolution()} using
#'   \code{deconv_methods}. Ignored if \code{dt} is supplied.
#' @param dt (Optional) A precomputed cell-subgroup object, typically the output of
#'   \code{multideconv::compute.deconvolution.analysis()} (or the \code{Processed_deconvolution} element
#'   returned by a previous \code{CellTFusion()} run). If supplied, cell-type deconvolution and the
#'   deconvolution analysis step are both skipped and the pipeline proceeds straight to cell group
#'   construction using this object. It must have the same samples as \code{raw_counts}.
#' @param tfs (Optional) A precomputed TF activity matrix (samples as rows, TFs as columns), typically
#'   the \code{TFs_matrix} element returned by a previous \code{CellTFusion()} run or by
#'   \code{compute.TFs.activity()}. When \code{batch} is given, a named list of such matrices (one per cohort).
#'   If supplied, TF activity inference is skipped.
#' @param pathways (Optional) A precomputed pathway activity matrix, typically the \code{Pathways_scores}
#'   element returned by a previous \code{CellTFusion()} run or by \code{compute.pathway.activity()}.
#'   If supplied, pathway activity inference is skipped.
#' @param normalized Logical; if TRUE (default), \code{raw_counts} are raw counts and are normalized to
#'   log2(TPM + 1) for the TF activity and pathway activity steps (for deconvolution they are normalized to
#'   TPM only). If FALSE, \code{raw_counts} must already be TPM-normalized and \strong{not} log-transformed: the
#'   TPM values are used as they are for deconvolution and log2(TPM + 1) is computed internally for the other
#'   steps. Do not pass log-transformed values, as they would be log-transformed a second time (a warning is
#'   given when the values do not look like TPM, i.e. no value above 20 or negative values).
#' @param batch (Optional) A vector with the cohort of each sample (e.g. \code{coldata$Cohort}), in the same order
#'   as the columns of \code{raw_counts}; if it is named, it is matched to the samples by name instead. It is
#'   always treated as categorical. If given, the multi-cohort mode is run: TF and pathway activities are computed per
#'   cohort, TF modules are built from a consensus network across cohorts (see \code{consensus_quantile}), and the cohort is controlled for in the
#'   deconvolution analysis and cell group construction. Default is \code{NULL} (single cohort); a vector with
#'   only one cohort is treated the same way.
#' @param deconv_methods A character vector of deconvolution methods to apply. Default is
#'   \code{c("Quantiseq", "CBSX", "Epidish", "DeconRNASeq", "DWLS")}. \code{"CBSX"} (CIBERSORTx) is only run
#'   when \code{cbsx_mail} and \code{cbsx_token} are given; otherwise it is skipped with a message.
#' @param cbsx_mail (Optional) Email credential for CIBERSORTx. Needed to run "CBSX"; without it, CIBERSORTx is skipped.
#' @param cbsx_token (Optional) Token credential for CIBERSORTx. Needed to run "CBSX"; without it, CIBERSORTx is skipped.
#' @param do_parallel Logical; whether to run the deconvolution methods that support it (CIBERSORTx and DWLS) in
#'   parallel, one signature per worker (passed to \code{multideconv::compute.deconvolution()} as \code{doParallel}). Default is FALSE.
#' @param workers (Optional) Number of parallel workers when \code{do_parallel = TRUE}. If \code{NULL} (default),
#'   the number of available cores minus one is used.
#' @param file_name (Optional) Name used in the output files saved in the "Results/" directory. If \code{NULL}
#'   and \code{return = TRUE}, \code{"CellTFusion"} is used.
#' @param TF.collection Character. The source of the TF-target network. Options are `"CollecTRI"` (default), `"Dorothea"`, or `"ARACNE"`.
#' - `"CollecTRI"` uses the prebuilt collection from OmnipathR and `"Dorothea"` the one from the `dorothea` package.
#' - `"ARACNE"` reads a network file (tab-separated with `Regulator` and `Target` columns, and a `mor` column with
#'   the sign of each interaction; see \code{compute.TFs.activity()}) from
#'   \code{input/ARACNE/<cancer_type>/network/network.txt}.
#' @param method Character. The \code{decoupleR} statistic used to score the TFs and the PROGENy pathways: \code{"ulm"}
#'   (univariate linear model, default), \code{"mlm"} (multivariate linear model), \code{"wsum"} (weighted sum) or
#'   \code{"consensus"} (consensus of the three). Passed to \code{compute.TFs.activity()} and \code{compute.pathway.activity()}.
#' @param min_targets_size Integer. Minimum number of target genes per regulon, passed to
#'   \code{compute.TFs.activity()}. Default is 5.
#' @param universe Optional. A user-specified data frame of TF-target interactions (columns \code{source}, \code{target} and
#'   \code{mor}, the mode of regulation: 1 activation, -1 repression). If not provided, the function will fetch the relevant network based on the `TF.collection` argument.
#' @param paths Optional. A user-specified data frame of pathways gene sets. If not provided, the function will fetch the relevant pathways based on `PROGENy`.
#' @param gene_sets Optional. A named list of custom gene sets (character vectors of gene symbols) passed to
#'   \code{compute.pathway.activity()}'s \code{gene_sets} argument for GSVA-based scoring. If \code{NULL}, only PROGENy is used.
#' @param min_module_size Integer; minimum module size for WGCNA module detection. Default is 5.
#' @param corr_mod Numeric; correlation threshold for merging TF modules. Default is 0.9.
#' @param soft_power (Optional) Soft-thresholding power used to build the TF network, passed to
#'   \code{compute.WTCNA()}: one value, or one value per cohort (in the order the cohorts first appear in
#'   \code{batch}) in multi-cohort mode. If \code{NULL} (default), it is chosen automatically.
#' @param consensus_quantile Numeric between 0 and 1, used in multi-cohort mode only: quantile of the topological
#'   overlaps of each pair of TFs across cohorts used as their consensus (0, the default, is the minimum: two TFs are
#'   close only if they are close in every cohort; higher values are less strict). Passed to \code{compute.WTCNA()}.
#' @param clustering_method Linkage method of the hierarchical clusterings (\code{stats::hclust()}): of the TFs, to
#'   build the TF modules (\code{compute.WTCNA()}), and of the deconvolution features of each TF module, to build the
#'   cell groups (\code{construct_cell_groups()}). Default is \code{"ward.D2"}.
#' @param corr Numeric; correlation threshold used in the deconvolution analysis. Default is 0.7.
#' @param corr_type Correlation type used in every step that computes correlations: \code{"spearman"} (default) or
#'   \code{"pearson"}. Passed to \code{multideconv::compute.deconvolution.analysis()} (deconvolution subgroups),
#'   \code{compute.WTCNA()} (TF network and module merging) and \code{construct_cell_groups()} (correlation between TF
#'   modules and deconvolution features).
#' @param cells_extra (Optional) A character vector with the cell names to consider that are not included in the nomenclature of multideconv (see R package).
#'   This includes the group names created with \code{multideconv::aggregate_cell_groups()} on \code{deconv} (e.g. \code{"Lymphocytes"}):
#'   if they are not listed here they are discarded.
#' @param pval Numeric; p-value threshold used when building cell groups (selection of the deconvolution
#'   features significantly correlated with each TF module). Default is 0.05.
#' @param padj Logical; if TRUE, the p-values used when building cell groups are corrected for multiple testing
#'   (Benjamini-Hochberg, across the deconvolution features for each TF module) and compared with \code{pval}.
#'   Passed to \code{construct_cell_groups()}. Default is FALSE.
#' @param min_cluster_size Integer; minimum number of features of a cell group (minimum cluster size of the
#'   dynamic tree cut), passed to \code{construct_cell_groups()}. Default is 3.
#' @param compute_factors Logical. If TRUE (default), the NMF latent factors (\code{compute.latent_factors()}) and the
#'   cell niches (\code{compute_cells_niches()}) are computed from the cell group scores. If FALSE, they are skipped
#'   (\code{Latent_spaces} and \code{Cells_niches} are \code{NULL}), e.g. in the cross-validation folds of
#'   \code{prepare_celltfusion_folds()}, where only the cell group scores are used.
#' @param top_n Integer. Number of top-contributing cell groups per NMF factor in the cell niches (see
#'   \code{compute_cells_niches()}). Default is 5.
#' @param return Logical; if TRUE, intermediate matrices and plots are written to the "Results/" folder. Default is TRUE.
#' @param cancer_type Character. Cancer type abbreviation (e.g. \code{"skcm"}), only used when
#'   \code{TF.collection = "ARACNE"}, to locate the ARACNe network at
#'   \code{input/ARACNE/<cancer_type>/network/network.txt}. If \code{NULL}, the network is auto-detected when only
#'   one exists under \code{input/ARACNE/}.
#' @param verbose Logical; whether to print the progress and the console output and messages of every step
#'   (warnings are always shown). Default is TRUE.
#' @param seed Integer. Random seed of the pipeline (deconvolution methods with random steps and NMF), so that results
#'   are reproducible. Default is 123.
#'
#' @details
#' The pipeline runs with a fixed random seed (\code{seed}) so results are reproducible; the caller's random
#' number generator state is restored when the function returns. With \code{normalized = TRUE}, sample
#' names are converted with \code{make.names()} (e.g. \code{"TCGA-XX-1234"} becomes \code{"TCGA.XX.1234"}),
#' as \code{multideconv::compute.deconvolution()} also does; \code{deconv} rows are matched to the samples
#' in either form.
#'
#' @return A list containing:
#' \describe{
#'   \item{Deconvolution}{A matrix with cell-type proportions (samples as rows, cell types as columns); \code{NULL} if \code{dt} was supplied without \code{deconv}.}
#'   \item{TFs_matrix}{A matrix with TF activity scores (samples as rows, TFs as columns), or a list of matrices (one per cohort) when \code{batch} is given.
#'     Single cohort: raw scores (not scaled), which do not depend on the other samples (e.g. for machine learning features,
#'     to be scaled with the training samples). With \code{batch}: scores z-scored within each cohort (correction of cohort shifts).}
#'   \item{TF_network}{The TF module network returned by \code{compute.WTCNA()}. Its module matrix is always
#'     in the sample order of \code{raw_counts} (also in multi-cohort mode).}
#'   \item{Pathways_scores}{Pathway activity scores returned by \code{compute.pathway.activity()}: raw in single-cohort mode,
#'     z-scored within each cohort with \code{batch} (as \code{TFs_matrix}).}
#'   \item{Processed_deconvolution}{The processed deconvolution subgroups (output of \code{multideconv::compute.deconvolution.analysis()}).}
#'   \item{Cell_groups}{Output of \code{construct_cell_groups()}: cell group scores, compositions and projection weights.}
#'   \item{Latent_spaces}{NMF latent factors returned by \code{compute.latent_factors()} (\code{NULL} with \code{compute_factors = FALSE}).}
#'   \item{Cells_niches}{Top cell groups of each latent factor, with their deconvolution features and CCA coefficients (see \code{compute_cells_niches()}; \code{NULL} with \code{compute_factors = FALSE}).}
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
#'   raw_counts = raw.counts.tuto,
#'   normalized = TRUE,
#'   deconv_methods = c("Quantiseq", "DeconRNASeq"),
#'   file_name = "TestRun",
#'   min_targets_size = 15,
#'   min_module_size = 20,
#'   corr_mod = 0.25,
#'   corr = 0.7,
#'   pval = 0.05
#' )
#'
#' # Re-run with previously computed features, skipping straight to cell group construction
#' res2 <- CellTFusion(
#'   raw_counts = raw.counts.tuto,
#'   dt = res$Processed_deconvolution,
#'   tfs = res$TFs_matrix,
#'   pathways = res$Pathways_scores,
#'   normalized = TRUE,
#'   file_name = "TestRun_rerun",
#'   min_module_size = 20,
#'   corr_mod = 0.25,
#'   pval = 0.05
#' )
#'}
#'
CellTFusion = function(raw_counts, deconv = NULL, dt = NULL, tfs = NULL, pathways = NULL, normalized = TRUE, batch = NULL, deconv_methods = c("Quantiseq", "CBSX", "Epidish", "DeconRNASeq", "DWLS"), cbsx_mail = NULL, cbsx_token = NULL, do_parallel = FALSE, workers = NULL, file_name = NULL,
                       TF.collection = "CollecTRI", method = "ulm", min_targets_size = 5, universe = NULL, paths = NULL, gene_sets = NULL, min_module_size = 5, corr_mod = 0.9, soft_power = NULL, consensus_quantile = 0, clustering_method = "ward.D2", corr = 0.7, corr_type = "spearman", cells_extra = NULL, pval = 0.05, padj = FALSE, min_cluster_size = 3, compute_factors = TRUE, top_n = 5, cancer_type = NULL, return = TRUE, verbose = TRUE, seed = 123){

  withr::local_preserve_seed()
  set.seed(seed)

  # Run a step, silencing its console output and messages when verbose = FALSE
  quiet = function(expr) if(verbose) expr else { utils::capture.output(res <- suppressMessages(expr)); res }

  #Normalize counts
  if(normalized == T){
    counts.norm = quiet(data.frame(ADImpute::NormalizeTPM(raw_counts, log = T)))
  }else{
    if(max(raw_counts, na.rm = T) <= 20 || min(raw_counts, na.rm = T) < 0){ # log2(TPM + 1) cannot exceed 20
      warning("raw_counts does not look like TPM values (range ", signif(min(raw_counts, na.rm = T), 3), " to ", signif(max(raw_counts, na.rm = T), 3),
              "): with normalized = FALSE, TPM-normalized values that are not log-transformed are expected, and they are log-transformed here")
    }
    counts.norm = log2(raw_counts + 1) # raw_counts are TPM: log2(TPM + 1) as above (deconvolution uses the TPM as they are)
  }

  if(return && is.null(file_name)){
    file_name = "CellTFusion"
  }

  #Cohort of each sample (always categorical), in the sample order of the counts
  if(!is.null(batch)){
    if(is.data.frame(batch) && ncol(batch) == 1) batch = batch[[1]] # One-column data frame or tibble
    if(is.data.frame(batch) || length(batch) != ncol(counts.norm)){
      stop("batch must be a vector with the cohort of each sample (", ncol(counts.norm), " samples in raw_counts)")
    }
    if(!is.null(names(batch))){ # Named: matched to the samples by name
      idx = match(make.names(colnames(counts.norm)), make.names(names(batch)))
      if(anyNA(idx)) stop("The names of batch do not match the sample names of raw_counts")
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


  #Deconvolution (skipped entirely if precomputed cell subgroups `dt` are supplied)
  if(is.null(dt)){
    if(is.null(deconv)){
      if(verbose){
        cat("Calculating cell type deconvolution............................................................\n")
      }
      #CIBERSORTx needs credentials: without them it is skipped
      if("CBSX" %in% deconv_methods && (is.null(cbsx_mail) || is.null(cbsx_token))){
        message("CIBERSORTx (CBSX) is skipped: it needs the credentials cbsx_mail and cbsx_token")
        deconv_methods = setdiff(deconv_methods, "CBSX")
        if(length(deconv_methods) == 0) stop("No deconvolution method to run: CBSX needs cbsx_mail and cbsx_token")
      }
      deconv = quiet(compute.deconvolution(raw_counts, normalized = normalized, methods = deconv_methods, credentials.mail = cbsx_mail, credentials.token = cbsx_token, doParallel = do_parallel, workers = workers, file_name = file_name, return = return))
    }

    # Validate sample alignment between counts and deconv
    if (!isTRUE(all.equal(make.names(colnames(counts.norm), unique = TRUE), make.names(rownames(deconv), unique = TRUE)))) {
      stop("Sample mismatch between counts columns and deconvolution rows.")
    }
    rownames(deconv) = colnames(counts.norm) # multideconv converts sample names with make.names() (e.g. "TCGA-XX" -> "TCGA.XX")
  } else { # Precomputed cell subgroups: same samples as raw_counts
    if(!is.list(dt) || is.data.frame(dt) || is.null(dt[["Deconvolution matrix"]])){
      stop("dt must be the output of multideconv::compute.deconvolution.analysis() (element 'Deconvolution matrix')")
    }
    dt_samples = make.names(rownames(dt[["Deconvolution matrix"]]))
    count_samples = make.names(colnames(counts.norm))
    if(!setequal(dt_samples, count_samples)){
      stop("dt and raw_counts must have the same samples: ", sum(!dt_samples %in% count_samples), " samples of dt are not in raw_counts and ",
           sum(!count_samples %in% dt_samples), " samples of raw_counts are not in dt")
    }
  }

  #TF activity (skipped entirely if a precomputed `tfs` matrix is supplied)
  if(is.null(tfs)){
    if (verbose) {
      cat("\nCalculating TF activity............................................................\n")
    }

    if (is.null(batch_vec)) {  ## No batch -> single matrix
      tfs <- quiet(compute.TFs.activity(counts.norm, TF.collection = TF.collection, min_targets_size = min_targets_size, universe = universe, cancer_type = cancer_type, method = method, scale = FALSE, return = return, file_name = file_name))
    }else {
      cohorts <- split(seq_len(ncol(counts.norm)), batch_vec)

      tfs <- lapply(names(cohorts), function(cohort) {
        idx <- cohorts[[cohort]]
        quiet(compute.TFs.activity(counts.norm[, idx, drop = FALSE], TF.collection = TF.collection, min_targets_size = min_targets_size, universe = universe, cancer_type = cancer_type, method = method, scale = TRUE, file_name = paste0(file_name, "_", cohort), return = return))
      })

      names(tfs) <- names(cohorts)
    }
  }

  # 1. TFs network construction
  if(verbose){
    cat("\nConstructing TF network............................................................\n")
  }
  tfs_network = tfs
  if(!is.null(batch_vec)){ # Per-cohort TF activities joined (TFs shared by all cohorts), in the sample order of the counts
    common_tfs = Reduce(intersect, lapply(tfs, colnames))
    if(length(common_tfs) == 0) stop("No common TFs found across cohorts")
    tfs_network = do.call(rbind, unname(lapply(tfs, function(x) x[, common_tfs, drop = FALSE])))
    idx = match(make.names(colnames(counts.norm)), make.names(rownames(tfs_network)))
    if(anyNA(idx)) stop("The TF activities of some samples are missing in tfs")
    tfs_network = tfs_network[idx, , drop = FALSE]
  }
  network = quiet(compute.WTCNA(tfs = tfs_network, batch = batch_vec, network_type = "signed", clustering_method = clustering_method, min_module_size = min_module_size, corr_mod = corr_mod, corr_type = corr_type, soft_power = soft_power, consensus_quantile = consensus_quantile, verbose = verbose, return = return, file_name = file_name))

  # 2. Pathways activity inference (skipped if a precomputed `pathways` matrix is supplied); returned for downstream analyses
  if(is.null(pathways)){
    if(verbose){
      cat("\nCalculating pathway activities............................................................\n")
    }
    pathways = quiet(compute.pathway.activity(counts.norm, gene_sets = gene_sets, paths = paths, method = method, scale = !is.null(batch_vec), batch = batch_vec, return = return, file_name = file_name))
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
  cell.groups = quiet(construct_cell_groups(network, dt, batch = batch_vec, pval = pval, padj = padj, corr_type = corr_type, clustering_method = clustering_method, file_name = file_name, return = return, min_cluster_size = min_cluster_size))

  latent_spaces <- NULL
  cells_niches <- NULL
  if(compute_factors){
    # 5. Cell groups latent spaces
    if(verbose){
      cat("\nLatent spaces calculation............................................................\n")
    }
    latent_spaces <- quiet(compute.latent_factors(cell.groups[[1]], rank = NULL, seed = seed, file_name = file_name, return = return))

    # 6. Extract cell niches
    if(verbose){
      cat("\nExtract cell niches from latent spaces............................................................\n")
    }
    cells_niches = compute_cells_niches(latent_spaces, cell.groups, top_n = top_n, return = return, file_name = file_name)
  }

  if(verbose){
    cat("\nEverything done!", if(return) "Results are saved in Results/ folder", "............................................................\n")
  }

  res = list("Deconvolution" = deconv, "TFs_matrix" = tfs, "TF_network" = network, "Pathways_scores" = pathways,
             "Processed_deconvolution" = dt, "Cell_groups" = cell.groups, "Latent_spaces" = latent_spaces, 
             "Cells_niches" = cells_niches)

  return(res)

}

#' Compute cell group scores from deconvolution and TF module network
#'
#' This function identifies cell groups based on dendrogram cuts, computes composite scores for each group
#' using deconvolution features and TF module networks, and optionally exports the results.
#'
#' @param deconv A data frame with deconvolution features (samples as rows, cell-type features as columns).
#' This is usually the first element returned by \code{multideconv::compute.deconvolution.analysis()}. It must
#' contain all the features of \code{cell_dendrograms}.
#'
#' @param cell_dendrograms A named list of dendrogram (\code{hclust}) objects, each corresponding to a TF module, typically
#' returned by \code{identify.cell.groups()}.
#'
#' @param network A list containing network information of transcription factor (TF) modules, as
#' obtained from \code{compute.WTCNA()} (the module scores, element \code{TFs module matrix}, are used).
#' @param batch Optional vector indicating batch assignment for samples.
#' @param return Logical; if TRUE (default), saves a PDF of the colored cell-group dendrograms
#'   (\code{Results/Dendrogram_color_clusters_<file_name>.pdf}) and CSV files with cell group compositions and scores
#'   to the "Results/" folder (the folder is created if needed).
#' @param file_name Optional character. Suffix of the output file names (see \code{return}).
#' @param min_cluster_size Integer. Minimum number of features of a cell group (minimum cluster size of the
#'   dynamic tree cut of each dendrogram). Default is 3.
#'
#' @return An unnamed list of three elements:
#' \itemize{
#'   \item \code{[[1]]}: A data frame with the composite scores of all identified cell groups across samples. A
#'     positive score means that the sample is on the active side of the TF module of that cell group.
#'   \item \code{[[2]]}: A named list of vectors indicating the composition (original features) of each cell group.
#'   \item \code{[[3]]}: A named list of CCA projection parameters (\code{xcoef}, \code{train_means}, \code{train_sds}; only
#'     \code{xcoef} with \code{batch}) for each cell group.
#' }
#' Cell groups are named \code{Dendrogram_<module>_group_<n>}, where \code{n} is the cluster number of the
#' dendrogram of that module (the number shown in the dendrogram PDF); the three elements use the same names.
#' An error is raised if no cluster is found in any module. Features that the tree cut does not assign to any cluster
#' are not part of any cell group: a message reports how many there are per module.
#' If \code{return=TRUE}, two CSV files will be created:
#' \itemize{
#'   \item \code{Results/Cell.groups.composition_<file_name>.csv}: A table showing the composition of each cell group.
#'   \item \code{Results/Cell.groups.scores_<file_name>.csv}: A matrix of cell group scores across samples.
#' }
#' Without \code{file_name}, they are named \code{Cell.groups.composition.csv} and \code{Cell.groups.scores.csv}.
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
#'   deconv = deconv_results[[1]],
#'   cell_dendrograms = dendrograms,
#'   network = tf_network,
#'   return = TRUE
#' )
#' }
cell.groups.computation = function(deconv, cell_dendrograms, network, batch = NULL, return = TRUE, file_name = NULL, min_cluster_size = 3){

  if(length(cell_dendrograms) == 0){
    stop("No dendrograms given: no TF module has deconvolution features to cluster")
  }
  missing_features = setdiff(unlist(lapply(cell_dendrograms, function(x) as.hclust(x)$labels)), colnames(deconv))
  if(length(missing_features) > 0){
    stop("Features of the dendrograms that are not columns of deconv: ", paste(utils::head(missing_features, 5), collapse = ", "),
         if(length(missing_features) > 5) paste0(", ... (", length(missing_features), " in total)"))
  }
  if(return){
    dir.create("Results", showWarnings = FALSE)
  }

  cell.groups = list()
  unassigned = c()
  cuts_per_dendrogram = calculate_dendrogram_cuts(cell_dendrograms, min_cluster_size = min_cluster_size)

  if(return){
    plot_dendrogram_clusters(cell_dendrograms, cuts_per_dendrogram, file_name = file_name)
  }

  for(k in 1:length(cell_dendrograms)){
    clusters = cuts_per_dendrogram[[k]]
    names(clusters) = as.hclust(cell_dendrograms[[k]])$labels # Same order as the cluster labels (original feature order, not the leaf order of the plot)

    cluster_ids = sort(unique(clusters[clusters > 0]))  # exclude unassigned (label 0)
    unassigned[names(cell_dendrograms)[k]] = sum(clusters == 0)

    y = list()
    x = list()
    z = list()

    for (j in seq_along(cluster_ids)) {
      cells = names(clusters)[clusters == cluster_ids[j]]
      y[[j]] = cells
      ###################################################Compute score for each cell group
      pca_group = deconv[, colnames(deconv) %in% cells, drop = F]
      color = names(cell_dendrograms)[k]
      score = compute_composite_score(pca_group, color, network, batch = batch)
      x[[j]] = score[[1]]
      z[[j]] = score[[2]]
      y[[j]] = rownames(score[[2]]$xcoef) # Features with a weight (the CCA can drop collinear features)
    }

    if(length(x) != 0){
      ###################################################Naming of cell groups: cluster number of the dendrogram (same number as in the dendrogram plot)
      names(x) = names(y) = names(z) = paste0("Dendrogram_", names(cell_dendrograms)[k], "_group_", cluster_ids)
      cell.groups[[k]] = list(x, y, z)
    }else{
      cell.groups[[k]] <- NA #To avoid the problem when putting NULL (list ignore it and not consider it an element)
    }
  }

  names(cell.groups) = names(cell_dendrograms) #Named before removing NA to avoid omitting dendrogram names
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

  cell.groups = list(unlist(unname(cell.groups.values), recursive = F), unlist(unname(cell.groups.names), recursive = F),  unlist(unname(cell.groups.loadings), recursive = F))
  cell.groups[[1]] = cell.groups[[1]] %>%
    data.frame()

  ###################################################Export clusters in a table and save
  clusters = data.frame(matrix(nrow = length(cell.groups[[2]]), ncol = 2))
  for (i in 1:length(cell.groups[[2]])) {
    clusters[i,1] = names(cell.groups[[2]])[[i]]
    clusters[i,2] = paste(cell.groups[[2]][[i]], collapse ="\n")
  }
  colnames(clusters) = c("Cell groups", "Methods-signatures")

  if(return){
    file_suffix = if(!is.null(file_name)) paste0("_", file_name)
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
#' @param scores A numeric matrix or data frame of TF module scores across samples.
#'        Typically the first element (`TFs module matrix`) of the output from `compute.WTCNA()`. Rows represent samples, columns represent TF modules.
#' @param coldata A data frame containing clinical traits (both categorical and numerical) for the same samples.
#'        Rows must be the same samples, in the same order, as the rows of `scores` (checked by name; an
#'        error is raised otherwise). Numeric columns are treated as continuous: convert number-coded categories
#'        (e.g. a stage coded 1 to 4) to factor to test them as groups.
#' @param pval A numeric threshold (default = 0.05) to determine statistical significance (p-value, or adjusted p-value,
#'   \code{<= pval}).
#'        Only associations with (adjusted) p-values below this threshold are annotated in the heatmap or drawn as boxplots.
#' @param corr_type Correlation used for the numeric traits: `"spearman"` (default) or `"pearson"`.
#' @param file_name Optional character suffix of the names of the saved plots (\code{NULL} by default: no suffix).
#' @param ncol_grid Integer. Number of columns in the grid of association boxplots. If \code{NULL} (default), it is set
#'   to the ceiling of the square root of the number of panels.
#' @param plot_grid Logical; if TRUE, tests categorical traits with ANOVA and saves boxplot grids. Default is FALSE.
#' @param padj Logical; if TRUE, the p-values are corrected for multiple testing with the Benjamini-Hochberg method,
#'        across the modules for each trait, and the adjusted p-values are compared with `pval`. Default is FALSE.
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
#'   \item \code{p_adjusted}: matrix of their Benjamini-Hochberg adjusted p-values (equal to \code{p_value} with \code{padj = FALSE}).
#'   \item \code{anova}: data frame with the ANOVA of every module for every categorical trait (\code{trait},
#'     \code{module}, \code{F}, \code{p}, \code{p_adjusted}); \code{NULL} if \code{plot_grid = FALSE}.
#' }
#' The elements are \code{NULL} when there is no trait of that type. The following files are saved to the `Results/` directory:
#' \itemize{
#'   \item \code{TF.modules_metadata_<file_name>.pdf}: a heatmap of module-trait correlations
#'     for numeric traits (only associations with p-value <= \code{pval} are annotated, with the correlation and
#'     its significance: * < 0.05, ** < 0.01, *** < 0.001, **** < 0.0001).
#'   \item \code{ANOVA_boxplot_summary_<file_name>_<trait>.svg} (if \code{plot_grid = TRUE}): one boxplot
#'     grid per categorical trait with the significant modules.
#' }
#'
#' @examples
#'
#' data("network.tuto")
#' data("traitdata.tuto")
#'
#' # Categorical traits (response, Enrollment.IC) are tested with plot_grid = TRUE; os and censOS are numeric
#' res <- compute.metadata.association(
#'   scores = network.tuto[[1]],
#'   coldata = traitdata.tuto,
#'   pval = 0.05,
#'   file_name = 'Tutorial',
#'   plot_grid = TRUE
#' )
#' head(res$anova)
#'
compute.metadata.association <- function(
    scores,
    coldata,
    pval = 0.05,
    corr_type = "spearman",
    file_name = NULL,
    ncol_grid = NULL,
    plot_grid = FALSE,
    padj = FALSE
) {
  if(!corr_type %in% c("spearman", "pearson")){
    stop("corr_type must be 'spearman' or 'pearson'")
  }
  if(!isTRUE(padj) && !isFALSE(padj)){
    stop("padj must be TRUE or FALSE")
  }
  if(!identical(make.names(rownames(scores)), make.names(rownames(coldata)))){
    stop("Samples must be the same and in the same order in scores and coldata rows")
  }

  ### Association with categorical variables (text, factor and TRUE/FALSE columns)
  coldata_categorical <- coldata %>%
    dplyr::select(dplyr::where(is.character) | dplyr::where(is.factor) | dplyr::where(is.logical))

  anova_res <- NULL
  if(ncol(coldata_categorical) != 0) {
    if(plot_grid){
      # Use the new wrapper function to generate boxplot summaries
      anova_res <- compute.metadata.association.boxplot_summary(
        scores = scores,
        coldata = coldata_categorical,
        pval = pval,
        file_name = file_name,
        ncol_grid = ncol_grid,
        padj = padj
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
      scores,
      coldata_quantitative,
      method = corr_type,
      use = "pairwise.complete.obs"
    )
    nObs <- t(!is.na(scores)) %*% (!is.na(coldata_quantitative)) # Samples with both values, per module-trait pair
    moduleTraitPvalue <- WGCNA::corPvalueStudent(moduleTraitCor, nObs)

    # Multiple testing correction across the modules, for each trait
    moduleTraitPadj <- moduleTraitPvalue
    for (i in seq_len(ncol(moduleTraitPadj))) {
      moduleTraitPadj[, i] <- stats::p.adjust(moduleTraitPvalue[, i], method = if(padj) "BH" else "none")
    }

    #### Text of the significant cells: correlation and significance label
    breaks <- c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf)
    labels <- c("****", "***", "**", "*", "")
    stars <- as.character(cut(moduleTraitPadj, breaks = breaks, labels = labels, right = FALSE))
    significant <- !is.na(moduleTraitPadj) & moduleTraitPadj <= pval

    textMatrix <- ifelse(significant, paste0(signif(moduleTraitCor, 2), ifelse(stars != "", paste0("\n", stars), "")), "")
    dim(textMatrix) <- dim(moduleTraitCor)

    p_label <- if(padj) "BH adjusted p-value" else "p-value"
    heatmap_title <- c(paste0("Clinical associations ", file_name),
                       paste0(if(corr_type == "spearman") "Spearman" else "Pearson", " correlation"),
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
      filename = paste0("Results/TF.modules_metadata", if(!is.null(file_name)) paste0("_", file_name), ".pdf")
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
#' @param counts_norm A numeric matrix of normalized gene expression values.
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
#' @param cancer_type Optional character. Cancer type used to locate the ARACNe network
#'        (`input/ARACNE/<cancer_type>/network/network.txt`), only used when `TF.collection = "ARACNE"`. If `NULL`,
#'        the network is auto-detected when only one exists under `input/ARACNE/`.
#' @param top_n Integer. Number of pathways shown in the dot plot of each module. Default is 10. All the
#'        pathways are in the returned tables.
#' @param file_name Optional character suffix for the names of the saved plots.
#'
#' @return Invisibly, a list with two elements, each a named list with one table per module (the result table
#'         of `ReactomePA::enrichPathway()`, or `NULL` for a module without enrichment):
#' \itemize{
#'   \item \code{enrichment}: the module-specific pathways (adjusted p-value < 0.05, and not significant in any other module).
#'   \item \code{enrichment_all}: all the significant pathways of each module, before removing those shared between modules.
#' }
#'         For each module with module-specific pathways, a dot plot is also saved in the `Results/` directory as a PDF file named
#'         `Enrichment_Reactome_Module_<color>_<file_name>.pdf` (`Enrichment_Reactome_Module_<color>.pdf` without `file_name`).
#'         If no enrichment is found for a module, a message is given and no file is saved for that module.
#'         An error is raised if none of the hub TFs is in the TF-target network, or if none of their targets is
#'         found in the row names of `counts_norm` (they must be gene symbols).
#'
#' @details
#' The targets of the hub TFs are taken from the TF-target network (for CollecTRI and Dorothea, cached in
#' `Results/TF_target_collection_<TF.collection>.csv`), and `ReactomePA::enrichPathway()` is used for ORA
#' with the Reactome database. The ORA input is the 20% most variable targets of the module's hub TFs; a module with
#' fewer than 5 such genes is skipped with a message. The background is all the genes of `counts_norm`. Any set of TFs has
#' targets enriched in the pathways those TFs regulate, so the results describe what the hub TFs of each module
#' regulate, not evidence that the module itself is special.
#' Pathways enriched in more than one module are removed, so only module-specific results are kept.
#' TF names are matched in their syntactic form too (e.g. `NKX2.1`, as returned by `compute.TFs.activity()`, matches `NKX2-1`).
#'
#'
#' @examples
#' \dontrun{
#' # Identify the hub TFs of a TF network from compute.WTCNA()
#' hub_tfs <- identify_hub_TFs(wgcna_network, mm_thresh = 0.8, degree_thresh = 0.9)
#'
#' # Perform module-level Reactome enrichment
#' enr <- compute.modules.enrichment(counts.norm, hub_tfs)
#' head(enr$enrichment[[1]])
#' }
compute.modules.enrichment <- function(counts_norm, hub_tfs, TF.collection = "CollecTRI", universe = NULL, cancer_type = NULL, top_n = 10, file_name = NULL){
  if (!requireNamespace("ReactomePA", quietly = TRUE) ||
      !requireNamespace("enrichplot", quietly = TRUE)) {
    stop("Packages 'ReactomePA' and 'enrichplot' are required for compute.modules.enrichment(). Install them with BiocManager::install(c('ReactomePA', 'enrichplot'))")
  }

  if(length(hub_tfs[[1]]) == 0){
    stop("hub_tfs contains no module with hub TFs (see identify_hub_TFs())")
  }

  #TF-target network: the collection used to compute the TF activity (as in compute.TFs.activity())
  if(!is.null(universe)){
    net = universe
  }else if(TF.collection == "ARACNE"){
    if(is.null(cancer_type)){
      # auto-discover when only one network exists
      candidates <- list.files("input/ARACNE", pattern = "^network\\.txt$", recursive = TRUE, full.names = TRUE)
      if(length(candidates) != 1)
        stop("TF.collection = 'ARACNE': specify 'cancer_type' (", length(candidates), " network.txt files found under input/ARACNE/)")
      aracne.network <- candidates[1]
    } else {
      aracne.network <- file.path("input/ARACNE", cancer_type, "network", "network.txt")
      if(!file.exists(aracne.network))
        stop("ARACNe network not found for cancer type '", cancer_type, "': ", aracne.network)
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
  targets_found = sapply(hub_tfs[[1]], function(tfs) sum(rownames(counts_norm) %in% net$target[make.names(net$source) %in% make.names(tfs)]))
  if(all(tfs_found == 0)){
    stop("None of the hub TFs is in the TF-target network: use the TF.collection (or universe) the TF activity was computed with")
  }
  if(all(targets_found == 0)){
    stop("No target gene of the hub TFs was found in the row names of counts_norm: they must be gene symbols")
  }
  no_targets = names(targets_found)[targets_found == 0]
  if(length(no_targets) > 0){
    message("No target gene of the hub TFs found in counts_norm for module(s): ", paste(no_targets, collapse = ", "))
  }

  message("Enrichment computed for the ", length(hub_tfs[[1]]), " module(s) with hub TFs: ", paste(names(hub_tfs[[1]]), collapse = ", "))

  res = list()
  pathways = list()

  #Pathway enrichment using target genes
  for (i in 1:length(hub_tfs[[1]])) {
    color = names(hub_tfs[[1]])[i]
    res[i] = list(module_enrich(as.matrix(counts_norm), color, hub_tfs, net)) #list() keeps modules without enrichment (NULL)
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
        pdf(paste0("Results/Enrichment_Reactome_Module_", color, if(!is.null(file_name)) paste0("_", file_name), ".pdf"))
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
#'   levels, partial correlations (\code{ppcor::pcor.test()}, Pearson or Spearman according to \code{corr_type}) are
#'   computed controlling for batch. Batch is treated as categorical (one dummy variable per level beyond the first).
#' @param pval A numeric threshold for the p-value to define significance (p-value, or adjusted p-value with
#'   \code{padj = TRUE}, <= \code{pval}). Default is 0.05.
#' @param padj Logical; if TRUE, the p-values are corrected for multiple testing with the Benjamini-Hochberg method,
#'   across the columns of \code{matB} for each column of \code{matA}, and the adjusted p-values are used to define
#'   significance (in \code{significant} and in the figures). Default is FALSE.
#' @param corr_type Type of correlation to compute: "spearman" (default) or "pearson". It also sets
#'   the type of the partial correlations computed with \code{batch}.
#' @param vertical Logical; if TRUE, modules are on the x-axis and the \code{matB} features on the y-axis.
#'   Otherwise (default), the \code{matB} features are on the x-axis and modules on the y-axis.
#' @param return Logical; if TRUE, saves the heatmap plot as a PDF (\code{file_name} is needed). Default is TRUE.
#' @param plot_grid Logical; if TRUE, generates per-pair scatter grid plots for significant associations, saved to
#'   "Results/<file_name>_scatter_grid.svg" (the 100 strongest pairs at most). With Spearman, the ranks of the values
#'   are shown, as the correlation is computed on them. With \code{batch}, the values shown are
#'   corrected for batch, as the correlations are. Default is FALSE.
#' @param ncol_grid Integer number of columns used in scatter grid layout. If NULL (default), it is set to
#'   the ceiling of the square root of the number of panels.
#'
#' @return Invisibly, a list with:
#' \itemize{
#'   \item \code{correlation}: the correlation matrix (columns of \code{matA} as rows, columns of \code{matB} as columns).
#'   \item \code{significant}: a named list with, for each column of \code{matA}, the columns of \code{matB} significantly
#'     correlated with it (p-value, or adjusted p-value with \code{padj = TRUE}, <= \code{pval}).
#'   \item \code{p_value}: the matrix of p-values.
#'   \item \code{p_adjusted}: the matrix of Benjamini-Hochberg adjusted p-values (equal to \code{p_value} with
#'     \code{padj = FALSE}).
#' }
#' If `return = TRUE` and `file_name` is given, a heatmap of the correlations is saved to
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
#'                                     return = FALSE,
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
compute.modules.relationship <- function(matA, matB, file_name = NULL, batch = NULL, pval = 0.05, padj = FALSE, corr_type = "spearman", vertical = FALSE, return = TRUE, plot_grid = FALSE, ncol_grid = NULL){

  matA = data.frame(matA, check.names = FALSE)
  matB = data.frame(matB, check.names = FALSE) # keep feature names unchanged (e.g. "-") so they match the deconvolution columns

  if(!corr_type %in% c("spearman", "pearson")){
    stop("corr_type must be 'spearman' or 'pearson'")
  }
  if(!isTRUE(padj) && !isFALSE(padj)){
    stop("padj must be TRUE or FALSE")
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
  if(corr_type == "spearman"){
    matA_plot[] = lapply(matA, rank, na.last = "keep")
    matB_plot[] = lapply(matB, rank, na.last = "keep")
  }
  cor_method = corr_type

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
    cor_label = paste0(if(corr_type == "spearman") "Spearman" else "Pearson", " partial correlation (corrected for batch)")
  } else {
    # ---------- STANDARD CORRELATION (samples with both values, for each pair of columns) ----------
    moduleTraitCor = WGCNA::cor(matA, matB, method = corr_type, use = "pairwise.complete.obs")
    nObs = t(!is.na(matA)) %*% (!is.na(matB))
    moduleTraitPvalue = WGCNA::corPvalueStudent(moduleTraitCor, nObs)
    dimnames(moduleTraitPvalue) = dimnames(moduleTraitCor)
    cor_label = paste0(if(corr_type == "spearman") "Spearman" else "Pearson", " correlation")
  }

  # Multiple testing correction (Benjamini-Hochberg) across the features of matB, for each column of matA
  moduleTraitPadj = moduleTraitPvalue
  if(padj){
    for(i in seq_len(nrow(moduleTraitPadj))){
      moduleTraitPadj[i, ] = stats::p.adjust(moduleTraitPvalue[i, ], method = "BH")
    }
  }
  p_label = if(padj) "adjusted p-value" else "p-value"

  #Significant features of matB for each column of matA
  sig = lapply(seq_len(nrow(moduleTraitPadj)), function(i) colnames(moduleTraitPadj)[which(moduleTraitPadj[i,] <= pval)])
  names(sig) = rownames(moduleTraitPadj)

  if((return || plot_grid) && is.null(file_name)){
    message("No file_name given: the figures are not saved")
  }

  if (plot_grid && !is.null(file_name)) {
    plot_module_scatter_grid(
      matA = matA_plot,
      matB = matB_plot,
      cor_mat = moduleTraitCor,
      p_mat = moduleTraitPadj,
      file_name = file_name,
      pval = pval,
      ncol = ncol_grid,
      only_sig = TRUE,
      axis_suffix = if(corr_type == "spearman") " (rank)" else "",
      p_label = if(padj) "padj" else "p"
    )
  }

  ###Heatmap of the correlations (size computed from the number of features)
  if(return && !is.null(file_name)){
    significant = !is.na(moduleTraitPadj) & moduleTraitPadj <= pval
    stars = as.character(cut(moduleTraitPadj, breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf), labels = c("****", "***", "**", "*", ""), right = FALSE))
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

    heatmap_title = c(cor_label, paste0("shown for ", p_label, " <= ", pval), if(show_numbers) "* < 0.05, ** < 0.01, *** < 0.001, **** < 0.0001")
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

  invisible(list(correlation = moduleTraitCor, significant = sig, p_value = moduleTraitPvalue, p_adjusted = moduleTraitPadj))
}

#' Compute pathway activity scores
#'
#' This function computes pathway activity scores from log-normalized gene expression data with \code{decoupleR}
#' (Badia-i-Mompel et al., 2022), based on the PROGENy resource (Schubert et al., 2018). By default the univariate
#' linear model (\code{"ulm"}) is used.
#' Optionally, it also performs Gene Set Variation Analysis (GSVA) using hallmark signatures or any user-provided gene sets.
#'
#' @param counts_norm A numeric matrix of log-normalized gene expression values (e.g. log2(TPM + 1)) with genes (symbols) as
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
#' @param batch Optional vector with the cohort of each sample, in the column order of \code{counts_norm} (if it is named,
#'   it is matched to the samples by name instead). It is always treated as categorical. If given, the samples are split
#'   by cohort, the pathway activities (PROGENy and, if requested, GSVA) are computed and scaled (see \code{scale})
#'   separately in each cohort, and the results are joined in the sample order of \code{counts_norm}; the output has the same
#'   form as without \code{batch}. Pathways or gene sets that cannot be scored in every cohort are left out (with a message).
#'   Default is \code{NULL} (single cohort); a vector with only one cohort is treated the same way.
#' @param return Logical; if TRUE, saves the score matrices in Results/ folder. Default is TRUE.
#' @param file_name Optional character suffix used when writing output CSV files
#'   (\code{Results/Pathway_matrix_PROGENy_<file_name>.csv} and \code{Results/Pathway_matrix_GSVA_<file_name>.csv};
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
compute.pathway.activity <- function(counts_norm, gene_sets = NULL, paths = NULL, method = "ulm", scale = TRUE, batch = NULL, return = TRUE, file_name = NULL) {

  if (!method %in% c("ulm", "mlm", "wsum", "consensus")) {
    stop("method must be 'ulm', 'mlm', 'wsum' or 'consensus'")
  }
  if (!is.null(gene_sets) && (!is.list(gene_sets) || is.null(names(gene_sets)) || any(is.na(names(gene_sets)) | names(gene_sets) == "") || !all(sapply(gene_sets, is.character)))) {
    stop("gene_sets must be a named list of character vectors (the genes of each gene set)")
  }

  rn <- rownames(counts_norm)
  counts_norm <- apply(as.matrix(counts_norm), 2, as.numeric)
  rownames(counts_norm) <- rn
  results_list <- list()

  # Cohort of each sample (always categorical), in the column order of counts_norm; if named, matched to the samples by name
  if (!is.null(batch)) {
    if (length(batch) != ncol(counts_norm)) {
      stop("batch must have one value per sample of counts_norm (", ncol(counts_norm), " samples)")
    }
    if (!is.null(names(batch))) {
      idx <- match(make.names(colnames(counts_norm)), make.names(names(batch)))
      if (anyNA(idx)) stop("The names of batch do not match the sample names of counts_norm")
      batch <- batch[idx]
    }
    if (anyNA(batch)) stop("batch contains missing values")
    if (length(unique(batch)) < 2) batch <- NULL # A single cohort: same as no batch
  }

  # With batch, these checks are done within each cohort
  if (is.null(batch) && max(counts_norm, na.rm = TRUE) > 20) { # log2(TPM + 1) cannot exceed 20
    warning("counts_norm has values above 20, so it does not look log-transformed: log-normalized expression (e.g. log2(TPM + 1)) is expected")
  }
  if (is.null(batch) && scale && ncol(counts_norm) < 3) {
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

  ###### Batch: pathway activities computed (and scaled) separately in each cohort, then joined in the sample order of counts_norm
  if (!is.null(batch)) {
    cohorts <- split(seq_len(ncol(counts_norm)), factor(batch, levels = unique(batch)))
    per_cohort <- lapply(cohorts, function(idx) {
      res <- compute.pathway.activity(counts_norm[, idx, drop = FALSE], gene_sets = gene_sets, paths = paths, method = method, scale = scale, return = FALSE)
      if (is.data.frame(res)) list(PROGENy = res) else res
    })
    for (name in names(per_cohort[[1]])) {
      common <- Reduce(intersect, lapply(per_cohort, function(res) colnames(res[[name]])))
      not_common <- setdiff(unique(unlist(lapply(per_cohort, function(res) colnames(res[[name]])))), common)
      if (length(not_common) > 0) {
        message(name, " scores left out because they could not be computed in every cohort: ", paste(not_common, collapse = ", "))
      }
      joined <- do.call(rbind, unname(lapply(per_cohort, function(res) res[[name]][, common, drop = FALSE])))
      results_list[[name]] <- joined[colnames(counts_norm), , drop = FALSE]
    }

    if (return) {
      file_suffix <- if (!is.null(file_name)) paste0("_", file_name) else ""
      for (name in names(results_list)) {
        utils::write.csv(results_list[[name]], paste0("Results/Pathway_matrix_", name, file_suffix, ".csv"))
      }
    }
    return(if (length(results_list) == 1) results_list[[1]] else results_list)
  }

  progeny <- decoupleR::decouple(
    mat             = counts_norm,
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
  sample_acts_progeny <- sample_acts_progeny[colnames(counts_norm), , drop = FALSE] # samples in the order of the input

  not_scored <- setdiff(unique(paths$source), colnames(sample_acts_progeny))
  if (length(not_scored) > 0) {
    message("Pathways left out (fewer than 5 of their genes in the data): ", paste(not_scored, collapse = ", "))
  }
  if (scale) sample_acts_progeny <- base::scale(sample_acts_progeny)

  results_list$PROGENy <- as.data.frame(sample_acts_progeny)

  ###### GSVA (optional)
  if (!is.null(gene_sets)) {
    gsva_results <- GSVA::gsva(
      GSVA::gsvaParam(counts_norm, gene_sets, kcdf = "Gaussian", minSize = 1, maxDiff = TRUE),
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
    file_suffix <- if (!is.null(file_name)) paste0("_", file_name) else ""
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
#' @param counts_norm A gene expression matrix with genes (symbols) as rows and samples as columns, with unique sample
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
#'   `input/ARACNE/<cancer_type>/network/network.txt` (relative to the working directory). The mode of
#'   regulation of each edge is read from a `mor` column of the file (1 activation, -1 repression), which should be
#'   the sign of the Spearman correlation between TF and target computed on the cohort used to build the network.
#'   If the file has no `mor` column, a warning is given and the signs are computed as the sign of the Spearman
#'   correlation between TF and target in the samples of `counts_norm`; they then depend on the samples given (e.g.
#'   they can differ between cohorts, or between a training and a test set).
#' @param min_targets_size Integer. Minimum number of target genes per regulon (passed to
#'   \code{decoupleR::decouple()} as \code{minsize}). Default is 5.
#' @param universe Optional. A user-specified data frame of TF-target interactions (columns \code{source},
#'   \code{target}, \code{mor}). The three columns are required: \code{mor} is the mode of regulation
#'   (1 activation, -1 repression; use 1 for every edge if the network has no sign). If not provided, the
#'   network is fetched based on `TF.collection`. Ignored when \code{TF.collection = "ARACNE"}.
#' @param cancer_type Optional character. TCGA cancer type abbreviation used to locate the ARACNe network
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
#' @param file_name Optional character suffix used when writing the TF activity matrix to disk
#'   (\code{Results/TF_matrix_<file_name>.csv}; without it, \code{Results/TF_matrix.csv}).
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
compute.TFs.activity <- function(counts_norm, TF.collection = "CollecTRI", min_targets_size = 5, universe = NULL, cancer_type = NULL, method = "ulm", scale = TRUE, return = TRUE, file_name = NULL){

  if(!method %in% c("ulm", "mlm", "wsum", "consensus")){
    stop("method must be 'ulm', 'mlm', 'wsum' or 'consensus'")
  }
  if(!TF.collection %in% c("CollecTRI", "Dorothea", "ARACNE") && is.null(universe)){
    stop("TF.collection must be 'CollecTRI', 'Dorothea' or 'ARACNE'")
  }
  if(!is.null(universe) && TF.collection != "ARACNE" && !all(c("source", "target", "mor") %in% colnames(universe))){
    stop("universe must have the columns 'source', 'target' and 'mor' (mode of regulation: 1 activation, -1 repression; use 1 for every edge if the network has no sign)")
  }
  if(is.null(colnames(counts_norm)) || anyDuplicated(colnames(counts_norm)) > 0){
    stop("counts_norm must have unique sample names as column names")
  }
  if(max(counts_norm, na.rm = TRUE) > 20){ # log2(TPM + 1) cannot exceed 20
    warning("counts_norm has values above 20, so it does not look log-transformed: log-normalized expression (e.g. log2(TPM + 1)) is expected")
  }
  if(scale && ncol(counts_norm) < 3){
    message("Fewer than three samples: the TF activity scores are not scaled")
    scale <- FALSE
  }

  tf_cache_file <- paste0("Results/TF_target_collection_", TF.collection, ".csv")

  if(TF.collection == "ARACNE"){

    if(is.null(cancer_type)){
      # auto-discover when only one network exists
      candidates <- list.files("input/ARACNE", pattern = "^network\\.txt$",
                               recursive = TRUE, full.names = TRUE)
      if(length(candidates) == 0)
        stop("TF.collection = 'ARACNE' requires a 'cancer_type' or a network.txt under input/ARACNE/")
      if(length(candidates) > 1)
        stop("Multiple ARACNe networks found. Specify 'cancer_type' (e.g. cancer_type = 'skcm'):\n",
             paste(dirname(dirname(candidates)), collapse = "\n"))
      aracne.network <- candidates[1]
      message("Auto-detected ARACNe network: ", aracne.network)
    } else {
      aracne.network <- file.path("input/ARACNE", cancer_type, "network", "network.txt")
      if(!file.exists(aracne.network))
        stop("ARACNe network not found for cancer type '", cancer_type, "': ", aracne.network)
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
        dplyr::filter(target %in% rownames(counts_norm) & !is.na(mor) & mor != 0)
    } else {
      warning("The ARACNe network file has no 'mor' column (sign of each interaction): the signs are computed as the sign of the Spearman correlation ",
              "between each TF and its target in the samples given, so they depend on these samples. Add a 'mor' column (1 activation, -1 repression) ",
              "computed on the cohort used to build the network: ", aracne.network)
      aracne_net <- aracne_net %>%
        dplyr::select(source = Regulator, target = Target) %>%
        dplyr::filter(source %in% rownames(counts_norm) & target %in% rownames(counts_norm))

      # Compute Spearman correlation between every TF and its targets in one matrix op
      # (this is exactly what TFmode1 does internally)
      all_tfs     <- unique(aracne_net$source)
      all_targets <- unique(aracne_net$target)
      cor_mat <- suppressWarnings(
        stats::cor(t(counts_norm[all_tfs, , drop = FALSE]),
                   t(counts_norm[all_targets, , drop = FALSE]),
                   method = "spearman")
      )

      # mor = sign of Spearman correlation - +1 activation, -1 repression
      universe <- aracne_net %>%
        dplyr::mutate(mor = cor_mat[cbind(source, target)]) %>%
        dplyr::mutate(mor = sign(mor)) %>%
        dplyr::filter(!is.na(mor) & mor != 0)
    }

    message("Computing TF activities...")

    sample_acts <- decoupleR::decouple( mat     = counts_norm,
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

    sample_acts <- decoupleR::decouple(mat     = counts_norm,
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

  sample_acts <- sample_acts[colnames(counts_norm), , drop = FALSE]
  sample_acts <- if (scale) base::scale(sample_acts) else sample_acts

  result <- data.frame(sample_acts)
  colnames(result) <- make.names(colnames(result))

  if(return){
    utils::write.csv(result, paste0("Results/TF_matrix", if(!is.null(file_name)) paste0("_", file_name), ".csv"))
  }

  return(result)

}

#' Compute Weighted TF-coactivity Network Analysis (WTCNA)
#'
#' Construct a weighted signed or unsigned network using TF activity to cluster protein regulators
#' into modules that share similar activity patterns. Each TF module will have a sample-level score
#' represented by the eigengene of the module.
#'
#' @param tfs Matrix of TF activity (samples x TFs), scaled or not: the TFs are standardized (within each cohort with
#'   \code{batch}) to compute the module eigengenes, and the other steps use correlations, which do not depend on the
#'   scaling. With \code{batch}, the TF activities of all the cohorts, each computed within its cohort (as
#'   \code{CellTFusion()} does).
#' @param batch (Optional) A vector with the cohort of each sample, in the row order of \code{tfs}; if it is
#'   named, it is matched to the samples by name instead. It is always treated as categorical. If given, consensus TF
#'   modules are built across the cohorts (see Details) and module eigengenes are computed and scaled within each
#'   cohort. Default is \code{NULL} (single cohort); a vector with only one cohort is treated the same way.
#' @param network_type Network type: "signed" (default) or "unsigned".
#' @param clustering_method Linkage method of the hierarchical clustering of the TFs (\code{stats::hclust()}), applied to the
#'   TOM dissimilarity (1 - topological overlap) between TFs, or to the consensus TOM dissimilarity in batch mode. Default
#'   is "ward.D2" (WGCNA's own default is "average": two groups of TFs are joined according to the average dissimilarity
#'   between their TFs).
#' @param min_module_size Minimum number of TFs per module. Default is 5.
#' @param corr_mod Correlation threshold (0-1): modules whose eigengenes have a correlation (of type \code{corr_type})
#'   above it are merged. In batch mode, the correlation is computed on the eigengenes of all the samples, scaled within
#'   each cohort. Default is 0.9.
#' @param corr_type Correlation used to pick the soft-threshold power, build the adjacency matrix, merge modules and
#'   compute the proportion of variance explained: "spearman" (default) or "pearson", in single-cohort and batch mode. The
#'   module eigengene is always the first principal component of the module's TF activities.
#' @param soft_power Optional numeric value specifying the soft-thresholding power used to build the adjacency
#'   matrix (one value per cohort, in the order the cohorts first appear in \code{batch}, in batch mode). If \code{NULL}, the lowest power whose scale-free fit
#'   \eqn{R^2} reaches 0.9 is chosen automatically (the power with the best fit if none reaches 0.9, with a message).
#' @param consensus_quantile Numeric between 0 and 1 (batch mode only). For each pair of TFs, the consensus topological
#'   overlap is this quantile of its (calibrated) topological overlaps across the cohorts. With 0 (default), the
#'   minimum: two TFs are close only if they are close in every cohort, the strictest choice, which removes more TF
#'   relationships as more cohorts are added. Higher values are less strict: with \code{K} cohorts,
#'   \code{consensus_quantile = k / (K - 1)} takes the \code{(k + 1)}-th lowest value, so two TFs are close if they are
#'   close in all cohorts but \code{k} (e.g. 0.25 with 5 cohorts: in at least 4 of them); 0.5 is the median.
#' @param verbose Logical; whether to print the function messages. Default is TRUE.
#' @param file_name Optional character suffix used when writing WTCNA outputs (CSV files and, in single-cohort
#'   mode, the diagnostic plots \code{Soft_Threshold_<file_name>.pdf},
#'   \code{Gene_dendrogram_and_module_colors_<file_name>.pdf} and
#'   \code{Gene_dendrogram_and_module_colors_after_merging_<file_name>.pdf}).
#' @param return Logical, whether to save output plots, the module list (\code{Results/TFs_modules_<file_name>.csv})
#'   and the module matrix (\code{Results/TF_module_matrix_<file_name>.csv}, with the sample names as row names)
#'   to "Results/" (without \code{file_name}, \code{TFs_modules.csv} and \code{TF_module_matrix.csv}). Default is TRUE.
#'
#' @return A named list with:
#' \itemize{
#'   \item \code{TFs module matrix}: Scaled module eigengenes (samples x modules), in the row order of \code{tfs}.
#'   \item \code{TFs colors}: Vector of module colors assigned to each TF (\code{"grey"} for TFs not assigned to any module).
#'   \item \code{TFs per module}: List of TF names in each module.
#'   \item \code{Proportion of variance}: Variance explained per module (single-cohort mode only).
#'   \item \code{TFs_matrix}: The TF activity matrix used (a list of per-cohort matrices in batch mode, in the order the cohorts first appear in \code{batch}).
#'   \item \code{Network parameters}: A list with the soft-thresholding power(s) used (\code{soft_power}, one per cohort in
#'     batch mode), \code{network_type} and \code{corr_type}, used by \code{identify_hub_TFs()}, and in batch mode
#'     \code{consensus_quantile}.
#' }
#'
#' @details
#' TFs that WGCNA does not assign to any module (color \code{"grey"}) are not a module: they are left out of the
#' module matrix and of \code{TFs per module}, with a message.
#'
#' In batch mode, the modules are built with the same steps as in single-cohort mode, on a consensus network:
#' \enumerate{
#'   \item the samples are split by cohort, and a soft-thresholding power is chosen for each cohort;
#'   \item the adjacency and the topological overlap matrix (TOM) are computed in each cohort, with its own power;
#'   \item the TOMs are calibrated to a common scale as in \code{WGCNA::blockwiseConsensusModules()}: each TOM is raised to
#'     the power that makes its 0.95 quantile equal to that of the first cohort;
#'   \item the consensus TOM is, for each pair of TFs, the \code{consensus_quantile} quantile of its calibrated TOMs across
#'     cohorts;
#'   \item the TFs are clustered on 1 - consensus TOM with \code{clustering_method} and the tree is cut with the same dynamic
#'     tree cut as in single-cohort mode;
#'   \item the module eigengenes are computed and scaled within each cohort, and the modules are merged as in single-cohort
#'     mode (see \code{corr_mod}); merged scores are scaled within each cohort again.
#' }
#' Unlike \code{WGCNA::blockwiseConsensusModules()}, the TFs of each module are not filtered or reassigned afterwards
#' by their module membership (as in single-cohort mode). The diagnostic plots are only saved in single-cohort mode.
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
compute.WTCNA <- function(tfs, batch = NULL, network_type = "signed", clustering_method = "ward.D2",
                          min_module_size = 5, corr_mod = 0.9, corr_type = "spearman", verbose = TRUE, file_name = NULL,
                          soft_power = NULL, consensus_quantile = 0, return = TRUE) {

  if(verbose){
    cat("Creating weighted TF-coactivity network......................................................................\n\n")
  }

  if(!network_type %in% c("signed", "unsigned")){
    stop("network_type must be 'signed' or 'unsigned'")
  }
  if(!corr_type %in% c("spearman", "pearson")){
    stop("corr_type must be 'spearman' or 'pearson'")
  }
  if(!all(sapply(as.data.frame(tfs), is.numeric))){
    stop("tfs must contain only numeric TF activity values")
  }

  # Cohort of each sample (always categorical), in the row order of tfs; if named, matched to the samples by name
  if(!is.null(batch)){
    if(length(batch) != nrow(tfs)){
      stop("batch must have one value per sample (row) of tfs (", nrow(tfs), " samples)")
    }
    if(!is.null(names(batch))){
      idx = match(make.names(rownames(tfs)), make.names(names(batch)))
      if(anyNA(idx)) stop("The names of batch do not match the sample names of tfs")
      batch = batch[idx]
    }
    if(anyNA(batch)) stop("batch contains missing values")
    if(length(unique(batch)) < 2) batch = NULL # A single cohort: same as no batch
  }

  if(!is.null(batch)){
    if(!is.numeric(consensus_quantile) || length(consensus_quantile) != 1 || is.na(consensus_quantile) || consensus_quantile < 0 || consensus_quantile > 1){
      stop("consensus_quantile must be a number between 0 and 1")
    }
    cor_method = corr_type

    # Step 1: split the samples by cohort (cohorts in order of first appearance); the module scores are returned in the row order of tfs
    cohort_idx <- split(seq_len(nrow(tfs)), factor(batch, levels = unique(batch)))
    tfs <- lapply(cohort_idx, function(i) tfs[i, , drop = FALSE])
    nCohorts <- length(tfs)
    if(verbose) cat("Running consensus WGCNA across", nCohorts, "cohorts...\n")

    # Step 2: pick soft-threshold per cohort if not provided
    powers = c(1:10, 12:20)
    if(is.null(soft_power)){
      soft_power <- numeric(nCohorts)
      for(i in 1:nCohorts){
        if(verbose) cat("Picking soft threshold for cohort", i, "...\n")
        utils::capture.output(sft <- WGCNA::pickSoftThreshold(tfs[[i]], powerVector = powers, verbose = 0, networkType = network_type,
                                                              corOptions = list(use = "p", method = cor_method))) # it prints its table otherwise
        target = 0.9
        fit = -sign(sft$fitIndices[,3]) * sft$fitIndices[,2] # signed scale-free fit R^2
        reached = which(fit >= target)
        soft_power[i] <- if(length(reached) > 0) powers[reached[1]] else powers[which.max(fit)] # lowest power reaching the target
        if(length(reached) == 0) message("Cohort ", i, ": the scale-free fit never reaches ", target, " (best ", round(max(fit, na.rm = TRUE), 2), " at power ", soft_power[i], "), this power is used")
        if(verbose) cat("Cohort", i, "soft_power =", soft_power[i], "\n")
      }
    } else {
      if(length(soft_power) != nCohorts){
        stop("soft_power must be a vector with one value per cohort")
      }
    }

    # Step 3: TOM of each cohort (with its own power), calibrated to the scale of the first cohort as in
    # WGCNA::blockwiseConsensusModules() ("single quantile": the 0.95 quantile of each TOM is matched to that of the first cohort)
    TOMs <- lapply(seq_len(nCohorts), function(i){
      adjacency = WGCNA::adjacency(tfs[[i]], power = soft_power[i], type = network_type, corFnc = "cor", corOptions = list(use = "p", method = cor_method))
      WGCNA::TOMsimilarity(adjacency, TOMType = network_type, verbose = 0)
    })
    q_ref <- stats::quantile(TOMs[[1]][upper.tri(TOMs[[1]])], probs = 0.95, type = 8)
    for(i in seq_len(nCohorts)[-1]){
      q <- stats::quantile(TOMs[[i]][upper.tri(TOMs[[i]])], probs = 0.95, type = 8)
      TOMs[[i]] <- TOMs[[i]]^(log(q_ref)/log(q))
    }

    # Step 4: consensus TOM: for each pair of TFs, the consensus_quantile quantile of its calibrated TOMs across cohorts (0 = minimum)
    consTOM <- do.call(WGCNA::pquantile, c(list(prob = consensus_quantile), unname(TOMs)))
    dissTOM <- 1 - consTOM

    # Step 5: clustering and dynamic tree cut of the consensus dissimilarity, as in single-cohort mode
    geneTree = stats::hclust(stats::as.dist(dissTOM), method = clustering_method)
    dynamicMods = dynamicTreeCut::cutreeDynamic(dendro = geneTree, distM = dissTOM,
                                                deepSplit = 2, pamRespectsDendro = FALSE,
                                                minClusterSize = min_module_size, verbose = 0);
    dynamicColors = WGCNA::labels2colors(dynamicMods)
    rm(TOMs, consTOM, dissTOM)
    gc()

    # Step 6: module eigengenes computed and scaled within each cohort, then stacked cohort by cohort (grey = TFs without module, not a module)
    MEs = do.call(rbind, lapply(tfs, function(x){
      ME = WGCNA::moduleEigengenes(x, colors = dynamicColors, scale = TRUE)$eigengenes ## TFs standardized within the cohort (no change if they are already scaled by compute.TFs.activity())
      rownames(ME) = rownames(x)
      scale(ME)
    }))
    MEs = WGCNA::orderMEs(MEs)
    MEs = MEs[, colnames(MEs) != "MEgrey", drop = FALSE]
    if(ncol(MEs) == 0) stop("No TF module found: no TF was assigned to a consensus module. Try a smaller min_module_size")

    # Step 7: merging of the modules as in single-cohort mode, on the stacked eigengenes; merged scores are scaled within each cohort again
    merge = mergeModules(MEs, dynamicColors, corr_mod, corr_type = cor_method)
    MEs = merge[[1]]
    dynamicColors = merge[[2]]
    cohort_rows = rep(seq_len(nCohorts), sapply(tfs, nrow))
    for(i in seq_len(nCohorts)){
      MEs[cohort_rows == i, ] = scale(MEs[cohort_rows == i, , drop = FALSE])
    }

    # Step 8: Map TFs to modules (same as single-cohort mode)
    modtfs = list()
    tf_names = colnames(tfs[[1]])
    modules = c(substring(names(MEs), 3))
    for(i in seq_along(modules)){
      modtfs[[i]] = tf_names[dynamicColors == modules[i]]
    }
    names(modtfs) = modules
    colnames(MEs) <- gsub("ME", "", colnames(MEs))
    MEs = MEs[order(unlist(cohort_idx)), , drop = FALSE] # Back to the row order of tfs

    output <- list(
      "TFs module matrix" = MEs,
      "TFs colors" = dynamicColors,
      "TFs per module" = modtfs,
      "TFs_matrix" = tfs,
      "Network parameters" = list(soft_power = soft_power, network_type = network_type, corr_type = corr_type, consensus_quantile = consensus_quantile)
    )


  } else {


    #####Choose parameter for scale-free network topology
    powers = c(c(1:10), seq(from = 12, to=20, by=1))
    utils::capture.output(sft <- WGCNA::pickSoftThreshold(tfs, powerVector = powers, verbose = 0, networkType = network_type,
                                                          corOptions = list(use = "p", method = corr_type))) # it prints its table otherwise

    if(return){
      pdf(paste0("Results/Soft_Threshold", if(!is.null(file_name)) paste0("_", file_name), ".pdf"))
      plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
           xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
           main = paste("Scale independence"))
      graphics::text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=0.9,col="red");
      graphics::abline(h=0.90,col="red")
      dev.off()
    }

    # Automatic or user-defined soft-threshold selection
    if (is.null(soft_power)) {
      target = 0.9 # Target SFT.R.sq value
      fit = -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2] # signed scale-free fit R^2
      reached = which(fit >= target)
      soft_power = if(length(reached) > 0) powers[reached[1]] else powers[which.max(fit)] # lowest power reaching the target
      if(length(reached) == 0) message("The scale-free fit never reaches ", target, " (best ", round(max(fit, na.rm = TRUE), 2), " at power ", soft_power, "): this power is used")

      if (verbose) {
        cat("Automatically choosing", soft_power, "as soft-threshold......................................................................\n\n")
      }
    } else {
      if (verbose) {
        cat("Using user-defined soft-threshold =", soft_power, "......................................................................\n\n")
      }
    }

    if(verbose){
      #####Co-expression matrix using nodes adjacency and topological overlapping nodes
      cat("Calculating nodes adjacency and topological overlapping nodes.................................................\n\n")
    }

    adjacency = WGCNA::adjacency(tfs, power =soft_power, type=network_type, corFnc = "cor", corOptions = list(use = "p", method = corr_type))
    TOM = WGCNA::TOMsimilarity(adjacency, TOMType = network_type, verbose = 0)
    dissTOM = 1-TOM

    #####Unsupervised hierarchical clustering using dissimilarity matrix
    geneTree = stats::hclust(stats::as.dist(dissTOM), method = clustering_method) # Tree built on the TOM dissimilarity itself, as in WGCNA
    dynamicMods = dynamicTreeCut::cutreeDynamic(dendro = geneTree, distM = dissTOM,
                                                deepSplit = 2, pamRespectsDendro = FALSE,
                                                minClusterSize = min_module_size, verbose = 0);
    dynamicColors = WGCNA::labels2colors(dynamicMods)

    #####Remove variables and clean garbage
    rm(adjacency, TOM, dissTOM)
    gc()

    if(return){
      pdf(paste0("Results/Gene_dendrogram_and_module_colors", if(!is.null(file_name)) paste0("_", file_name), ".pdf"))
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

    MEList = WGCNA::moduleEigengenes(tfs, colors = dynamicColors, scale = TRUE) ## TFs standardized (no change if they are already scaled by compute.TFs.activity())
    MEs = scale(MEList$eigengenes) # Scale MEs scores (not the same as above) and this is to have standard ranges and same input with or without batch (for this case)
    MEs = WGCNA::orderMEs(MEs)
    MEs = MEs[, colnames(MEs) != "MEgrey", drop = FALSE] # grey = TFs without module, not a module
    if(ncol(MEs) == 0) stop("No TF module found: no TF was assigned to a module. Try a smaller min_module_size")

    if(verbose){
      print(paste0("Merging modules significantly correlated with ", corr_mod, "........"))
    }

    merge = mergeModules(MEs, dynamicColors, corr_mod, corr_type = corr_type)
    MEs = merge[[1]]
    dynamicColors = merge[[2]]
    modtfs = list()
    for(i in 1:ncol(MEs)){
      tf_names = colnames(tfs)
      modules = c(substring(names(MEs), 3))
      inModule = is.finite(match(dynamicColors,modules[i]))
      modtfs[[i]] = tf_names[inModule]
    }
    names(modtfs) = modules

    if(return){
      pdf(paste0("Results/Gene_dendrogram_and_module_colors_after_merging", if(!is.null(file_name)) paste0("_", file_name), ".pdf"))
      WGCNA::plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                                 dendroLabels = FALSE, hang = 0.03,
                                 addGuide = TRUE, guideHang = 0.05,
                                 main = "Gene dendrogram and module colors")
      dev.off()
    }

    TFspropVar = WGCNA::propVarExplained(tfs[, dynamicColors != "grey", drop = FALSE], dynamicColors[dynamicColors != "grey"], MEs, corFnc = "cor", corOptions = if(corr_type == "pearson") "use = 'p'" else "use = 'p', method = 'spearman'")

    colnames(MEs) <- gsub("ME", "", colnames(MEs))

    output = list(MEs, dynamicColors, modtfs, TFspropVar, tfs, list(soft_power = soft_power, network_type = network_type, corr_type = corr_type))
    names(output) = c("TFs module matrix", "TFs colors", "TFs per module", "Proportion of variance", "TFs_matrix", "Network parameters")

  }

  unassigned = colnames(if(!is.null(batch)) tfs[[1]] else tfs)[dynamicColors == "grey"]
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
    file_suffix = if(!is.null(file_name)) paste0("_", file_name) else ""
    utils::write.csv(tfs_modules, paste0("Results/TFs_modules", file_suffix, ".csv"), row.names = F)
    utils::write.csv(MEs, file = paste0("Results/TF_module_matrix", file_suffix, ".csv")) # Sample names as row names
  }

  return(output)
}

#' Identify hub TFs
#'
#' Identifies hub TFs per module using values of module membership and degree.
#' TFs with high module membership (correlation of their activity with the module score > \code{mm_thresh}, with the
#' correlation used to build the network, \code{corr_type} of \code{compute.WTCNA()})
#' and a degree at or above the \code{degree_thresh} quantile of their module are considered hub TFs.
#' The degree of a TF is its intramodular connectivity, as in WGCNA: the sum of its network adjacencies with the
#' other TFs of the same module. The adjacencies are computed as in \code{compute.WTCNA()}, with the same
#' soft-thresholding power, network type and correlation (for a \code{"signed"} network, \eqn{((1 + r)/2)^{power}}, so a
#' negative correlation adds almost nothing instead of being subtracted).
#'
#' The TF activities are taken from \code{network} (element \code{TFs_matrix}) and matched to the module scores
#' by sample name. For a multi-cohort network (\code{compute.WTCNA()} with \code{batch}), the per-cohort TF activities
#' are stacked to compute the membership, and the degree is the mean of the per-cohort degrees (each with the power
#' of its cohort).
#'
#' @param network TF network obtained from \code{compute.WTCNA()}.
#' @param mm_thresh Threshold for module membership. Default is 0.8.
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
#' Networks built before the soft-thresholding power was stored in \code{network} (element \code{Network parameters})
#' have their power chosen again with the rule of \code{compute.WTCNA()} (with a message), assuming a signed network
#' and Pearson correlation.
#'
#' @export
#'
#' @examples
#'
#' data("network.tuto")
#'
#' hub_tfs <- identify_hub_TFs(network.tuto, mm_thresh = 0.8, degree_thresh = 0.9)
#'
identify_hub_TFs <- function(network, mm_thresh = 0.8, degree_thresh = 0.9) {

  if(!is.list(network) || is.data.frame(network) || is.null(network[["TFs module matrix"]]) || is.null(network[["TFs_matrix"]])){
    stop("network must be the output of compute.WTCNA(): the TF activities are taken from it (the datExpr argument was removed)")
  }
  moduleEigengenes_df = network[["TFs module matrix"]]
  moduleColors = network[["TFs colors"]]
  modules = intersect(unique(moduleColors), colnames(moduleEigengenes_df)) # TFs without module (grey) are not a module

  # TF activities: one matrix (samples x TFs), or one per cohort for a multi-cohort network
  tfs_list = network[["TFs_matrix"]]
  if(is.data.frame(tfs_list) || is.matrix(tfs_list)) tfs_list = list(tfs_list)
  tfs_list = lapply(tfs_list, as.matrix)
  tf_names = colnames(tfs_list[[1]])
  if(length(moduleColors) != length(tf_names)){
    stop("network is inconsistent: ", length(moduleColors), " module colors for ", length(tf_names), " TFs")
  }

  # Samples matched to the module scores by name (cohort matrices stacked in batch mode). TF activities standardized within
  # each cohort, as the module scores, so that the membership does not depend on whether they were scaled
  tfs_all = do.call(rbind, unname(lapply(tfs_list, function(x){ z = scale(x); z[is.na(z)] = 0; z })))
  idx = match(make.names(rownames(moduleEigengenes_df)), make.names(rownames(tfs_all)))
  if(is.null(rownames(moduleEigengenes_df)) || is.null(rownames(tfs_all)) || anyNA(idx)){
    stop("The samples of the module matrix and of the TF activities in network do not match (sample names are needed)")
  }
  tfs_all = tfs_all[idx, , drop = FALSE]

  # Network parameters of compute.WTCNA() (chosen again for networks built before they were stored)
  params = network[["Network parameters"]]
  if(is.null(params)){
    powers = c(1:10, 12:20)
    soft_power = sapply(tfs_list, function(x){
      utils::capture.output(sft <- WGCNA::pickSoftThreshold(x, powerVector = powers, verbose = 0, networkType = "signed"))
      fit = -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
      reached = which(fit >= 0.9)
      if(length(reached) > 0) powers[reached[1]] else powers[which.max(fit)]
    })
    params = list(soft_power = soft_power, network_type = "signed", corr_type = "pearson")
    message("network does not store its soft-thresholding power (built with an older version of compute.WTCNA()): power ",
            paste(soft_power, collapse = ", "), " chosen again, signed network and Pearson correlation assumed")
  }
  cor_method = params$corr_type

  # Identify hub genes
  hubGenesList <- list()
  hubGenesData <- data.frame()  # Initialize an empty data frame

  for (module in modules) {
    inModule <- which(moduleColors == module)

    # Module membership: correlation of each TF with the module score
    memberships <- stats::cor(tfs_all[, inModule, drop = FALSE], moduleEigengenes_df[, module], method = cor_method)

    # Degree: intramodular connectivity as in WGCNA (sum of the adjacencies with the other TFs of the module), mean over cohorts
    degrees_per_cohort <- sapply(seq_along(tfs_list), function(i){
      adjacency <- WGCNA::adjacency(tfs_list[[i]][, inModule, drop = FALSE], power = params$soft_power[i], type = params$network_type,
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
    highMMGenes <- which(memberships > mm_thresh)
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
#' @param relationship The list returned by \code{compute.modules.relationship()}. Its first two elements are used
#'   (accessed by position):
#'   \itemize{
#'     \item \code{[[1]]} (\code{correlation}): A matrix of correlations between TF modules (rows, named) and deconvolution features (columns).
#'     \item \code{[[2]]} (\code{significant}): A list, named by TF module, of the features significantly correlated with each module
#'       (p-value, or adjusted p-value with \code{padj = TRUE}, <= the \code{pval} given to \code{compute.modules.relationship()}).
#'   }
#'   The correlations of each module are taken from the row with its name.
#' @param clustering_method Linkage method of the hierarchical clustering (\code{stats::hclust()}) of the features
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
#' rel <- compute.modules.relationship(network.tuto[[1]], deconv_subgroups.tuto[[1]], return = FALSE)
#' dendrograms <- identify.cell.groups(rel, return = FALSE)
#' names(dendrograms)
#'
identify.cell.groups = function(relationship, clustering_method = "ward.D2", file_name = NULL, return = TRUE){

  if(!is.list(relationship) || length(relationship) < 2 || !is.matrix(relationship[[1]]) || !is.numeric(relationship[[1]]) || !is.list(relationship[[2]]) ||
     is.null(rownames(relationship[[1]])) || is.null(names(relationship[[2]]))){
    stop("relationship must be the list returned by compute.modules.relationship(): a correlation matrix with TF modules as row names, ",
         "and a list of significant features named by TF module")
  }
  moduleTraitCor = relationship[[1]]
  features_vec = relationship[[2]]
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
    dendrogram <- stats::hclust(d, method = clustering_method)
    dendrogram$correlation <- TFmoduleTraitcor[1, dendrogram$labels] # Kept for calculate_dendrogram_cuts(): each cell group must be one range of these values
    if(return){
      # Labels with the correlation of each feature; size from the number of features and the length of their labels
      dendrogram_plot <- dendrogram
      dendrogram_plot$labels <- paste0(dendrogram$labels, " (r = ", sprintf("%.2f", TFmoduleTraitcor[1, dendrogram$labels]), ")")
      label_lines <- max(nchar(dendrogram_plot$labels)) * 0.62
      pdf(paste0("Results/Dendrogram_cell_types_", module, if(!is.null(file_name)) paste0("_", file_name), ".pdf"),
          width = 4 + label_lines * 0.2 * 0.8, height = max(4, 1.5 + 0.16 * length(dendrogram$labels)))
      par(mar = c(5, 1, 4, label_lines), cex = 0.8) #bottom, left, top, right
      plot(as.dendrogram(dendrogram_plot), horiz = T, xlab = paste0("Merge height (", clustering_method, ") on differences in correlation"),
           main = paste0("TF module ", module, ": ", length(dendrogram$labels), " significantly correlated features"),
           adj = 0, cex.main = 1) # title starts at the left edge of the tree and runs over the labels, so it is not cut
      dev.off()
    }
    lis.dendrogram[[module]] = dendrogram
  }

  return(lis.dendrogram)

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
#'     score of its TF module (\code{cell.groups.computation()}).
#' }
#' The deconvolution matrix and \code{network} must have the same samples (matched by name; their order can differ): the
#' module scores of \code{network} are put in the sample order of the deconvolution matrix.
#'
#' @param network A TF module network as returned by \code{compute.WTCNA()} (its module scores are used).
#' @param dt Deconvolution subgroups as returned by \code{multideconv::compute.deconvolution.analysis()}; its element
#'   \code{Deconvolution matrix} (samples x features) is used. It must have the same samples as \code{network}.
#' @param batch Optional vector with the cohort of each sample: in the row order of the deconvolution matrix, or named by
#'   sample. Required when \code{network} was built with \code{compute.WTCNA()} with \code{batch}. With two or more cohorts,
#'   correlations and scores are corrected for cohort.
#' @param pval Numeric. P-value threshold used to select the deconvolution features significantly
#'   correlated with each TF module (adjusted p-value with \code{padj = TRUE}). Default: 0.05.
#' @param padj Logical. If TRUE, the p-values of each TF module are corrected for multiple testing across the
#'   deconvolution features (Benjamini-Hochberg, see \code{compute.modules.relationship()}) before selecting the
#'   significant features. Default: FALSE.
#' @param corr_type Correlation used between the TF modules and the deconvolution features (passed to
#'   \code{compute.modules.relationship()}): \code{"spearman"} (default) or \code{"pearson"}. It is used to
#'   select the significant features and to cluster them.
#' @param clustering_method Linkage method of the hierarchical clustering of the deconvolution features of each module
#'   on their correlation with it (see \code{identify.cell.groups()}). Default: "ward.D2".
#' @param file_name Optional character. Suffix of the output file names (see \code{return}).
#' @param return Logical. If TRUE (default), saves a PDF of the colored cell-group dendrograms
#'   (\code{Results/Dendrogram_color_clusters_<file_name>.pdf}) and the composition and the scores of the cell groups
#'   (\code{Results/Cell.groups.composition_<file_name>.csv} and \code{Results/Cell.groups.scores_<file_name>.csv}).
#' @param min_cluster_size Integer. Minimum number of features of a cell group (minimum cluster size of the
#'   dynamic tree cut of each dendrogram). Default is 3.
#'
#' @return A named list of 3 elements:
#' \describe{
#'   \item{Cell_groups}{A data frame with the projected cell group scores (samples x groups). A positive score
#'     means that the sample is on the active side of the TF module of that cell group, a negative score that
#'     it is on the inactive side.}
#'   \item{Composition}{A named list where each element is a character vector of the deconvolution features in each group.}
#'   \item{Weights}{A named list of CCA projection parameters (\code{xcoef}, \code{train_means}, \code{train_sds}) for each group.
#'     With \code{batch}, only \code{xcoef} is kept: test data are z-scored on their own means and SDs, as each training
#'     cohort was.}
#' }
#' Cell groups are named \code{Dendrogram_<module>_group_<n>}, where \code{n} is the cluster number shown in the
#' dendrogram PDF of that module. Features that the tree cut does not assign to any cluster are not part of any
#' cell group (a message reports how many).
#' @examples
#' data("network.tuto")
#' data("deconv_subgroups.tuto")
#' cell_groups <- construct_cell_groups(network.tuto, deconv_subgroups.tuto)
#' cell_groups$Cell_groups[1:3, 1:3]
#'
#' @export
construct_cell_groups = function(network, dt, batch = NULL, pval = 0.05, padj = FALSE, corr_type = "spearman", clustering_method = "ward.D2", file_name = NULL, return = TRUE, min_cluster_size = 3){

  if(!is.list(network) || is.data.frame(network) || is.null(network[["TFs module matrix"]])){
    stop("network must be the output of compute.WTCNA()")
  }
  if(!is.list(dt) || is.data.frame(dt) || is.null(dt[["Deconvolution matrix"]])){
    stop("dt must be the output of multideconv::compute.deconvolution.analysis() (element 'Deconvolution matrix')")
  }
  deconv = dt[["Deconvolution matrix"]]
  samples = make.names(rownames(deconv))

  # Multi-cohort network (compute.WTCNA() with batch): TF activities are one matrix per cohort
  multi_cohort = is.list(network[["TFs_matrix"]]) && !is.data.frame(network[["TFs_matrix"]])
  if(multi_cohort && is.null(batch)){
    stop("network was built with compute.WTCNA() with batch: give the cohort of each sample in batch")
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

  # Module scores in the sample order of the deconvolution matrix (matched by name): both must have the same samples
  network_samples = make.names(rownames(network[["TFs module matrix"]]))
  idx_modules = match(samples, network_samples)
  missing_samples = rownames(deconv)[is.na(idx_modules)]
  extra_samples = rownames(network[["TFs module matrix"]])[!network_samples %in% samples]
  if(length(missing_samples) > 0 || length(extra_samples) > 0){
    stop("The deconvolution matrix and the TF network must have the same samples: ",
         if(length(missing_samples) > 0) paste0(length(missing_samples), " samples of the deconvolution matrix are not in the network (e.g. ", paste(utils::head(missing_samples, 3), collapse = ", "), ")"),
         if(length(missing_samples) > 0 && length(extra_samples) > 0) "; ",
         if(length(extra_samples) > 0) paste0(length(extra_samples), " samples of the network are not in the deconvolution matrix (e.g. ", paste(utils::head(extra_samples, 3), collapse = ", "), ")"))
  }
  network[["TFs module matrix"]] = network[["TFs module matrix"]][idx_modules, , drop = F]

  corr_modules = compute.modules.relationship(network[["TFs module matrix"]], deconv, batch = batch, return = FALSE, pval = pval, padj = padj, corr_type = corr_type)
  cell_dendrograms = identify.cell.groups(corr_modules, clustering_method = clustering_method, return = F)
  if(is.null(cell_dendrograms)){
    stop("No cell groups found: no TF module is significantly correlated (", if(padj) "adjusted p-value" else "p-value", " <= ", pval, ") with at least two deconvolution features")
  }
  cell.groups = cell.groups.computation(deconv, cell_dendrograms, network, batch = batch, return = return, file_name = file_name, min_cluster_size = min_cluster_size)

  names(cell.groups) = c("Cell_groups", "Composition", "Weights")
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
#' @param dt Deconvolution subgroups of the training set, as returned by
#'   \code{multideconv::compute.deconvolution.analysis()} (the elements \code{Deconvolution matrix} and
#'   \code{Deconvolution subgroups composition} are used).
#' @param cell_groups The cell groups of the training set, as returned by \code{construct_cell_groups()}: its elements
#'   \code{Composition} (features of each cell group) and \code{Weights} (CCA projection parameters \code{xcoef},
#'   \code{train_means}, \code{train_sds} of each cell group) are used. Older unnamed lists, with these as second and
#'   third elements, are also accepted.
#' @param group_names Optional character vector with the names of the cell groups to score. If \code{NULL} (default), all
#'   cell groups are scored. Names that are not cell groups are ignored with a message; an error is raised if none
#'   of them is a cell group.
#' @param deconv_test A data frame or matrix with the deconvolution features of the test set before
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
#' If the cell groups were built with batch correction (\code{batch} in \code{construct_cell_groups()}), the features of
#' the test samples are z-scored on their own means and SDs (as each training cohort was), so \code{deconv_test} should contain a
#' single cohort with several samples: the scores then describe each sample relative to the other test samples.
#' A warning is given when such cell groups are scored with fewer than 10 test samples.
#'
#' @examples
#' data("network.tuto")
#' data("deconv_subgroups.tuto")
#' data("deconv.tuto")
#' cell_groups <- construct_cell_groups(network.tuto, deconv_subgroups.tuto)
#' test_scores <- compute.test.set(deconv_subgroups.tuto, cell_groups,
#'                                 deconv_test = deconv.tuto[1:20, ])
#' test_scores[1:3, 1:3]
#'
compute.test.set = function(dt, cell_groups, group_names = NULL, deconv_test){

  # Features and CCA parameters of each cell group (by name; older unnamed lists by position)
  unnamed = is.list(cell_groups) && is.null(names(cell_groups)) && length(cell_groups) >= 3
  composition = if(is.list(cell_groups) && !is.null(cell_groups[["Composition"]])) cell_groups[["Composition"]] else if(unnamed) cell_groups[[2]] else NULL
  weights = if(is.list(cell_groups) && !is.null(cell_groups[["Weights"]])) cell_groups[["Weights"]] else if(unnamed) cell_groups[[3]] else NULL
  if(!is.list(composition) || !is.list(weights) || length(composition) == 0 || length(weights) != length(composition) || is.null(names(composition))){
    stop("cell_groups must be the output of construct_cell_groups() (elements 'Composition' and 'Weights')")
  }
  if(is.null(group_names)) group_names = names(composition)
  unknown = setdiff(group_names, names(composition))
  if(length(unknown) > 0 && length(unknown) < length(group_names)){
    message("Not cell groups of cell_groups, ignored: ", paste(utils::head(unknown, 10), collapse = ", "), if(length(unknown) > 10) paste0(", ... (", length(unknown), " in total)"))
  }

  ################################################################################Replicate cell subgroups (training features, in the same order)
  deconv_test = replicate_deconvolution_subgroups(dt, deconv_test)

  # Compute composite scores
  idx = which(names(composition) %in% group_names)
  if(length(idx) == 0){
    stop("None of group_names match the cell group names")
  }
  cell_dendrogram = c()
  names = c()
  skipped = c()
  for (i in 1:length(idx)) {
    group = names(composition)[idx[i]]
    cells = composition[[idx[i]]]
    pca_cells = deconv_test[,colnames(deconv_test) %in% cells, drop = F]
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
         "column names of deconv_test are the deconvolution features used for training)")
  }
  if(length(skipped) > 0){
    warning("Cell groups skipped because some of their features are missing or have missing values in the test set (their score would not ",
            "be comparable with the training scores): ",
            paste(utils::head(skipped, 10), collapse = ", "), if(length(skipped) > 10) paste0(", ... (", length(skipped), " in total)"), call. = FALSE)
  }

  colnames(cell_dendrogram) = names

  # Cell groups built with batch: the test samples are their own reference (z-scored on their own means and SDs)
  batch_trained = any(sapply(names, function(group) is.null((if(!is.null(names(weights)) && group %in% names(weights)) weights[[group]] else weights[[match(group, names(composition))]])$train_means)))
  if(batch_trained && nrow(deconv_test) < 10){
    warning("The cell groups were built with batch correction, so the test samples are z-scored on their own means and SDs and their scores depend ",
            "on the other test samples. With fewer than 10 test samples (here ", nrow(deconv_test), ") the scores are unreliable",
            if(nrow(deconv_test) == 1) " (with one sample, every score is 0)", ": give a whole cohort as test set.", call. = FALSE)
  }

  return(data.frame(cell_dendrogram, check.names = FALSE))

}

#' Merge highly correlated TF modules
#'
#' Merges TF modules whose scores (eigengenes) are too similar, repeating until no pair of modules has a correlation
#' above \code{corr_mod}: at each step the most correlated pair above \code{corr_mod} is merged, and the correlations are
#' computed again. The merged module keeps the name and color of the first module of the pair; its score is the mean
#' of the scores of all the original modules it contains, scaled (mean 0, standard deviation 1).
#'
#' @param data A data frame of TF module eigengenes (samples x modules), with columns named \code{ME<color>}.
#' @param colors A character vector with the module color label of each TF.
#' @param corr_mod Numeric. Correlation threshold above which two modules are merged.
#' @param corr_type Character. Correlation used: "pearson" or "spearman" (default).
#'
#' @return A list of two elements:
#' \itemize{
#'   \item \code{[[1]]}: Data frame of merged module eigengenes (samples x modules).
#'   \item \code{[[2]]}: Updated character vector of module color labels after merging.
#' }
#'
#' @keywords internal
mergeModules = function(data, colors, corr_mod, corr_type = "spearman"){
  original = data
  members = stats::setNames(as.list(colnames(data)), colnames(data)) # Original modules contained in each module

  repeat{
    if(ncol(data) < 2) break
    cc = stats::cor(data, method = corr_type, use = "pairwise.complete.obs")
    cc[lower.tri(cc, diag = TRUE)] = NA # Each pair once
    if(all(is.na(cc)) || max(cc, na.rm = TRUE) <= corr_mod) break

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
#'   \item the dendrogram is cut with dynamic tree cutting (\code{dynamicTreeCut::cutreeDynamic()}, \code{method = "tree"},
#'     \code{deepSplit = TRUE}, which splits clusters more finely than \code{FALSE});
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
#' @param cell_dendrograms A list of \code{hclust} objects (or objects convertible to \code{hclust}),
#'   one per TF module, as returned by \code{identify.cell.groups()}.
#' @param min_cluster_size Integer. Minimum number of features per cluster: passed to
#'   \code{dynamicTreeCut::cutreeDynamic()}'s \code{minClusterSize} argument, and used in step 2. Default is 3.
#'
#' @return A list of integer cluster label vectors, one per dendrogram, in the
#'   same order as \code{cell_dendrograms} and named like it, each in the order of the dendrogram labels
#'   (label \code{0} marks unassigned elements).
#'
#' @keywords internal
calculate_dendrogram_cuts = function(cell_dendrograms, min_cluster_size = 3){

  cuts = list()
  no_values = c()
  for (i in 1:length(cell_dendrograms)) {
    hc <- if(inherits(cell_dendrograms[[i]], "hclust")) cell_dendrograms[[i]] else as.hclust(cell_dendrograms[[i]])

    labels <- dynamicTreeCut::cutreeDynamic(
      dendro         = hc,
      method         = "tree",
      deepSplit      = TRUE,
      minClusterSize = min_cluster_size,
      verbose        = 0
    )

    # Each cluster must be one range of correlation values: split clusters into their separated pieces,
    # and join pieces smaller than min_cluster_size to their closest neighbour
    values <- hc$correlation
    if(is.null(values) || length(values) != length(labels)){
      no_values <- c(no_values, if(is.null(names(cell_dendrograms))) i else names(cell_dendrograms)[i])
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
  names(cuts) <- names(cell_dendrograms)

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
#' @param cell_dendrograms A named list of \code{hclust} objects (or objects convertible to \code{hclust}),
#'   one per TF module, as returned by \code{identify.cell.groups()} (the correlations they store are used to order the
#'   leaves; without them, the leaves keep the order of the tree).
#' @param cuts_per_dendrogram A list of integer cluster label vectors, one per dendrogram,
#'   as returned by \code{calculate_dendrogram_cuts()}.
#' @param file_name Optional character. Suffix of the output file
#'   (\code{Results/Dendrogram_color_clusters_<file_name>.pdf}; \code{Results/Dendrogram_color_clusters.pdf} if
#'   \code{NULL}, the default).
#'
#' @return Called for its side effect (saves a PDF); returns \code{NULL} invisibly.
#'
#' @keywords internal
plot_dendrogram_clusters = function(cell_dendrograms, cuts_per_dendrogram, file_name = NULL) {
  plots <- list()

  for (i in 1:length(cell_dendrograms)) {
    hc     <- if (inherits(cell_dendrograms[[i]], "hclust")) cell_dendrograms[[i]] else as.hclust(cell_dendrograms[[i]])
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
      ggplot2::ggtitle(paste("Dendrogram:", names(cell_dendrograms)[i]),
                       subtitle = paste0("Numbers: clusters (cell groups Dendrogram_", names(cell_dendrograms)[i], "_group_<number>)")) +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text   = ggplot2::element_blank(),
                     axis.ticks  = ggplot2::element_blank(),
                     plot.title  = ggplot2::element_text(hjust = 0.5, face = "bold"),
                     plot.subtitle = ggplot2::element_text(hjust = 0.5)) +
      ggplot2::labs(x = NULL, y = "Height")

    plots[[i]] <- p
  }

  out_path <- file.path("Results", paste0("Dendrogram_color_clusters", if (!is.null(file_name)) paste0("_", file_name), ".pdf"))
  pdf(out_path, width = 10, height = 12)
  for (p in plots) print(p)
  dev.off()
  message("Dendrogram plots saved to: ", out_path)

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
#' @param counts_norm A numeric matrix of normalized expression values (genes x samples), with gene symbols as row names.
#' @param module_color Character. Name of the TF module (color label) to analyze.
#' @param hub_tfs A list as returned by \code{identify_hub_TFs()}, where the first
#'   element maps module names to vectors of hub TF gene symbols.
#' @param universe A data frame of TF-target interactions with at minimum columns
#'   \code{source} (TF) and \code{target} (gene).
#' @param min_genes Integer. Minimum number of genes to test (20% most variable targets); with fewer, no enrichment is
#'   computed and a message is given. Default is 5.
#'
#' @return A \code{ReactomePA} enrichResult object restricted to pathways with
#'   adjusted p-value < 0.05, or \code{NULL} if no enrichment is found, if none of the
#'   targets is in the row names of \code{counts_norm}, or if fewer than \code{min_genes} genes would be tested.
#'
#' @keywords internal
module_enrich = function(counts_norm, module_color, hub_tfs, universe, min_genes = 5){
  targets = universe[make.names(universe$source) %in% make.names(hub_tfs[[1]][[module_color]]),] #Extract targets from TFs (syntactic names too, e.g. NKX2.1 for NKX2-1)
  targets = unique(targets$target) #Keep only unique targets from TFs

  targets_genes = counts_norm[rownames(counts_norm)%in%targets, , drop = F] #Extract gene expression from targets
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
  background = suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = rownames(counts_norm), columns = "ENTREZID", keytype = "SYMBOL")) #Change to EntrezID
  unmapped_tested = length(unique(entrz$SYMBOL[is.na(entrz$ENTREZID)]))
  unmapped_universe = length(unique(background$SYMBOL[is.na(background$ENTREZID)]))
  if(unmapped_tested > 0 || unmapped_universe > 0){
    message("Module ", module_color, ": ", unmapped_tested, " of the ", n_tested, " tested genes and ", unmapped_universe, " of the ", nrow(counts_norm),
            " background genes have no Entrez ID and are left out")
  }

  reac <- suppressMessages(ReactomePA::enrichPathway(gene    = unique(stats::na.omit(entrz$ENTREZID)),
                                                     organism     = 'human',
                                                     universe = unique(stats::na.omit(background$ENTREZID)),
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
#' between the deconvolution features of a cell group and the score of its TF module (module eigengene).
#' The score is the weighted sum of the (scaled) features of the cell group, with the weights of the first
#' canonical component: the combination of the features that correlates most with the TF module score (with a single
#' module score, this is the linear regression of the module score on the features, up to a scaling factor). As the sign of CCA weights is arbitrary, the score is oriented so that it is aligned
#' with the TF module: a positive score means that the sample is on the active side of the module for that cell
#' group. The score therefore increases with the features of the cell group when they are positively related
#' to the module, and decreases with them when they are negatively related.
#'
#' @param cell_group A numeric matrix of cell deconvolution features for a cell group (samples x features), with sample
#'   names as row names. The module score of \code{network} is matched to these samples by name; every sample must
#'   be in the network.
#' @param module Character. Name (color) of the TF module the cell group was built from; must match a
#'   column name of the module matrix in \code{network} exactly.
#' @param network Output of compute.WTCNA() (single- or multi-cohort).
#' @param batch Optional vector with the cohort of each sample, in the row order of \code{cell_group}. It is treated as
#'   categorical: the cell group features are z-scored within each cohort (mean 0, SD 1 per cohort; a feature constant
#'   in a cohort is 0 there) before the CCA. The module score is used as it is: \code{compute.WTCNA()} already scales it
#'   within each cohort.
#'
#' @return An unnamed list of two elements:
#' \itemize{
#'   \item \code{[[1]]}: Numeric matrix (samples x 1) with the composite score, i.e. the scaled cell group
#'     features projected onto the first canonical component, oriented to be aligned with the TF module. With
#'     \code{batch}, it has mean 0 within each cohort.
#'   \item \code{[[2]]}: Projection parameters used to score new samples: \code{xcoef} (canonical weights of
#'     the first component), \code{train_means} and \code{train_sds} (column means/SDs used for scaling).
#'     The sign of \code{xcoef} follows the orientation of the score.
#'     When batch correction was applied, only \code{xcoef} is returned: new samples are then z-scored on their own
#'     means and SDs, as each training cohort was.
#' }
#'
#' @details
#' The CCA is reliable only when the number of features of the cell group is well below the number of samples. When
#' it reaches the number of samples, the first canonical correlation is 1 whatever the data, the weights are not unique
#' and the score is decided by numerical noise. Cell groups with many features and few samples (e.g. cross-validation
#' folds) are the most affected.
#'
#' @examples
#' data("network.tuto")
#' data("deconv_subgroups.tuto")
#' cell_groups <- construct_cell_groups(network.tuto, deconv_subgroups.tuto)
#' group <- names(cell_groups$Composition)[1]
#' module <- sub("^Dendrogram_(.*)_group_.*$", "\\1", group) # TF module of the cell group
#' features <- deconv_subgroups.tuto[[1]][, cell_groups$Composition[[group]], drop = FALSE]
#' score <- compute_composite_score(features, module, network.tuto)
#' head(score[[1]])
#'
#' @export
#'
compute_composite_score = function(cell_group, module, network, batch = NULL){

  modules = network[["TFs module matrix"]]
  if(is.null(modules)){
    stop("network must be the output of compute.WTCNA()")
  }
  if(!module %in% colnames(modules)){
    stop("module '", module, "' is not a TF module of network")
  }
  if(!is.null(batch) && length(batch) != nrow(cell_group)){
    stop("batch must have one value per sample of cell_group (", nrow(cell_group), ")")
  }

  # Module score in the sample order of cell_group (matched by name)
  if(is.null(rownames(cell_group))) stop("cell_group must have the sample names as row names")
  idx_modules = match(make.names(rownames(cell_group)), make.names(rownames(modules)))
  if(anyNA(idx_modules)){
    missing_samples = rownames(cell_group)[is.na(idx_modules)]
    stop(length(missing_samples), " samples of cell_group are not in network: ", paste(utils::head(missing_samples, 5), collapse = ", "))
  }
  modules <- modules[idx_modules, , drop = F]

  tf_module_matrix = modules[, colnames(modules) == module, drop = F]

  # ---- Correct for batch if provided ----
  batch_corrected = !is.null(batch) && length(unique(batch)) > 1
  if(batch_corrected){
    batch <- as.factor(batch)  # Categorical
    # Cell group features z-scored within each cohort (mean 0, SD 1 per cohort), as the TF activities and module scores;
    # a feature constant in a cohort is 0 there
    cell_group <- apply(cell_group, 2, function(x) stats::ave(x, batch, FUN = function(v) if(isTRUE(stats::sd(v) > 0)) (v - mean(v)) / stats::sd(v) else v - mean(v)))
  }

  cell_group_matrix = if(batch_corrected) as.matrix(cell_group) else scale(cell_group) # With batch, the features are already z-scored within each cohort

  # CCA analysis: Note 'cell_group_matrix' is already scale
  cca_result <- stats::cancor(cell_group_matrix, tf_module_matrix) # Between the features of the cell group and the TF module score

  # Apply canonical weights to each matrix
  cell_group <- cell_group[, rownames(cca_result$xcoef), drop = F] #Ensure no features have been discard on the way and if yes, subset the matrix (CCA can discard collinear or zero variance features)

  scaled_obj <- if(batch_corrected) as.matrix(cell_group) else scale(as.matrix(cell_group)) ### Why do you scale again? Because the CCA weights are based on the scale matrix
  train_means <- colMeans(cell_group)
  train_sds   <- apply(cell_group, 2, stats::sd)

  ### xcoef[,1] corresponds to the most correlated linear component (we used the scale cell_group_matrix because the coef came from the scale matrix)
  xcoef = cca_result$xcoef[, 1, drop = F]
  weighted_cell_group_matrix <- scaled_obj %*% xcoef

  # The CCA weights are only defined up to their sign (w and -w are equally good solutions, and which one is returned
  # depends e.g. on the sample order): orient the score so that it is aligned with the TF module. It increases with the
  # features of the cell group when they are positively related to the module, and decreases with them when negatively related
  # Positive CellTFusion score = higher TF module activity
  # ---------------------------------------------------------
  module_score <- as.numeric(tf_module_matrix[, 1])

  orientation_cor <- stats::cor(
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
  if(batch_corrected) projection_params = list(xcoef = xcoef) # Batch-corrected: new cohorts are z-scored on their own means and SDs, as each training cohort was

  return(list(weighted_cell_group_matrix, projection_params))

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
#'     \item{train_means}{Named numeric vector of training column means.}
#'     \item{train_sds}{Named numeric vector of training column standard deviations.}
#'   }
#'   \code{train_means} and \code{train_sds} are absent when the training used batch correction: the test data are
#'   then z-scored on their own column means and SDs (a constant feature is 0), as each training cohort was.
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
  if(is.null(train_means)) train_means = colMeans(cell_group) # Batch-corrected training: z-score the test cohort on its own means and SDs
  if(is.null(train_sds)){
    train_sds = apply(cell_group, 2, stats::sd)
    train_sds[is.na(train_sds) | train_sds == 0] = 1 # A constant feature is 0, as in training
  }
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
#' @param padj Logical. If TRUE, the p-values are adjusted with the Benjamini-Hochberg method.
#'
#' @return A list with \code{result} (see the \code{scores.*} tests), \code{sig} (indices of the significant
#'   features) and \code{adjusted_text} (adjusted p-value of each feature for the plot subtitles, empty without
#'   adjustment).
#'
#' @keywords internal
finish_trait_test = function(input, statistic, p, test, statistic_name, pval, padj){
  if (!isTRUE(padj) && !isFALSE(padj)) stop("padj must be TRUE or FALSE")
  scores = input$scores
  features = colnames(input$X)
  p_adjusted = stats::p.adjust(p, method = if (padj) "BH" else "none")
  sig = which(!is.na(p_adjusted) & p_adjusted <= pval)

  statistics = data.frame(feature = features, test = test, statistic_name = statistic_name, statistic = statistic,
                          p_value = p, p_adjusted = p_adjusted, significant = seq_along(features) %in% sig,
                          stringsAsFactors = FALSE)

  result = list()
  result[[1]] = scores[[1]][, sig, drop = FALSE]
  if (length(scores) >= 2) result[[2]] = scores[[2]][sig]
  if (length(scores) >= 3) result[[3]] = scores[[3]][sig]
  if (!is.null(names(scores))) names(result) = names(scores)[seq_along(result)]
  result[["Statistics"]] = statistics

  threshold = paste0(if (padj) "adjusted p-value (BH)" else "p-value", " <= ", pval)
  if (length(sig) == 0) {
    message("No significant feature after ", test, " (", threshold, ")")
  } else {
    message(length(sig), " significant feature(s) after ", test, " (", threshold, "): ",
            paste(utils::head(features[sig], 10), collapse = ", "), if (length(sig) > 10) ", ...")
  }

  # Adjusted p-value for the plot subtitles (empty without adjustment)
  adjusted_text = if (!padj) rep("", length(p)) else
    paste0(" | adjusted p (BH) ", ifelse(p_adjusted < 0.001, "< 0.001", paste0("= ", round(p_adjusted, 4))))

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
#' @param padj Logical; if TRUE, the p-values are corrected for multiple testing across the features with the
#'   Benjamini-Hochberg method, and the adjusted p-values are compared with the threshold. Default FALSE.
#'
#' @return A list with the significant features: its first element is the score matrix restricted to them (no
#'   column if none is significant), followed by the matching parts of the other elements of \code{scores} (e.g.
#'   \code{Composition} and \code{Weights} of \code{construct_cell_groups()}), with their names.
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
scores.ttest <- function(scores, coldata, trait, pval = 0.05, padj = FALSE) {
  input = prepare_trait_test(scores, coldata, trait, binary = TRUE, test = "Welch's t-test")
  X = input$X

  tests = lapply(seq_len(ncol(X)), function(j) tryCatch(stats::t.test(X[, j] ~ input$trait, var.equal = FALSE), error = function(e) NULL))
  res = finish_trait_test(input, sapply(tests, function(t) if (is.null(t)) NA else unname(t$statistic)),
                          sapply(tests, function(t) if (is.null(t)) NA else t$p.value), "Welch's t-test", "t", pval, padj)

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
scores.kruskal.test <- function(scores, coldata, trait, pval = 0.05, padj = FALSE) {
  input = prepare_trait_test(scores, coldata, trait, test = "Kruskal-Wallis test")
  X = input$X

  tests = lapply(seq_len(ncol(X)), function(j) tryCatch(stats::kruskal.test(X[, j], input$trait), error = function(e) NULL))
  res = finish_trait_test(input, sapply(tests, function(t) if (is.null(t)) NA else unname(t$statistic)),
                          sapply(tests, function(t) if (is.null(t)) NA else t$p.value), "Kruskal-Wallis test", "H", pval, padj)

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
scores.wilcox.test <- function(scores, coldata, trait, pval = 0.05, padj = FALSE) {
  input = prepare_trait_test(scores, coldata, trait, binary = TRUE, test = "Wilcoxon test")
  X = input$X

  tests = lapply(seq_len(ncol(X)), function(j) tryCatch(stats::wilcox.test(X[, j] ~ input$trait, exact = FALSE), error = function(e) NULL))
  res = finish_trait_test(input, sapply(tests, function(t) if (is.null(t)) NA else unname(t$statistic)),
                          sapply(tests, function(t) if (is.null(t)) NA else t$p.value), "Wilcoxon rank-sum test", "W", pval, padj)

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
scores.anova.test = function(scores, coldata, trait, pval = 0.05, padj = FALSE){
  input = prepare_trait_test(scores, coldata, trait, test = "ANOVA")
  X = input$X

  tests = lapply(seq_len(ncol(X)), function(j) tryCatch(summary(stats::aov(X[, j] ~ input$trait))[[1]][1, c("F value", "Pr(>F)")], error = function(e) NULL))
  res = finish_trait_test(input, sapply(tests, function(t) if (is.null(t)) NA else t[["F value"]]),
                          sapply(tests, function(t) if (is.null(t)) NA else t[["Pr(>F)"]]), "one-way ANOVA", "F", pval, padj)

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
scores.fisher.test = function(scores, coldata, trait, pval = 0.05, padj = FALSE){
  input = prepare_trait_test(scores, coldata, trait, test = "Fisher's exact test")
  X = input$X

  # Each score split at its median into High/Low
  levels_list = lapply(seq_len(ncol(X)), function(j) factor(ifelse(X[, j] > stats::median(X[, j]), "High", "Low"), levels = c("High", "Low")))
  tests = lapply(levels_list, function(l) tryCatch(stats::fisher.test(table(l, input$trait), workspace = 2e7), error = function(e) NULL))
  res = finish_trait_test(input, sapply(tests, function(t) if (is.null(t) || is.null(t$estimate)) NA else unname(t$estimate)),
                          sapply(tests, function(t) if (is.null(t)) NA else t$p.value), "Fisher's exact test", "odds ratio", pval, padj)

  for (j in res$sig) {
    p_value = res$result$Statistics$p_value[j]
    df = data.frame("score_level" = levels_list[[j]], "Trait" = input$trait)
    p_label <- ifelse(p_value < 0.001, "< 0.001", paste0("= ", round(p_value, 4)))

    # Percentage of High/Low samples within each trait group
    counts <- df %>% dplyr::count(Trait, score_level) %>%
      dplyr::group_by(Trait) %>% dplyr::mutate(proportion = n / sum(n)) %>% dplyr::ungroup()
    n_labels  <- df %>% dplyr::count(Trait) %>%
      dplyr::mutate(label = paste0(as.character(Trait), "\n(n=", n, ")"))
    label_map <- setNames(n_labels$label, as.character(n_labels$Trait))

    pdf(paste0("Results/Fisher_", trait, "_", colnames(X)[j], ".pdf"), width = 8, height = 6)
    print(
      ggplot2::ggplot(counts, ggplot2::aes(x = Trait, y = proportion, fill = score_level)) +
        ggplot2::geom_col(width = 0.7, color = "white") +
        ggplot2::geom_text(ggplot2::aes(label = paste0(round(100 * proportion), "%")),
                           position = ggplot2::position_stack(vjust = 0.5), size = 4) +
        ggplot2::scale_x_discrete(labels = label_map) +
        ggplot2::scale_y_continuous(labels = scales::percent) +
        ggplot2::scale_fill_brewer(palette = "Set2") +
        ggplot2::labs(
          title    = colnames(X)[j],
          subtitle = paste0("Fisher's exact test | p ", p_label, " | ", trait, res$adjusted_text[j]),
          x        = trait,
          y        = "Percentage",
          fill     = "Score level"
        ) +
        ggplot2::theme_classic(base_size = 13) +
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
#'   the score matrix; optional second and third elements (e.g. the output of
#'   \code{construct_cell_groups()}) are subsetted to the significant features and returned with their names.
#' @param coldata A data frame (or tibble) containing clinical or experimental metadata for samples.
#'   Must include the column specified in `trait`. Its rows are matched to the samples of the scores by row name;
#'   without sample names, they must be in the order of the scores. Samples with a missing trait are left out.
#' @param trait Character. The name of the column in `coldata` representing the clinical
#'   or experimental trait to test against (e.g., response, subtype, etc.). It must be categorical.
#' @param test Character. Statistical test to perform. One of:
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
#' @param padj Logical; if TRUE, the p-values are corrected for multiple testing across the features with the
#'   Benjamini-Hochberg method, and the adjusted p-values are compared with the threshold. Default FALSE.
#'
#' @details
#' The function automatically calls the corresponding statistical test function
#' based on the `test` argument:
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
#'   is significant). Optional second and third elements (from the input list, e.g. the composition and weights
#'   of the cell groups) are subsetted accordingly and keep their names. The last element,
#'   \code{Statistics}, is a data frame with the test of every feature (\code{feature}, \code{test},
#'   \code{statistic_name}, \code{statistic}, \code{p_value}, \code{p_adjusted}, \code{significant}).
#'
#' @examples
#' data("network.tuto")
#' data("deconv_subgroups.tuto")
#' data("traitdata.tuto")
#' cell_groups <- construct_cell_groups(network.tuto, deconv_subgroups.tuto)
#' sig <- scores.stat.analysis(cell_groups, traitdata.tuto, trait = "Enrollment.IC",
#'                             test = "kruskal", pval = 0.05, padj = TRUE)
#' head(sig$Statistics)
#' names(sig$Composition)
#'
#' @export
#'
scores.stat.analysis <- function(scores, coldata, trait,
                                 test = c("fisher", "wilcox", "anova", "kruskal", "ttest"),
                                 pval = 0.05, padj = FALSE) {
  test <- match.arg(test)

  result <- switch(test,
                   fisher   = scores.fisher.test(scores, coldata, trait, pval, padj),
                   wilcox   = scores.wilcox.test(scores, coldata, trait, pval, padj),
                   anova    = scores.anova.test(scores, coldata, trait, pval, padj),
                   kruskal  = scores.kruskal.test(scores, coldata, trait, pval, padj),
                   ttest    = scores.ttest(scores, coldata, trait, pval, padj))

  return(result)
}

#' Compute latent factors from cell group scores using NMF
#'
#' Decomposes the cell group scores into non-negative latent factors with Non-negative Matrix Factorization (NMF,
#' \code{RcppML::nmf()}), applied to the positive part of the scores (see Details). The number of factors is chosen
#' by bootstrap NMF if \code{rank} is not given.
#'
#' @param scores Numeric matrix of size samples x cell groups (signed CCA composite scores).
#' @param rank Integer; number of NMF factors. If NULL (default), chosen among \code{ranks} with bootstrap NMF (see
#'   Details). If given, the stability of its factors is still computed (\code{stability}).
#' @param ranks Candidate ranks when \code{rank = NULL}. Default 2:15, limited to the number of cell groups (rank 1 with
#'   a single cell group).
#' @param n_boot Integer; number of bootstrap NMF runs per rank. Default 100.
#' @param recovery_thresh Numeric; minimum recovery frequency (proportion of bootstrap runs in which the factor is
#'   found) for a factor to be robust. Default 0.7.
#' @param similarity_thresh Numeric; minimum cosine similarity for a bootstrap factor to count as the same factor, and
#'   minimum mean similarity for a factor to be robust. Two factors of the same rank more similar than this are
#'   considered redundant. Default 0.8. Both thresholds are practical starting points, not universal cutoffs.
#' @param seed Random seed used for the NMF fits and the bootstrap samples. Default 123. The caller's random number
#'   generator state is restored when the function returns.
#'
#' @param file_name Optional character suffix for the saved patient-mixture plot.
#' @param return Logical. If TRUE (default), saves the patient-mixture barplot to
#'   \code{Results/NMF_patient_mixture_<file_name>.pdf} (\code{Results/NMF_patient_mixture.pdf} if
#'   \code{file_name} is \code{NULL}), and the stability plot to \code{Results/NMF_rank_stability_<file_name>.pdf}:
#'   the lowest recovery, lowest mean similarity and highest similarity between factors of each candidate rank (with
#'   the thresholds and the selected rank; only when the rank is selected), and the recovery and mean similarity of
#'   each factor of the selected rank.
#'
#' @return A named list with:
#' \describe{
#'   \item{Z}{Sample-level NMF factor scores (samples x rank). Non-negative. Each factor (column) sums to 1 over the
#'     samples: the size of each factor is in \code{nmf_model$d} (the NMF is \code{W \%*\% diag(d) \%*\% t(Z)}). Tests
#'     on one factor at a time do not depend on this scale; to compare factors within a sample, multiply each column
#'     by its \code{d}, as done for \code{patient_mixture}.}
#'   \item{W}{Weights of the cell groups per factor (cell groups x rank). Non-negative.}
#'   \item{nmf_input}{The matrix fed to NMF: the positive part of the cell group scores (samples x cell groups).}
#'   \item{nmf_model}{The \code{RcppML::nmf()} model object (includes the scaling vector \code{d}).}
#'   \item{patient_mixture}{Long-format data frame of per-sample factor proportions used for the mixture plot: the
#'     factor scores multiplied by the size of each factor (\code{nmf_model$d}), divided by their sum in each sample.}
#'   \item{stability}{One row per factor of the selected rank: \code{recovery} (proportion of bootstrap runs in which a
#'     matching factor is found, cosine similarity >= \code{similarity_thresh}), \code{mean_similarity} (mean cosine
#'     similarity of the matched bootstrap factors), \code{max_similarity_other} (highest cosine similarity with another
#'     factor of the same rank: a high value suggests redundancy) and \code{robust} (recovery >= \code{recovery_thresh}
#'     and mean similarity >= \code{similarity_thresh}).}
#'   \item{rank_selection}{One row per candidate rank (\code{NULL} if \code{rank} is given): the number of robust
#'     factors, the lowest recovery and mean similarity of its factors, the highest similarity between two of its factors,
#'     and \code{stable_distinct} (all factors robust and none redundant).}
#' }
#'
#' @details
#' Only the positive part of the signed cell group scores is used: \code{max(score, 0)}, how far the sample is on the
#' active side of the TF module of the cell group (negative scores, on the inactive side, become 0). The factors
#' therefore describe the high-activity direction of the cell groups. The scores must not contain missing values.
#'
#' Rank selection by bootstrap NMF (consensus of the factors across resamples of the patients). The profile of a factor
#' is its column of \code{W} (the weights of the cell groups). For each candidate rank:
#' \enumerate{
#'   \item NMF is fitted on all the samples (reference factors) and on \code{n_boot} bootstrap samples of the patients
#'     (sampled with replacement);
#'   \item the factors of each bootstrap run are matched one to one to the reference factors by the cosine similarity of
#'     their profiles (most similar pairs first), as NMF factors can come in any order;
#'   \item each reference factor gets its recovery frequency and mean similarity across runs, and the factors of the rank
#'     are compared with each other to detect redundancy (two factors with almost the same profile).
#' }
#' The selected rank is the largest one whose factors are all robust and none redundant: as many factors as the data
#' support reliably (low ranks capture the dominant patterns and are almost always stable, so the smallest stable rank
#' would favour 2). If no rank qualifies, the rank with the most robust factors is used (the smallest in case of ties),
#' with a message; check \code{rank_selection} and whether the main factors persist at the neighbouring ranks.
#'
#' @examples
#' data("network.tuto")
#' data("deconv_subgroups.tuto")
#' cell_groups <- construct_cell_groups(network.tuto, deconv_subgroups.tuto)
#' latent <- compute.latent_factors(cell_groups$Cell_groups, n_boot = 20, return = FALSE)
#' head(latent$Z)
#' latent$nmf_model$d     # size of each factor
#' latent$rank_selection  # stability and redundancy of each candidate rank
#' latent$stability       # stability of the factors of the selected rank
#'
compute.latent_factors <- function(scores, rank = NULL, ranks = 2:15, n_boot = 100, recovery_thresh = 0.7, similarity_thresh = 0.8,
                                   seed = 123, file_name = NULL, return = TRUE) {

  if (!requireNamespace("RcppML", quietly = TRUE))
    stop("Install RcppML: install.packages('RcppML')")

  withr::local_preserve_seed()
  set.seed(seed)
  scores <- as.matrix(scores)
  if (anyNA(scores)) stop("scores has missing values: NMF needs complete cell group scores")

  # -- Positive part of the scores: high-activity side of the TF module of each cell group (negative scores set to 0)
  nmf_input <- pmax(scores, 0)

  A <- t(nmf_input)   # features x samples (RcppML convention)

  # -- Rank selection: bootstrap NMF -----------------------------------------
  # Profile of a factor = its column of W (weights of the cell groups)

  # Cosine similarity between the profiles (columns) of a and b
  cosine <- function(a, b) {
    a <- sweep(a, 2, sqrt(colSums(a^2)), "/")
    b <- sweep(b, 2, sqrt(colSums(b^2)), "/")
    sim <- crossprod(a, b)
    sim[!is.finite(sim)] <- 0 # factor with all weights 0
    sim
  }

  # One-to-one matching of reference factors (rows) to bootstrap factors (columns), most similar pairs first:
  # similarity of the bootstrap factor matched to each reference factor
  match_factors <- function(sim) {
    matched <- rep(0, nrow(sim))
    while (any(!is.na(sim))) {
      ij <- which(sim == max(sim, na.rm = TRUE), arr.ind = TRUE)[1, ]
      matched[ij[1]] <- sim[ij[1], ij[2]]
      sim[ij[1], ] <- NA
      sim[, ij[2]] <- NA
    }
    matched
  }

  # Stability of the factors of rank k: reference fit on all the samples, n_boot fits on bootstrap samples of the patients
  rank_stability <- function(k) {
    ref <- RcppML::nmf(A, k = k, seed = seed, verbose = FALSE)
    sims <- matrix(unlist(lapply(seq_len(n_boot), function(b) {
      boot <- sample(ncol(A), replace = TRUE)
      fit <- RcppML::nmf(A[, boot, drop = FALSE], k = k, seed = seed + b, verbose = FALSE)
      match_factors(cosine(ref$w, fit$w))
    })), ncol = k, byrow = TRUE) # n_boot x k
    within <- cosine(ref$w, ref$w)
    diag(within) <- NA
    data.frame(factor = paste0("Factor", seq_len(k)),
               recovery = colMeans(sims >= similarity_thresh),
               mean_similarity = colMeans(sims),
               max_similarity_other = if (k > 1) apply(within, 1, max, na.rm = TRUE) else 0,
               stringsAsFactors = FALSE)
  }

  rank_selection <- NULL
  if (is.null(rank) && ncol(scores) == 1) {
    rank <- 1L
    message("A single cell group: using rank = 1")
  }
  if (is.null(rank)) {
    ranks <- ranks[ranks >= 1 & ranks <= ncol(scores)] # No more factors than cell groups
    if (length(ranks) == 0) stop("No candidate rank between 1 and the number of cell groups (", ncol(scores), ")")
    message("Selecting the rank with bootstrap NMF (", n_boot, " runs per rank, ranks ", paste(ranks, collapse = ", "), ")...")
    stab <- lapply(ranks, rank_stability)
    robust <- lapply(stab, function(s) s$recovery >= recovery_thresh & s$mean_similarity >= similarity_thresh)
    rank_selection <- data.frame(
      rank = ranks,
      robust_factors = vapply(robust, sum, numeric(1)),
      min_recovery = vapply(stab, function(s) min(s$recovery), numeric(1)),
      min_mean_similarity = vapply(stab, function(s) min(s$mean_similarity), numeric(1)),
      max_similarity_between_factors = vapply(stab, function(s) max(s$max_similarity_other), numeric(1))
    )
    rank_selection$stable_distinct <- vapply(robust, all, logical(1)) & rank_selection$max_similarity_between_factors < similarity_thresh
    if (any(rank_selection$stable_distinct)) {
      rank <- max(rank_selection$rank[rank_selection$stable_distinct])
      message("Selected rank: ", rank, " (largest rank with all factors robust and none redundant)")
    } else {
      rank <- rank_selection$rank[which.max(rank_selection$robust_factors)]
      message("No rank has all its factors robust and none redundant: rank ", rank, " has the most robust factors and is used. ",
              "Check rank_selection and whether its main factors persist at the neighbouring ranks")
    }
    stability <- stab[[match(rank, ranks)]]
  } else {
    stability <- rank_stability(rank)
  }
  stability$robust <- stability$recovery >= recovery_thresh & stability$mean_similarity >= similarity_thresh
  if (!all(stability$robust)) {
    message("Factors not robust across bootstrap runs: ", paste(stability$factor[!stability$robust], collapse = ", "))
  }

  # -- NMF fit ---------------------------------------------------------------
  message("Running RcppML NMF with rank = ", rank, "...")
  nmf_res <- RcppML::nmf(A, k = rank, seed = seed, verbose = FALSE)
  W <- as.matrix(nmf_res$w)    # features x rank
  Z <- t(as.matrix(nmf_res$h)) # samples x rank

  rownames(Z) <- rownames(scores)
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

    # -- Stability plot: rank selection (if done) and robustness of the factors of the selected rank --
    plots <- list()
    if (!is.null(rank_selection)) {
      metrics <- c("Lowest recovery" = "min_recovery", "Lowest mean similarity" = "min_mean_similarity",
                   "Highest similarity between factors" = "max_similarity_between_factors")
      rank_long <- do.call(rbind, lapply(names(metrics), function(m) data.frame(
        rank = rank_selection$rank, metric = m, value = rank_selection[[metrics[[m]]]],
        stable_distinct = ifelse(rank_selection$stable_distinct, "yes", "no"))))
      rank_long$metric <- factor(rank_long$metric, levels = names(metrics))
      thresholds <- data.frame(metric = factor(names(metrics), levels = names(metrics)),
                               value = c(recovery_thresh, similarity_thresh, similarity_thresh))
      plots$ranks <- ggplot2::ggplot(rank_long, ggplot2::aes(x = rank, y = value)) +
        ggplot2::geom_hline(data = thresholds, ggplot2::aes(yintercept = value), linetype = "dashed", colour = "grey40") +
        ggplot2::geom_vline(xintercept = rank, colour = "#D55E00") +
        ggplot2::geom_line(colour = "grey60") +
        ggplot2::geom_point(ggplot2::aes(fill = stable_distinct), shape = 21, size = 2.5) +
        ggplot2::scale_fill_manual(values = c(yes = "#0072B2", no = "white"), name = "All factors robust\nand none redundant") +
        ggplot2::facet_wrap(~ metric, nrow = 1) +
        ggplot2::scale_x_continuous(breaks = rank_selection$rank) +
        ggplot2::coord_cartesian(ylim = c(0, 1)) +
        ggplot2::labs(title = "Bootstrap NMF stability by rank",
                      subtitle = paste0("Dashed lines: thresholds (redundancy above the line in the third panel). Selected rank: ", rank, " (vertical line)"),
                      x = "Rank (number of factors)", y = NULL) +
        ggplot2::theme_bw(base_size = 10)
    }
    factor_long <- rbind(
      data.frame(factor = stability$factor, metric = "Recovery frequency", value = stability$recovery, threshold = recovery_thresh),
      data.frame(factor = stability$factor, metric = "Mean similarity", value = stability$mean_similarity, threshold = similarity_thresh))
    factor_long$metric <- factor(factor_long$metric, levels = c("Recovery frequency", "Mean similarity"))
    factor_long$factor <- factor(factor_long$factor, levels = stability$factor)
    factor_long$robust <- ifelse(stability$robust[match(factor_long$factor, stability$factor)], "yes", "no")
    plots$factors <- ggplot2::ggplot(factor_long, ggplot2::aes(x = factor, y = value, fill = robust)) +
      ggplot2::geom_col(width = 0.6) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = threshold), linetype = "dashed", colour = "grey40") +
      ggplot2::scale_fill_manual(values = c(yes = "#0072B2", no = "grey70"), name = "Robust") +
      ggplot2::facet_wrap(~ metric, nrow = 1) +
      ggplot2::coord_cartesian(ylim = c(0, 1)) +
      ggplot2::labs(title = paste0("Robustness of the ", rank, " factors across ", n_boot, " bootstrap runs"),
                    subtitle = "Dashed lines: thresholds", x = NULL, y = NULL) +
      ggplot2::theme_bw(base_size = 10)
    pdf(paste0("Results/NMF_rank_stability", if (!is.null(file_name)) paste0("_", file_name), ".pdf"), width = 11, height = 4)
    for (pl in plots) print(pl)
    dev.off()
  }

  return(list(
    Z         = Z,
    W         = W,
    nmf_input = nmf_input,
    nmf_model = nmf_res,
    patient_mixture = patient_mixture,
    stability = stability,
    rank_selection = rank_selection
  ))
}

#' Extract top-contributing features from NMF latent factors
#'
#' For each NMF factor (column of the W basis matrix), selects the \code{top_n} cell groups with the highest weight
#' (cell groups with zero weight are never selected).
#'
#' @param latent_factors A list returned by \code{compute.latent_factors()},
#'   containing at minimum a matrix \code{W} (features x factors) of NMF basis weights.
#' @param top_n Integer. Number of cell groups selected per factor. Default 5.
#'
#' @return A named list (one element per factor) of named numeric vectors, where
#'   names are cell group names and values are their NMF basis weights, sorted descending.
#'
#' @keywords internal
extract_contributing_features <- function(latent_factors, top_n = 5) {

  W <- latent_factors$W

  top_contributors <- list()

  for (factor_name in colnames(W)) {

    w <- W[, factor_name]
    w <- w[!is.na(w) & w > 0]

    if (length(w) == 0) {
      message("No features for ", factor_name)
      next
    }

    # top_n cell groups by weight, sorted descending
    top_contributors[[factor_name]] <- utils::head(sort(w, decreasing = TRUE), top_n)
  }

  return(top_contributors)
}

#' Identify cell niches from NMF latent factors
#'
#' For each NMF factor, selects the cell groups that contribute most to it (the \code{top_n} highest weights of the
#' factor in \code{W}, see \code{extract_contributing_features()}) and describes them by their deconvolution features
#' and CCA coefficients. Deconvolution features are kept as they are (subgroups are not expanded into their members).
#'
#' With \code{return = TRUE}, a tree network of each factor is saved, one page per factor: the factor on the left,
#' its top cell groups in the middle (filled with the colour of their TF module) and the deconvolution features of each
#' cell group listed from top to bottom on the right, next to their cell group. The edge between a cell
#' group and a feature is the CCA coefficient of the feature in the cell group score (\code{xcoef} in \code{Weights},
#' on the z-scored features): red if positive (the score increases with the feature), blue if negative, and wider for
#' a larger absolute coefficient (coefficients close to zero are hairlines).
#'
#' The factors describe the high-activity side of the TF modules of their cell groups (the NMF input is the positive
#' part of the cell group scores). A red edge therefore means that the feature is higher in samples where the cell
#' group is active, a blue edge that it is lower. The coefficients of features correlated with each other (e.g.
#' several estimates of the same cell type) are estimated jointly, so a single coefficient can be small or have the
#' opposite sign of the others of its cell group even when the feature itself is related to the TF module like them:
#' read the sign of a single edge together with the rest of its cell group.
#'
#' @param latent_factors A list returned by \code{compute.latent_factors()} (its \code{W} is used).
#' @param cell_groups The cell groups, as returned by \code{construct_cell_groups()} (its \code{Weights} are used).
#' @param top_n Integer. Number of cell groups selected per factor. Default 5.
#' @param return Logical. If \code{TRUE} (default), saves the networks in
#'   \code{Results/Cells_niches_<file_name>.pdf} (\code{Results/Cells_niches.pdf} without \code{file_name}), one page
#'   per factor, and a heatmap of the NMF weights (\code{W}) of all the cell groups (rows: cell groups, ordered by
#'   the factor where their weight is highest and by decreasing weight; columns: factors; outlined: the top cell
#'   groups of each factor) in
#'   \code{Results/Cells_niches_heatmap_<file_name>.pdf}.
#' @param file_name Character. Suffix appended to the output file name.
#'
#' @return A data frame with one row per factor, selected cell group and feature: \code{factor}, \code{cell_group},
#'   \code{module} (TF module of the cell group), \code{weight} (NMF weight of the cell group in the factor),
#'   \code{feature} and \code{coefficient} (CCA coefficient of the feature in the cell group score). Rows are ordered by
#'   factor and by decreasing weight.
#'
#' @examples
#' data("network.tuto")
#' data("deconv_subgroups.tuto")
#' cell_groups <- construct_cell_groups(network.tuto, deconv_subgroups.tuto, return = FALSE)
#' latent <- compute.latent_factors(cell_groups$Cell_groups, return = FALSE)
#' niches <- compute_cells_niches(latent, cell_groups, return = FALSE)
#' head(niches)
#'
#' @export
compute_cells_niches <- function(latent_factors, cell_groups, top_n = 5,
                                 return = TRUE, file_name = NULL) {

  if (!is.list(cell_groups) || is.null(cell_groups[["Weights"]])) {
    stop("cell_groups must be the output of construct_cell_groups()")
  }

  # -- top contributing cell groups per factor -------------------------------
  # named list: factor -> named numeric vector (names = cell groups, values = NMF weights, all > 0)
  top_features_list <- extract_contributing_features(latent_factors, top_n = top_n)

  # -- features of each selected cell group and their CCA coefficients ------
  niches <- do.call(rbind, lapply(names(top_features_list), function(factor_name) {
    top <- top_features_list[[factor_name]]
    do.call(rbind, lapply(names(top), function(g) {
      xcoef <- cell_groups[["Weights"]][[g]]$xcoef
      if (is.null(xcoef)) stop("Cell group ", g, " of the latent factors is not in cell_groups")
      data.frame(factor = factor_name, cell_group = g, module = sub("^Dendrogram_(.*)_group_[0-9]+$", "\\1", g),
                 weight = unname(top[g]), feature = rownames(xcoef), coefficient = unname(xcoef[, 1]),
                 stringsAsFactors = FALSE)
    }))
  }))
  if (is.null(niches)) {
    niches <- data.frame(factor = character(0), cell_group = character(0), module = character(0), weight = numeric(0),
                         feature = character(0), coefficient = numeric(0))
  }
  rownames(niches) <- NULL

  # -- tree network per factor ------------------------------------------------
  if (return && nrow(niches) > 0) {
    if (!dir.exists("Results")) dir.create("Results", recursive = TRUE)
    modules <- unique(niches$module)
    module_cols <- stats::setNames(ifelse(modules %in% grDevices::colors(), modules, "grey60"), modules)

    # Features listed top to bottom (one row each, with a gap of 0.6 rows between cell groups), names written
    # horizontally; the page is sized for the factor with most features and the longest feature name
    # (about 0.6 x text size mm per character)
    text_size <- 3.2
    n_rows <- max(tapply(niches$cell_group, niches$factor, function(g) length(g) + 0.6 * (length(unique(g)) - 1)))
    label_in <- 0.6 * text_size * max(nchar(niches$feature)) / 25.4 + 0.2
    pdf(paste0("Results/Cells_niches", if (!is.null(file_name)) paste0("_", file_name), ".pdf"),
        width = 4 + label_in, height = 2.6 + 0.22 * n_rows)
    for (factor_name in unique(niches$factor)) {
      e <- niches[niches$factor == factor_name, ]
      n <- nrow(e)
      # Factor at x = 0, cell groups at x = 1 (at the mean height of their features), features at x = 2
      e$y <- -cumsum(c(0, ifelse(e$cell_group[-1] != e$cell_group[-n], 1.6, 1)))
      g_y <- tapply(e$y, e$cell_group, mean)[unique(e$cell_group)]
      nodes <- data.frame(cell_group = names(g_y), yg = as.numeric(g_y))
      nodes$module <- e$module[match(nodes$cell_group, e$cell_group)]
      e$yg <- nodes$yg[match(e$cell_group, nodes$cell_group)]
      y0 <- mean(nodes$yg)
      e$sign <- factor(ifelse(e$coefficient > 0, "positive", "negative"), levels = c("positive", "negative"))

      p <- ggplot2::ggplot() +
        ggplot2::geom_segment(data = nodes, ggplot2::aes(x = 0, y = y0, xend = 1, yend = yg), colour = "grey75", linewidth = 0.4) +
        ggplot2::geom_segment(data = e, ggplot2::aes(x = 1, y = yg, xend = 2, yend = y, colour = sign, linewidth = abs(coefficient))) +
        ggplot2::scale_colour_manual(values = c(positive = "#B2182B", negative = "#2166AC"), drop = FALSE, name = "Coefficient",
                                     guide = ggplot2::guide_legend(order = 1)) +
        # Width proportional to |coefficient| from 0, so coefficients close to zero are hairlines
        ggplot2::scale_linewidth_continuous(range = c(0.05, 2.5), limits = c(0, max(abs(e$coefficient))), name = "|Coefficient|",
                                            guide = ggplot2::guide_legend(order = 2)) +
        ggplot2::geom_point(ggplot2::aes(x = 0, y = y0), size = 14, colour = "grey30") +
        ggplot2::geom_text(ggplot2::aes(x = 0, y = y0, label = factor_name), size = 3, colour = "white") +
        ggplot2::geom_point(data = nodes, ggplot2::aes(x = 1, y = yg, fill = module), size = 6, shape = 21, colour = "grey30") +
        ggplot2::scale_fill_manual(values = module_cols, name = "TF module", guide = ggplot2::guide_legend(order = 3, override.aes = list(size = 4))) +
        ggplot2::geom_point(data = e, ggplot2::aes(x = 2, y = y), size = 0.9, colour = "grey40") +
        ggplot2::geom_text(data = e, ggplot2::aes(x = 2.04, y = y, label = feature), hjust = 0, size = text_size) +
        ggplot2::scale_x_continuous(limits = c(-0.25, 2.05), expand = c(0, 0)) +
        ggplot2::scale_y_continuous(limits = c(-n_rows - 0.5, 0.5), expand = c(0, 0)) +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::labs(title = paste0(factor_name, ": top ", nrow(nodes), " cell groups"),
                      subtitle = "Middle: cell groups (colour: TF module). Right: their features, linked by their CCA coefficient") +
        ggplot2::theme_void(base_size = 10) +
        ggplot2::theme(legend.position = "bottom", legend.box = "vertical", legend.justification = "left",
                       plot.margin = ggplot2::margin(10, label_in * 72.27 + 10, 10, 10))
      print(p)
    }
    dev.off()

    # -- heatmap of the NMF weights of all the cell groups in all factors ------
    W <- latent_factors$W
    # Cell groups ordered by the factor where their weight is highest, then by decreasing weight
    best <- max.col(W, ties.method = "first")
    sel_groups <- rownames(W)[order(best, -W[cbind(seq_len(nrow(W)), best)])]
    heat <- data.frame(cell_group = rep(sel_groups, ncol(W)), factor = rep(colnames(W), each = length(sel_groups)),
                       weight = as.vector(W[sel_groups, , drop = FALSE]), stringsAsFactors = FALSE)
    heat$cell_group <- factor(heat$cell_group, levels = rev(sel_groups))
    heat$factor <- factor(heat$factor, levels = colnames(W))
    heat$top <- paste(heat$factor, heat$cell_group) %in% paste(niches$factor, niches$cell_group)
    p_heat <- ggplot2::ggplot(heat, ggplot2::aes(x = factor, y = cell_group, fill = weight)) +
      ggplot2::geom_tile(colour = "grey90") +
      ggplot2::geom_tile(data = heat[heat$top, ], colour = "black", linewidth = 0.6, fill = NA) +
      ggplot2::scale_fill_gradient(low = "white", high = "#B2182B", limits = c(0, NA), name = "NMF weight") +
      ggplot2::labs(title = "NMF weights of the cell groups", subtitle = paste0("Outlined: top ", top_n, " cell groups of each factor"),
                    x = NULL, y = NULL) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text.y = ggplot2::element_text(size = 6))
    ggplot2::ggsave(paste0("Results/Cells_niches_heatmap", if (!is.null(file_name)) paste0("_", file_name), ".pdf"), p_heat,
                    width = 3.5 + 0.5 * ncol(W) + 0.06 * max(nchar(sel_groups)), height = 1.5 + 0.11 * length(sel_groups),
                    limitsize = FALSE)
  }

  return(niches)
}


#' Project cell group scores onto trained NMF latent factors
#'
#' Keeps the positive part of the cell group scores of new samples (as in training) and projects them
#' onto the latent space learned during training using the Moore-Penrose
#' pseudoinverse of the scaled basis matrix \code{W \%*\% diag(d)}. Negative values of the projection (samples with
#' little of a factor; about 30\% of the values on the tutorial data) are set to 0: on the tutorial data, the result
#' correlates 0.998 with an exact non-negative least-squares projection, which reproduces the training \code{Z}. Cell
#' groups missing in \code{scores} are tolerated with a warning naming them: the projection uses the other cell groups
#' (on the tutorial data, with 24 of 171 cell groups missing, the factors still correlate 0.995 or more with the full
#' projection).
#'
#' @param latent_factors A list returned by \code{compute.latent_factors()},
#'   containing at minimum \code{W} (features x rank) and the RcppML model
#'   object with the \code{nmf_model$d} scaling vector.
#' @param scores A samples x cell groups matrix of signed CCA composite
#'   scores for the test cohort (for new samples, computed with \code{compute.test.set()}). Column names must match
#'   the cell groups used in training.
#'
#' @return A non-negative matrix of size samples x rank containing the
#'   projected NMF factor scores for the test cohort, on the scale of \code{Z} (projecting the training samples gives
#'   their \code{Z}, up to the clamping of negative values).
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
project_factors <- function(latent_factors, scores) {

  scores <- as.matrix(scores)

  # -- Step 1: positive part of the scores, as in training --------------------
  nmf_input_test <- pmax(scores, 0)   # samples x cell groups

  # -- Step 2: extract W and absorb d scaling ---------------------------------
  # RcppML stores: A ~= W %*% diag(d) %*% H
  # Must absorb d into W before pseudoinverse: W_scaled = W %*% diag(d)
  W <- latent_factors$W                          # features x rank
  d <- latent_factors$nmf_model$d                # length = rank

  W_scaled <- W %*% diag(d, nrow = length(d))  # features x rank, d absorbed

  # -- Step 3: restrict to features shared between training and test ----------
  common_features <- intersect(rownames(W_scaled), colnames(nmf_input_test))

  if (length(common_features) == 0)
    stop("No overlapping features between training NMF input and test scores. ",
         "Check that column names of scores match training cell group names.")

  missing_groups <- setdiff(rownames(W_scaled), common_features)
  if (length(missing_groups) > 0)
    warning(length(missing_groups), " cell group(s) of the training missing in the test scores, the projection uses the others: ",
            paste(utils::head(missing_groups, 10), collapse = ", "), if (length(missing_groups) > 10) ", ...", call. = FALSE)

  W_sub <- W_scaled[common_features, , drop = FALSE]   # common x rank
  Y_sub <- t(nmf_input_test[, common_features,
                              drop = FALSE])            # common x samples

  # -- Step 4: project via pseudoinverse -------------------------------------
  # H_test = pinv(W_sub) %*% Y_sub  ->  Z_test = t(H_test)
  Z_test <- t(MASS::ginv(W_sub) %*% Y_sub)             # samples x rank

  # The least-squares solution can be negative for factors a sample has little of: set to 0 (non-negative scores)
  Z_test <- pmax(Z_test, 0)

  rownames(Z_test) <- rownames(scores)
  colnames(Z_test) <- paste0("Factor", seq_len(ncol(Z_test)))

  return(Z_test)
}

#' Prepare CellTFusion cross-validation folds for pipeML
#'
#' Fold construction function for \code{pipeML::compute_features.training.ML()} (argument
#' \code{fold_construction_fun}) that computes \code{CellTFusion} cell group scores within each
#' cross-validation fold, without information leakage: cell groups are learned with \code{CellTFusion()} on the
#' training samples of each fold only, and the test samples of the fold are scored on them with
#' \code{compute.test.set()}. The features of the models are the cell group scores; the latent factors are a
#' downstream characterization: they are not computed in the folds, and are computed in final mode (unless
#' \code{compute_factors = FALSE} is passed) for the characterization of the final model; new samples can be projected
#' onto them with \code{project_factors()}. It works for classification and survival tasks.
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
#' @param raw_counts Optional count matrix (genes x samples; TPM values if \code{normalized = FALSE} is passed
#'   to \code{CellTFusion()}) with the samples of \code{data}. Its columns are matched to the samples of
#'   \code{deconv} by name (an error is raised if they are not the same samples); without column names, they must be
#'   in the same order as \code{data}.
#'   \code{pipeML} converts the column names of \code{features_train} with \code{make.names()} (e.g.
#'   \code{"HLA-A"} becomes \code{"HLA.A"}), so gene symbols in \code{data} may no longer match TF regulons and
#'   gene sets; passing \code{raw_counts} keeps the original gene symbols. If \code{NULL} (default),
#'   \code{t(data)} is used.
#' @param batch Optional vector with the cohort of each sample, in the same order as \code{data} (e.g.
#'   \code{coldata$Cohort}). If given, \code{CellTFusion()} runs in multi-cohort mode (see \code{CellTFusion()})
#'   and the test samples of each cohort are scored separately, z-scored on their own means and SDs. Default
#'   \code{NULL} (single cohort).
#' @param do_parallel Logical; whether to compute the folds in parallel. Default is FALSE (one after the other).
#' @param workers (Optional) Number of parallel workers when \code{do_parallel = TRUE}. If \code{NULL} (default),
#'   the number of available cores minus one is used.
#' @param return Logical; in final mode, whether \code{CellTFusion()} saves its outputs and plots in
#'   \code{Results/}. Outputs are never saved for the folds. Default TRUE.
#' @param verbose Logical; whether \code{CellTFusion()} prints its progress messages. Default TRUE.
#' @param ... Other arguments passed to \code{CellTFusion()} (e.g. \code{normalized}, \code{TF.collection},
#'   \code{min_targets_size}, \code{min_module_size}, \code{corr_mod}, \code{corr}, \code{pval}, \code{cancer_type}).
#'   \code{dt}, \code{tfs} and \code{pathways} are not allowed, as they must be computed within each fold.
#'
#' @return
#' \itemize{
#'   \item In fold mode (\code{bestune = NULL}): each fold is saved to \code{Results/fold_<fold name>.rds}, which
#'     \code{pipeML} reads back, as a list with \code{train_data} (cell group scores of the training samples plus
#'     the outcome columns), \code{test_data} (cell group scores of the test samples, plus \code{time} and
#'     \code{event} for survival), \code{obs_test} (observed outcome of the test samples), \code{rowIndex} (rows
#'     of the test samples) and \code{fold_name}. The list of folds is returned invisibly. If \code{CellTFusion()}
#'     fails in a fold (e.g. no cell groups are found on its training samples, or a cell group cannot be scored in
#'     the test samples), that fold is skipped with a
#'     warning and is not saved; an error is raised only if all the folds fail.
#'   \item In final mode (\code{bestune} provided): a list with (1) the cell group scores of all the samples plus
#'     the outcome columns, (2) the \code{CellTFusion()} result on all the samples, available as
#'     \code{Custom_output} in the \code{pipeML} result and used to compute the features of new samples with
#'     \code{compute.test.set()}, and (3) \code{bestune}.
#' }
#'
#' @details
#' Pass the fixed arguments (\code{deconv}, and optionally \code{raw_counts}, \code{batch}
#' and any \code{CellTFusion()} argument) through \code{fold_construction_args_fixed}. Tunable
#' arguments (\code{fold_construction_args_tunable}) are not supported. With \code{batch}, a warning names the folds and cohorts
#' with fewer than 10 test samples (their test samples are z-scored on their own means and SDs, which is unreliable with
#' few samples; folds that hold out whole cohorts avoid it, whereas \code{pipeML}'s \code{LODO = TRUE} stratifies the
#' folds by cohort, so every fold tests part of each cohort). With \code{do_parallel = TRUE}, the first fold is
#' computed before the others start in parallel, so that the TF and pathway collections cached in
#' \code{Results/} are downloaded only once; the package must be installed for the parallel workers.
#'
#' @examples
#' \dontrun{
#' raw_counts <- CellTFusion::raw.counts.tuto
#' traitdata  <- CellTFusion::traitdata.tuto
#' deconv     <- CellTFusion::deconv.tuto  # deconvolution of all the samples
#'
#' # Classification
#' ml_res <- pipeML::compute_features.training.ML(
#'   features_train = t(raw_counts),
#'   task_type      = "classification",
#'   target_var     = traitdata$Best.Confirmed.Overall.Response,
#'   trait.positive = "PD",
#'   metric         = "AUROC",
#'   k_folds        = 5,
#'   n_rep          = 1,
#'   fold_construction_fun        = prepare_celltfusion_folds,
#'   fold_construction_args_fixed = list(deconv = deconv, raw_counts = raw_counts)
#' )
#'
#' # Features of new samples: scored on the cell groups learned on all the training samples
#' final <- ml_res$Custom_output
#' features_test <- compute.test.set(final$Processed_deconvolution, final$Cell_groups,
#'                                   deconv_test = deconv_test)
#'
#' # Survival (time and event of each sample)
#' surv_res <- pipeML::compute_features.training.ML(
#'   features_train = t(raw_counts),
#'   task_type      = "survival",
#'   time_var       = traitdata$os,
#'   event_var      = traitdata$censOS,
#'   k_folds        = 5,
#'   n_rep          = 1,
#'   fold_construction_fun        = prepare_celltfusion_folds,
#'   fold_construction_args_fixed = list(deconv = deconv, raw_counts = raw_counts)
#' )
#' }
#'
#' @export
prepare_celltfusion_folds <- function(data, folds = NULL, bestune = NULL, deconv, raw_counts = NULL,
                                      batch = NULL, do_parallel = FALSE, workers = NULL,
                                      return = TRUE, verbose = TRUE, ...) {

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
  counts <- if (is.null(raw_counts)) t(data[, setdiff(colnames(data), outcome_cols), drop = FALSE]) else raw_counts
  if (is.data.frame(batch) && ncol(batch) == 1) batch <- batch[[1]]  # one-column data frame or tibble
  if (ncol(counts) != nrow(data) || nrow(deconv) != nrow(data) || (!is.null(batch) && (is.data.frame(batch) || length(batch) != nrow(data)))) {
    stop("raw_counts (columns), deconv (rows) and batch must have the same samples, in the same order, as data")
  }
  if (!is.null(batch)) batch <- as.character(batch)
  if (is.null(rownames(deconv))) rownames(deconv) <- rownames(data)
  if (!is.null(raw_counts) && !is.null(colnames(counts))) {
    # Counts matched to the deconvolution samples by name (not by position)
    idx <- match(make.names(rownames(deconv)), make.names(colnames(counts)))
    if (anyNA(idx)) {
      stop("raw_counts and deconv do not have the same samples: ", sum(is.na(idx)), " samples of deconv are not columns of raw_counts (e.g. ",
           paste(utils::head(rownames(deconv)[is.na(idx)], 3), collapse = ", "), ")")
    }
    counts <- counts[, idx, drop = FALSE]
  }
  colnames(counts) <- rownames(deconv)

  # CellTFusion on a subset of samples: cell groups are learned on these samples only. The latent factors are not
  # used as features: they are skipped in the folds, and computed in final mode unless compute_factors = FALSE
  fit_celltfusion <- function(idx, save, compute_factors) {
    args <- list(...)
    args$compute_factors <- compute_factors
    do.call(CellTFusion, c(list(raw_counts = counts[, idx, drop = FALSE], deconv = deconv[idx, , drop = FALSE],
                                batch = if (!is.null(batch)) batch[idx], return = save, verbose = verbose), args))
  }

  # Cell group scores of test samples projected onto the cell groups of a trained CellTFusion result
  # (with batch, each cohort is projected separately, z-scored on its own means and SDs)
  project_samples <- function(structure, idx) {
    cohorts <- if (!is.null(batch)) batch[idx] else rep("all", length(idx))
    projected <- lapply(split(idx, cohorts), function(i) withCallingHandlers(
      compute.test.set(structure$Processed_deconvolution, structure$Cell_groups, deconv_test = deconv[i, , drop = FALSE]),
      warning = function(w) if (grepl("built with batch correction", conditionMessage(w), fixed = TRUE)) invokeRestart("muffleWarning"))) # Already reported per fold and cohort
    do.call(rbind, unname(projected))[rownames(deconv)[idx], , drop = FALSE]
  }

  # -----------------------------
  # CASE 1: bestune provided - CellTFusion on all the training samples
  # -----------------------------
  if (!is.null(bestune)) {
    structure <- fit_celltfusion(seq_len(nrow(data)), save = return, compute_factors = !isFALSE(list(...)$compute_factors))
    train_data_final <- cbind(data.frame(structure$Cell_groups$Cell_groups), outcome)
    return(list(train_data_final, structure, bestune))
  }

  # -----------------------------
  # CASE 2: bestune NOT provided - compute folds
  # -----------------------------
  if (is.null(folds)) stop("Provide 'folds' (a list of training row indices) or 'bestune'.")
  fold_names <- if (is.null(names(folds))) paste0("Fold", seq_along(folds)) else names(folds)

  # With batch, the test samples of each cohort are z-scored on their own means and SDs: unreliable with few samples
  if (!is.null(batch)) {
    for (i in seq_along(folds)) {
      n_test <- table(batch[setdiff(seq_len(nrow(data)), folds[[i]])])
      small <- n_test[n_test < 10]
      if (length(small) > 0) {
        warning("Fold '", fold_names[i], "': fewer than 10 test samples in cohort(s) ", paste0(names(small), " (", small, ")", collapse = ", "),
                ". With batch, the test samples of each cohort are z-scored on their own means and SDs, so their scores are unreliable",
                if (any(small == 1)) " (with one sample, every score is 0)", ". Use leave-one-cohort-out folds or larger cohorts.", call. = FALSE)
      }
    }
  }

  run_fold <- function(i) {
    train_idx <- folds[[i]]
    test_idx  <- setdiff(seq_len(nrow(data)), train_idx)

    # TRAIN: cell groups learned on the training samples of the fold
    structure <- fit_celltfusion(train_idx, save = FALSE, compute_factors = FALSE)
    train_data <- cbind(data.frame(structure$Cell_groups$Cell_groups), outcome[train_idx, , drop = FALSE])

    # TEST: test samples scored on the cell groups learned on the training samples (same columns as train_data)
    test_data <- data.frame(project_samples(structure, test_idx))
    missing_groups <- setdiff(colnames(structure$Cell_groups$Cell_groups), colnames(test_data))
    if (length(missing_groups) > 0) stop(length(missing_groups), " cell groups could not be scored in the test samples: ", paste(utils::head(missing_groups, 3), collapse = ", "))
    test_data <- test_data[, colnames(structure$Cell_groups$Cell_groups), drop = FALSE]
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

  if (is.null(workers)) workers <- max(1, parallel::detectCores() - 1, na.rm = TRUE)
  if (do_parallel && workers > 1 && length(folds) > 1) {
    # The first fold runs first so the TF and pathway collections cached in Results/ are downloaded only once
    processed_folds <- list(try_fold(1))
    cl <- parallel::makeCluster(min(workers, length(folds) - 1))
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
#' @param p_label Character. Name of the values of \code{p_mat} written in each panel (e.g. \code{"padj"} for adjusted
#'   p-values). Default "p".
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
                                     axis_suffix = "",
                                     p_label = "p") {

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
             " | ", p_label, " = ", format.pval(pv, digits = 2)),
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
#' @param scores A numeric matrix or data frame of TF module scores
#'   (samples x modules), typically the first element (\code{TFs module matrix}) of the
#'   output from \code{compute.WTCNA()}.
#' @param coldata A data frame of sample metadata. Only character, factor and logical
#'   columns are used as traits. Traits with a single level, or with as many levels as
#'   samples, are skipped with a message; samples with a missing value are left out for that trait.
#' @param pval Numeric. ANOVA (adjusted) p-value threshold for significance (\code{<= pval}). Default 0.05.
#' @param file_name Character. Base name appended to output SVG file names
#'   (\code{Results/ANOVA_boxplot_summary_<file_name>_<trait>.svg}, or \code{Results/ANOVA_boxplot_summary_<trait>.svg}
#'   without \code{file_name}). In \code{<trait>}, characters other than letters, digits, ".", "_", "-" and spaces
#'   (e.g. "/") are replaced by "_".
#' @param ncol_grid Integer. Number of columns in the boxplot facet grid; if \code{NULL} (default), the ceiling of the
#'   square root of the number of panels. The size of the figure is
#'   computed from the number of panels, and the y-axis of each panel from its data.
#' @param padj Logical. If TRUE, the ANOVA p-values are corrected across the modules, for each trait, with the
#'   Benjamini-Hochberg method. Default FALSE.
#'
#' @return Invisibly, a data frame with the ANOVA of every module for every trait tested (\code{trait},
#'   \code{module}, \code{F}, \code{p}, \code{p_adjusted}), or \code{NULL} when no trait could be tested.
#'   Called mainly for its side effect (saves one SVG file per categorical trait with at
#'   least one significant module).
#'
#' @keywords internal
compute.metadata.association.boxplot_summary <- function(
    scores, coldata, pval = 0.05, file_name = NULL,
    ncol_grid = NULL, padj = FALSE
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

    df_long <- cbind(as.data.frame(scores), group = trait_values)[keep, , drop = FALSE] %>%
      as.data.frame() %>%
      tidyr::pivot_longer(
        cols = colnames(scores),
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
    anova_df$p_adjusted <- stats::p.adjust(anova_df$p, method = if(padj) "BH" else "none") # Across the modules, for this trait
    anova_all <- rbind(anova_all, data.frame(trait = tr, as.data.frame(anova_df), stringsAsFactors = FALSE))

    sig_pairs <- anova_df[anova_df$p_adjusted <= pval, ]
    if(nrow(sig_pairs) == 0){
      message("No module significantly associated with trait '", tr, "'")
      next
    }
    ncol_trait <- if(is.null(ncol_grid)) ceiling(sqrt(nrow(sig_pairs))) else ncol_grid # Columns of the grid (automatic if NULL)

    feature_labels <- sig_pairs %>%
      dplyr::mutate(
        feature_label = paste0(module, "\nF=", signif(F,3), ", ", if(padj) "p.adj=" else "p=", signif(.data[["p_adjusted"]],3))
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
      ggplot2::facet_wrap(~ module, ncol = ncol_trait, scales = "free_y",
                 labeller = ggplot2::labeller(module = feature_labels)) +
      ggplot2::scale_x_discrete(labels = group_labels) +
      ggplot2::scale_fill_brewer(palette = "Set2") +
      ggplot2::labs(
        x = NULL,
        y = "Module score",
        title = paste0("TF modules associated with ", tr),
        subtitle = paste0("One-way ANOVA, ", if(padj) "BH adjusted p-value" else "p-value", " <= ", pval,
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
      paste0("Results/ANOVA_boxplot_summary", if(!is.null(file_name)) paste0("_", file_name), "_",
             gsub("[^A-Za-z0-9._ -]", "_", tr), ".svg"), # Characters not allowed in file names (e.g. "/") replaced by "_"
      plot = p,
      device = grDevices::svg,
      width = max(5, 3.2 * min(ncol_trait, n_panels) + 0.5),
      height = 3.4 * ceiling(n_panels / ncol_trait) + 1.4,
      limitsize = FALSE
    )
  }

  invisible(anova_all)
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
#' @param coldata A data frame (or tibble) of clinical/survival metadata, with samples as rows.
#' @param time_column Character. Column name in \code{coldata} with the survival/follow-up time.
#' @param event_column Character. Column name in \code{coldata} with the event indicator
#'   (1 = event occurred, 0 = censored).
#' @param file_name Optional character. Suffix used when saving Kaplan-Meier plots to \code{Results/}.
#' @param scores Optional. A samples x features numeric matrix or data frame (e.g.
#'   \code{latent_spaces$Z}). If provided (and \code{group_column} is
#'   \code{NULL}), each feature is tested individually. Rows must be the same samples, in the same order, as the rows of
#'   \code{coldata} (checked by name). Mutually exclusive with \code{group_column}.
#' @param pval Numeric. Log-rank test threshold used to keep a feature as significant when \code{scores} is used:
#'   (adjusted) p-value \code{<= pval}. Default is 0.05.
#' @param quantile_cutoff Numeric between 0 and 1. Quantile cutoff used to split each feature into
#'   High/Low groups when \code{scores} is used (High = at or above the cutoff). Default is 0.5 (median split).
#' @param group_column Optional character. Column name in \code{coldata} defining a
#'   predefined categorical grouping (e.g. cell group membership, cluster, or trait class). Samples without a value
#'   are left out (with a message). Mutually exclusive with \code{scores}.
#' @param padj Logical; if TRUE, the log-rank p-values are corrected for multiple testing across the features with the
#'   Benjamini-Hochberg method, when \code{scores} is used. Default FALSE.
#'
#' @return
#' If \code{group_column} is provided, a list with:
#' \describe{
#'   \item{km_plot}{The \code{survminer::ggsurvplot} object.}
#'   \item{median_PFS}{Median survival time per group (\code{survminer::surv_median()}).}
#'   \item{p_value}{Log-rank test p-value across groups.}
#'   \item{median_follow_up}{Median follow-up time (reverse Kaplan-Meier).}
#' }
#' If \code{scores} is provided, a named list with one \code{Surv() ~ feature} formula per significant feature,
#' and a last element \code{Statistics}: a data frame with the test of every feature (\code{feature}, \code{cutoff},
#' \code{n_high}, \code{n_low}, \code{p_value}, \code{p_adjusted}, \code{significant}). When no feature is significant,
#' the list only contains \code{Statistics} (with a message).
#'
#' A Kaplan-Meier plot (with risk table) is saved as an SVG file to
#' \code{Results/SurvPlot_<group-or-feature>_<file_name>.svg} (\code{Results/SurvPlot_<group-or-feature>.svg} without
#' \code{file_name}; characters not allowed in file names, e.g. "/", replaced by "_"): always for
#' \code{group_column}, and one per significant feature for \code{scores}.
#'
#' @details
#' Requires the \code{survival}, \code{survminer}, and \code{gridExtra} packages (see \code{Suggests}).
#'
#' @examples
#' \dontrun{
#' # Overall survival of the tutorial data (os: months, censOS: 1 = death)
#' data("traitdata.tuto")
#'
#' # Predefined clinical groups
#' compute.survival.analysis(
#'   coldata       = traitdata.tuto,
#'   time_column   = "os",
#'   event_column  = "censOS",
#'   group_column  = "Best.Confirmed.Overall.Response",
#'   file_name     = "Tutorial"
#' )
#'
#' # Automatic screening of latent factor scores (res: output of CellTFusion() on the tutorial data)
#' compute.survival.analysis(
#'   coldata       = traitdata.tuto,
#'   time_column   = "os",
#'   event_column  = "censOS",
#'   scores        = res$Latent_spaces$Z,
#'   pval          = 0.05,
#'   quantile_cutoff = 0.5,
#'   file_name     = "Tutorial"
#' )
#' }
compute.survival.analysis = function(coldata, time_column, event_column, file_name = NULL,
                                     scores = NULL, pval = 0.05, quantile_cutoff = 0.5, group_column = NULL,
                                     padj = FALSE) {

  if (!isTRUE(padj) && !isFALSE(padj)) stop("padj must be TRUE or FALSE")
  if (!requireNamespace("survival", quietly = TRUE) ||
      !requireNamespace("survminer", quietly = TRUE) ||
      !requireNamespace("gridExtra", quietly = TRUE)) {
    stop("Packages 'survival', 'survminer', and 'gridExtra' are required for compute.survival.analysis().")
  }

  coldata <- as.data.frame(coldata) # Tibbles are accepted
  if (!all(c(time_column, event_column) %in% colnames(coldata))) {
    stop("coldata must have the columns '", time_column, "' and '", event_column, "'")
  }
  file_suffix <- if (!is.null(file_name)) paste0("_", file_name) else ""

  ###### CASE 1: Use predefined groups (like patient clusters or risk groups)
  if (!is.null(group_column)) {

    if (!group_column %in% colnames(coldata)) stop("group_column '", group_column, "' is not a column of coldata")
    groups <- coldata[[group_column]]
    keep <- !is.na(groups) # Samples without a group are left out (not a group of their own)
    if (any(!keep)) message(sum(!keep), " samples without a value for '", group_column, "' are left out")

    data_for_model <- data.frame("time" = coldata[[time_column]],
                                 "status" = coldata[[event_column]],
                                 "coxHL" = paste0("Group_", groups))[keep, , drop = FALSE]

    # Fit Kaplan-Meier survival curves across groups
    km_fit <- survival::survfit(survival::Surv(time, status) ~ coxHL, data = data_for_model)

    # Median PFS per group
    medians <- survminer::surv_median(km_fit)

    # Median follow-up (reverse KM)
    rev_km <- survival::survfit(survival::Surv(time, 1 - status) ~ 1, data = data_for_model)
    median_follow_up <- summary(rev_km)$table["median"]

    # Log-rank test for survival difference among groups
    logrank_p <- survminer::surv_pvalue(km_fit, data = data_for_model)$pval

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
      device = grDevices::svg
    )

    return(list(
      km_plot = p,
      median_PFS = medians,
      p_value = logrank_p,
      median_follow_up = median_follow_up
    ))

  ###### CASE 2: Use features to define high/low groups automatically
  } else if (!is.null(scores)) {
    if(!identical(make.names(rownames(scores)), make.names(rownames(coldata)))){
      stop("Samples must be the same and in the same order in scores and coldata rows")
    }
    scores <- as.data.frame(scores, check.names = FALSE)

    # Log-rank test of each feature split into High/Low at its quantile
    tests <- lapply(colnames(scores), function(feature) {
      values <- scores[[feature]]
      cutoff <- stats::quantile(values, quantile_cutoff, na.rm = TRUE)
      data_for_model <- data.frame("time" = coldata[[time_column]],
                                   "status" = coldata[[event_column]],
                                   "coxHL" = factor(ifelse(values >= cutoff, paste0("High_", feature), paste0("Low_", feature))))
      km_fit <- tryCatch(survival::survfit(survival::Surv(time, status) ~ coxHL, data = data_for_model), error = function(e) NULL)
      pval <- if (is.null(km_fit) || length(km_fit$strata) < 2) NA_real_ else survminer::surv_pvalue(km_fit, data = data_for_model)$pval
      list(data = data_for_model, fit = km_fit, cutoff = unname(cutoff), pval = pval,
           n_high = sum(values >= cutoff, na.rm = TRUE), n_low = sum(values < cutoff, na.rm = TRUE))
    })
    p_values <- vapply(tests, function(x) x$pval, numeric(1))
    p_adjusted <- stats::p.adjust(p_values, method = if (padj) "BH" else "none")
    sig <- which(!is.na(p_adjusted) & p_adjusted <= pval)

    statistics <- data.frame(feature = colnames(scores), cutoff = vapply(tests, function(x) x$cutoff, numeric(1)),
                             n_high = vapply(tests, function(x) x$n_high, numeric(1)), n_low = vapply(tests, function(x) x$n_low, numeric(1)),
                             p_value = p_values, p_adjusted = p_adjusted, significant = seq_along(p_values) %in% sig,
                             stringsAsFactors = FALSE)

    significant_combinations <- list()
    for (n in sig) {
      feature <- colnames(scores)[n]
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
        device = grDevices::svg
      )
    }

    threshold <- paste0(if (padj) "adjusted p-value (BH)" else "p-value", " <= ", pval)
    if (length(sig) == 0) {
      message("No feature with a significant log-rank test (", threshold, ")")
    } else {
      message(length(sig), " feature(s) with a significant log-rank test (", threshold, "): ", paste(utils::head(colnames(scores)[sig], 10), collapse = ", "))
    }
    significant_combinations[["Statistics"]] <- statistics
    return(significant_combinations)

  } else {
    stop("Please provide either 'group_column' for predefined groups or 'scores' for feature-based analysis.")
  }
}
