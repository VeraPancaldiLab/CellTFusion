
utils::globalVariables(c("Trait", "Value" ,"level", ".", "Cells_level", "PC1", "Features", "tf", "Module", "shap_value", "mean_shap", "direction", "Feature", "Impact", "values", "ind", "Freq", "new_column", ".data", "id", "value", "P", "p", "sig_p", "r", "i"))

#' Compute one-step CellTFusion
#'
#' @param raw.counts A matrix of raw gene expression counts (genes as rows, samples as columns).
#' @param deconv A data frame with deconvolution features (cell-type proportions as columns x samples as rows).
#' @param normalized Logical; if TRUE, normalize raw counts to log-transformed TPM for TF computation. For deconvolution they are going to be normalize just as TPM. Default is TRUE.
#' @param coldata (Optional) A data frame containing clinical metadata for association analysis with TF modules.
#' @param batch Logical; whether batch correction should be applied where supported. Default is FALSE.
#' @param batch_id Optional character indicating the column name in coldata containing batch identifiers.
#' @param deconv_methods A character vector of deconvolution methods to apply. Default includes:
#'   \code{c("Quantiseq", "Epidish", "DeconRNASeq", "DWLS", "CibersortX")}.
#' @param cbsx.mail (Optional) Email credential for CIBERSORTx. Required if "CibersortX" is among deconv_methods.
#' @param cbsx.token (Optional) Token credential for CIBERSORTx. Required if "CibersortX" is among deconv_methods.
#' @param file_name (Optional) Prefix for output files saved in the "Results/" directory.
#' @param task Analysis mode. Choose between \code{"supervised"} and \code{"unsupervised"}.
#' @param contrast Optional character indicating the condition column used for supervised DEG analysis.
#' @param ref_level Optional character indicating the reference level for supervised DEG analysis.
#' @param TF.collection Character. The source of the TF–target network. Options are `"CollecTRI"` (default), `"Dorothea"`, or `"ARACNE"`.
#' - `"CollecTRI"` and `"Dorothea"` use prebuilt collections from OmnipathR.
#' - `"ARACNE"` allows user input of a custom network file in a 3-column format: `regulator`, `target`, and `mutual information`.
#' @param min_targets_size Integer. Minimum number of target genes per regulon required for TF activity inference. Default is 5.
#' @param universe Optional. A user-specified data frame of TF–target interactions. If not provided, the function will fetch the relevant network based on the `TF.collection` argument.
#' @param paths Optional. A user-specified data frame of pathways gene sets. If not provided, the function will fetch the relevant pathways based on `PROGENy`.
#' @param gene_sets Optional. A data frame of custom gene sets passed to \code{compute.pathway.activity()}'s
#'   \code{gene_sets} argument for GSVA-based scoring. If \code{NULL}, only PROGENy is used.
#' @param minMod Integer; minimum module size for WGCNA module detection.
#' @param corr_mod Numeric; correlation threshold for merging TF modules.
#' @param corr Numeric; correlation threshold used in the deconvolution analysis.
#' @param corr_type Correlation type used in deconvolution analysis. Default is \code{"spearman"}.
#' @param cells_extra A string specifying the cells names to consider and that are not including in the nomenclature of multideconv (see R package)
#' @param pval Numeric; p-value threshold for statistical tests (e.g., metadata and relationship associations).
#' @param enrich_thresh Numeric. Minimum enrichment ratio (foreground/background cell-type frequency)
#'   required to include a cell type in a latent factor's niche. Default is 1.5.
#' @param quantile_cutoff Numeric between 0 and 1. Quantile threshold for selecting top-contributing
#'   cell groups per NMF factor. Default is 0.7.
#' @param return Logical; if TRUE, returns intermediate results from internal functions. Default is TRUE.
#' @param cancer_type Character. TCGA cancer type abbreviation (e.g., \code{"blca"}, \code{"brca"},
#'   \code{"cesc"}, \code{"chol"}, \code{"coad"}, \code{"skcm"}) used to load the corresponding
#'   TCGA meta-programs for TME state mapping. Must match an available meta-program file.
#' @param verbose Boolen value to whether print or no the function messages
#'
#' @return A list containing:
#' \describe{
#'   \item{Deconvolution}{A matrix with cell-type proportions (samples as rows, cell types as columns).}
#'   \item{TFs_matrix}{A matrix with TF activity scores (samples as rows, TFs as columns).}
#'   \item{TF_network}{A list representing the TF module network and related WGCNA output.}
#'   \item{Pathways_scores}{A matrix of pathway activity scores.}
#'   \item{Processed_deconvolution}{An object with the processed deconvolution analysis results.}
#'   \item{Cell_groups}{A matrix of scores representing the cell groups across samples.}
#' }
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#' data("raw.counts.tuto")
#' data("traitdata.tuto")
#'
#' res <- CellTFusion(
#'   raw.counts = raw.counts.tuto,
#'   normalized = TRUE,
#'   coldata = traitdata.tuto,
#'   deconv_methods = c("Quantiseq", "DeconRNASeq"),
#'   file_name = "TestRun",
#'   min_targets_size = 15,
#'   minMod = 20,
#'   corr_mod = 0.25,
#'   corr = 0.7,
#'   pval = 0.05
#' )
#'}
#'
CellTFusion = function(raw.counts, deconv = NULL, normalized = T, coldata = NULL, batch = F, batch_id = NULL, deconv_methods = c("Quantiseq", "CBSX", "Epidish", "DeconRNASeq", "DWLS"), cbsx.mail = NULL, cbsx.token = NULL, file_name = NULL, task = c("supervised", "unsupervised"),
                       contrast = NULL, ref_level = NULL, TF.collection = "CollecTRI", min_targets_size = 3, universe = NULL, paths = NULL, gene_sets = NULL, minMod = 3, corr_mod = 0.9, corr = 0.7, corr_type = "spearman", cells_extra = NULL, pval = 0.05, enrich_thresh = 1.5, quantile_cutoff = 0.7, cancer_type, return = T, verbose = T){

  set.seed(123)

  #Normalize counts
  if(normalized == T){
    counts.norm = data.frame(ADImpute::NormalizeTPM(raw.counts, log = T))
  }else{
    counts.norm = raw.counts
  }

  #Extract batch column if TRUE
  if(batch){
    batch_vec = coldata[,batch_id]
  }else{
    batch_vec = NULL
  }

  #Deconvolution
  if(is.null(deconv)){
    if(verbose){
      cat("Calculating cell type deconvolution............................................................\n")
    }
    if(("CBSX" %in% deconv_methods) == T){
      if(is.null(cbsx.mail)==T || is.null(cbsx.token)==T){
        stop("No CBSX credentials given!\n")
      }else{
        deconv = multideconv::compute.deconvolution(raw.counts, normalized = normalized, methods = deconv_methods, credentials.mail = cbsx.mail, credentials.token = cbsx.token, doParallel = T, workers = 3, file_name = file_name, return = return)
      }
    }else{
      deconv = multideconv::compute.deconvolution(raw.counts, normalized = normalized, methods = deconv_methods, file_name = file_name, return = return)
    }
  }

  #TF activity
  if (verbose) {
    cat("\nCalculating TF activity............................................................\n")
  }

  if(task == "supervised"){ ### missing to add batch variable option
    if (verbose) {
      cat("\nRunning supervised task............................................................\n")
    }
    if(!is.null(contrast) && !is.null(ref_level)){
      if(normalized){
        res_deg = run_deg_analysis(raw.counts, coldata, contrast, ref_level = ref_level) %>%
          dplyr::select(t) # compute DEGs for the contrast of interest
      }else{
        stop("For differential expression analysis, raw counts must be provided")
      }

      tfs_deg = compute.TFs.activity(res_deg, TF.collection, min_targets_size, universe, return = return, file.name = file_name) # compute TFs activity using DEGs as input
      tfs_mat <- compute.TFs.activity(counts.norm, TF.collection, min_targets_size, universe, return = return, file.name = file_name) # compute TFs activity using all genes as input
      tfs = tfs_mat[,colnames(tfs_mat) %in% colnames(tfs_deg)] # Subset the TFs matrix to keep only those TFs that are significant in the DEGs analysis
    }else{
      stop("For supervised analysis, contrast and ref_level must be provided")
    }
  }else{
    if (verbose) {
      cat("\nRunning unsupervised task............................................................\n")
    }
    if (!batch) {  ## No batch → single matrix
      tfs <- compute.TFs.activity(counts.norm, TF.collection, min_targets_size, universe, return = return, file.name = file_name)
    }else {
      if(is.null(coldata) || is.null(batch_id)) {
        stop("When batch = TRUE, coldata and batch_id must be provided")
      }

      batch_vec <- coldata[, batch_id]

      if (length(batch_vec) != ncol(counts.norm)) {
        stop("batch_id column must match number of samples")
      }

      batch_vec <- factor(batch_vec, levels = unique(batch_vec))
      cohorts <- split(seq_len(ncol(counts.norm)), batch_vec)

      tfs <- lapply(names(cohorts), function(cohort) {
        idx <- cohorts[[cohort]]
        compute.TFs.activity(counts.norm[, idx, drop = FALSE], TF.collection, min_targets_size, universe, file.name = file_name)
      })

      names(tfs) <- names(cohorts)
    }
  }

  # 1. TFs network construction
  if(verbose){
    cat("\nConstructing TF network............................................................\n")
  }
  network = compute.WTCNA(TFs.matrix = tfs, batch = batch, network.type = "signed", clustering.method = "ward.D2", minMod, corr_mod, cor_type = "p", return = return, file.name = file_name)

   # 2. Deconvolution analysis and cell groups construction
  # 1.2. Modules characterization
  # cat("\nPerforming TF module characterization............................................................\n")
  # hub_tfs = identify_hub_TFs(t(tfs), network, MM_thresh = 0.8, degree_thresh = 0.9)
  # compute.modules.enrichment(counts.norm, hub_tfs)

  # 2. Pathways activity inference: only needed for dictionary
  if(verbose){
    cat("\nCalculating pathway activities............................................................\n")
  }
  pathways = compute.pathway.activity(counts.norm, gene_sets = gene_sets, paths = paths, return = return, file.name = file_name)
  
  # 3. Deconvolution analysis
  if(verbose){
    cat("\nPerforming deconvolution analysis............................................................\n")
  }

  dt = compute.deconvolution.analysis(deconv, corr = corr, corr_type = corr_type, seed = 123, batch = batch_vec, cells_extra = cells_extra, file_name = file_name, return = return, verbose = FALSE)
  dt = deconvolution_dictionary(dt, pathways, batch_id = batch_vec) ## Apply dictionary of deconvolution (To be added in compute.deconvolution.analysis() soon)

  # 4. Cell groups construction and scores
  if(verbose){
    cat("\nCell groups identification............................................................\n")
  }
  cell.groups = construct_cell_groups(network, dt, batch = batch_vec, pval = pval, clustering.method = "ward.D2", dendrogram_file = file_name, return_dendrogram = return)

  # 5. Cell groups latent spaces
  if(verbose){
    cat("\nLatent spaces calculation............................................................\n")
  }
  latent_spaces <- compute.latent_factors(cell.groups[[1]], rank = NULL, seed = 123, file_name = file_name, return = return)

  # 6. Extract cell niches
  if(verbose){
    cat("\nExtract cell niches from latent spaces............................................................\n")
  }
  cells_niches = compute_cells_niches(latent_spaces, dt, cell.groups, enrich_thresh = enrich_thresh, quantile_cutoff = quantile_cutoff, cells_extra = cells_extra, return = return, file_name = file_name)

  # 7. Find TME states
  if(verbose){
    cat("\nCompute TME states............................................................\n")
  }
  gsea_results <- compute_factor_gsea(
    RNA.tpm      = counts.norm,          # genes x samples
    features_df  = latent_spaces$Z,      # samples x factors
    plot_dot     = return,
    top_n        = 10,
    file_name    = file_name
  )

  cat("\nMapping to metaprograms............................................................\n")
  metaprograms_mapping <- map_to_metaprograms(
    gsea_study   = gsea_results,
    cancer_type  = cancer_type,
    nes_thresh   = 1.0
  )

  if(verbose){
    cat("\nEverything done! Results are saved in Results/ folder............................................................\n")
  }

  res = list("Deconvolution" = deconv, "TFs_matrix" = tfs, "TF_network" = network, "Pathways_scores" = pathways,
             "Processed_deconvolution" = dt, "Cell_groups" = cell.groups, "Latent_spaces" = latent_spaces, 
             "Cells_niches" = cells_niches, "TME_states" = metaprograms_mapping)

  return(res)

}

#' Compute cell group scores from deconvolution and TF module network
#'
#' This function identifies cell groups based on dendrogram cuts, computes composite scores for each group
#' using deconvolution features and TF module networks, and optionally exports the results.
#'
#' @param deconvolution A data frame with deconvolution features (typically a cell-type or cluster x sample matrix).
#' This is usually the first element returned by \code{multideconv::compute.deconvolution.analysis()}.
#'
#' @param cell.dendrograms A named list of dendrogram objects, each corresponding to a TF module, typically
#' returned by \code{identify.cell.groups()}.
#'
#' @param tfs.module.network A list containing network information of transcription factor (TF) modules, as
#' obtained from \code{compute.WTCNA()}. It should contain at least one element with TF module membership or connectivity.
#' @param batch Optional vector indicating batch assignment for samples.
#'
#' @param return Logical; if TRUE (default), writes CSV files with cell group compositions and scores
#' to the "Results/" folder.
#'
#' @return A list of three elements:
#' \describe{
#'   \item{scores}{A data frame with the composite scores of all identified cell groups across samples.}
#'   \item{composition}{A list of vectors indicating the composition (original features) of each cell group.}
#'   \item{loadings}{A list of loadings (feature contributions) for each cell group.}
#' }
#' If \code{return=TRUE}, two CSV files will be created:
#' \itemize{
#'   \item{\code{Results/Cell.groups.composition.csv}: A table showing the composition of each cell group.}
#'   \item{\code{Results/Cell.groups.scores.csv}: A matrix of cell group scores across samples.}
#' }
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
cell.groups.computation = function(deconvolution, cell.dendrograms, tfs.module.network, batch = NULL, return = T, pval = 0.05, n_perm = 999, dendrogram_file = NULL, return_dendrogram = FALSE){

  cell.groups = list()
  cuts_per_dendrogram = calculate_dendrogram_cuts(cell.dendrograms)

  if(return_dendrogram){
    plot_dendrogram_clusters(cell.dendrograms, cuts_per_dendrogram, file_name = dendrogram_file)
  }

  for(k in 1:length(cell.dendrograms)){
    clusters = cuts_per_dendrogram[[k]]
    names(clusters) = labels(as.dendrogram(cell.dendrograms[[k]]))

    cluster_ids = sort(unique(clusters[clusters > 0]))  # exclude unassigned (label 0)

    y = list()
    x = list()
    z = list()

    for (j in seq_along(cluster_ids)) {
      cells = names(clusters)[clusters == cluster_ids[j]]
      y[[j]] = cells
      ###################################################Compute score for each cell group
      pca_group = deconvolution[, colnames(deconvolution) %in% cells, drop = F]
      color = names(cell.dendrograms)[k]
      score = compute_composite_score(pca_group, color, tfs.module.network, discard = T, batch = batch, pval = pval, n_perm = n_perm)
      x[[j]] = score[[1]]
      z[[j]] = score[[2]]
    }

    ###################################################Remove groups with no corr components
    idy = which(x == "NA")
    if(length(idy) > 0){
      y = y[-idy]
      x = x[-idy]
      z = z[-idy]
    }

    ###################################################Remove groups with one single feature (deprecated)
    if(length(x) != 0){
      groups_cut = remove_single_groups(x, y, z)
      cell.groups[[k]] = if(!is.null(groups_cut)) list(groups_cut[[1]], groups_cut[[2]], groups_cut[[3]]) else NA
    }else{
      cell.groups[[k]] <- NA #To avoid the problem when putting NULL (list ignore it and not consider it an element)
    }
  }

  ###################################################Naming of cell groups
  names(cell.groups) = names(cell.dendrograms) #Named before removing NA to avoid omitting dendrogram names
  cell.groups <- cell.groups[!sapply(cell.groups, function(x) all(is.na(x)))]

  if(length(cell.groups)==0){
    stop("No cell groups founds with significant features")
  }

  for (i in 1:length(cell.groups)) {
    for (j in 1:length(cell.groups[[i]][[1]])) {
      names(cell.groups[[i]][[1]])[j] = paste0("Dendrogram_", names(cell.groups)[i], "_group_", j)
      names(cell.groups[[i]][[2]])[j] = paste0("Dendrogram_", names(cell.groups)[i], "_group_", j)
    }
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
    utils::write.csv(clusters, "Results/Cell.groups.composition.csv", row.names = F)
    utils::write.csv(cell.groups[[1]], "Results/Cell.groups.scores.csv")
  }

  return(cell.groups)
}

#' Compute associations between TF module scores and clinical metadata
#'
#' This function tests for associations between transcription factor (TF) module scores
#' and available clinical traits. It uses Pearson correlation for continuous (numeric) traits,
#' and ANOVA for categorical traits. Results are visualized as a labeled heatmap and violin plots.
#' All plots are saved in the `Results/` directory.
#'
#' @param tfs.modules A numeric matrix or data frame of TF module scores across samples.
#'        Typically the output from `compute.WTCNA()`. Rows represent samples, columns represent TF modules.
#' @param coldata A data frame containing clinical traits (both categorical and numerical) for the same samples.
#'        Row names should match those of `tfs.modules`.
#' @param pval A numeric threshold (default = 0.05) to determine statistical significance.
#'        Only associations with p-values below this threshold are considered significant in the heatmap.
#' @param corr_method Character string specifying the correlation method to use for continuous variables.
#'        Options are `"p"` for pearson and `"s"` for spearman. Default is `"p"`.
#' @param file.name Character. Base file name for saving PDF plots of results.
#' @param width A numeric value indicating the width (in inches) of the output heatmap plot (default = 20).
#' @param height A numeric value indicating the height (in inches) of the output heatmap plot (default = 8).
#' @param ncol Integer. Number of columns in the grid of association boxplots.
#' @param y_min Numeric. Lower y-axis limit for grid boxplots.
#' @param y_max Numeric. Upper y-axis limit for grid boxplots.
#' @param plot_grid Logical; if TRUE, generates a grid of boxplot summaries.
#' @param width_grid Numeric width of the grid plot output.
#' @param height_grid Numeric height of the grid plot output.
#'
#' @return This function saves the following to the `Results/` directory:
#' \itemize{
#'   \item A labeled heatmap showing Pearson correlations and ANOVA test p-values.
#'   \item Individual violin plots for significant categorical trait associations.
#' }
#' The function does not return an object to the R environment.
#'
#' @export
#'
#' @examples
#'
#' data("network.tuto")
#' data("traitdata.tuto")
#'
#' compute.metadata.association(
#'   tfs.modules = network.tuto[[1]],
#'   coldata = traitdata.tuto,
#'   pval = 0.05,
#'   file.name = 'Tutorial',
#'   width = 15,
#'   height = 10
#' )
#'
compute.metadata.association <- function(
    tfs.modules,
    coldata,
    pval = 0.05,
    corr_method = "p",
    file.name,
    width = 20,
    height = 8,
    ncol = 5,
    y_min = 0,
    y_max = 0.5,
    plot_grid = F,
    width_grid = 18,
    height_grid = 10
) {
  ### Association with categorical variables
  coldata_categorical <- coldata %>%
    dplyr::select(dplyr::where(is.character) | dplyr::where(is.factor))

  if(ncol(coldata_categorical) != 0) {
    if(plot_grid){
      # Use the new wrapper function to generate boxplot summaries
      compute.metadata.association.boxplot_summary(
        tfs.modules = tfs.modules,
        coldata = coldata_categorical,
        pval = pval,
        file.name = file.name,
        ncol = ncol,
        y_min = y_min,
        y_max = y_max,
        width = width_grid,
        height = height_grid
      )
    }
  }

  ### Association with quantitative variables
  coldata_quantitative <- coldata %>%
    dplyr::select(dplyr::where(is.numeric))

  if(ncol(coldata_quantitative) != 0){
    moduleTraitCor <- WGCNA::cor(
      tfs.modules,
      coldata_quantitative,
      method = corr_method,
      use = "pairwise.complete.obs"
    )
    moduleTraitPvalue <- WGCNA::corPvalueStudent(moduleTraitCor, nrow(tfs.modules))

    #### Replace p-values for significance labels
    breaks <- c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf)
    labels <- c("****", "***", "**", "*", "")

    moduleTraitPvalue <- matrix(
      cut(moduleTraitPvalue, breaks = breaks, labels = labels, right = FALSE),
      nrow = nrow(moduleTraitPvalue),
      dimnames = dimnames(moduleTraitPvalue)
    )

    textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", moduleTraitPvalue, ")", sep = "")
    dim(textMatrix) <- dim(moduleTraitCor)

    # Set non-significant entries to NA
    idx <- which(moduleTraitPvalue == "" | moduleTraitPvalue > pval)
    for (i in idx) textMatrix[i] <- NA

    pdf(paste0("Results/TF.modules_metadata_", file.name), width = width, height = height)
    par(mar = c(25, 15, 3, 3))
    WGCNA::labeledHeatmap(
      Matrix = moduleTraitCor,
      xLabels = colnames(moduleTraitCor),
      yLabels = rownames(moduleTraitCor),
      ySymbols = rownames(moduleTraitCor),
      colorLabels = FALSE,
      colors = WGCNA::blueWhiteRed(50),
      textMatrix = textMatrix,
      setStdMargins = FALSE,
      cex.text = 0.7,
      zlim = c(-1,1),
      main = paste0(
        "Clinical associations ", file.name,
        "\nOnly showing significant associations (pvalue < ", pval, ")"
      )
    )
    dev.off()
  }
}

#' Compute TF module enrichment using directed target genes
#'
#' This function performs enrichment analysis for transcription factor (TF) modules
#' using known TF-target interactions. For each module, it identifies hub TFs and retrieves
#' their known targets from the CollecTRI database. It then conducts an over-representation
#' analysis (ORA) against the Reactome pathway database. To reduce redundancy, only unique
#' pathways per module are retained by filtering out overlaps between modules.
#'
#' @param RNA.tpm A numeric matrix of normalized gene expression values.
#'        Rows are gene symbols, columns are samples. This matrix is used to restrict the universe
#'        of genes used in enrichment and validate TF-target gene presence.
#' @param hub_tfs A list of hub TFs per module, typically obtained using `identify_hub_TFs()`.
#'        This should be a list of two elements: one with named hub TFs per module, and another
#'        with their module membership or additional metadata.
#'
#' @return No object is returned. For each module with significant enrichment (p-value < 0.05),
#'         a dot plot is saved in the `Results/` directory as a PDF file named
#'         `Module <color>.pdf`. If no enrichment is found for a module, a message is printed
#'         and no file is saved for that module.
#'
#' @details
#' The function uses the `decoupleR::get_collectri()` function to obtain TF-target relationships,
#' and `clusterProfiler::enrichPathway()` for ORA using the Reactome database.
#' Pathways shared between multiple modules are filtered using a Venn diagram-based comparison
#' to retain only module-specific results.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Identify hub TFs from a TF activity matrix and WGCNA network
#' hub_tfs <- identify_hub_TFs(t(tfs_activity), wgcna_network, MM_thresh = 0.8, degree_thresh = 0.9)
#'
#' # Perform module-level Reactome enrichment
#' compute.modules.enrichment(counts.norm, hub_tfs)
#' }
compute.modules.enrichment <- function(RNA.tpm, hub_tfs){
  tf_cache_file <- "Results/TF_target_collection.csv"
  if (file.exists(tf_cache_file)) {
    net = utils::read.csv(tf_cache_file, row.names = 1)
  } else {
    net = decoupleR::get_collectri(organism = 'human', split_complexes = F)
    utils::write.csv(net, tf_cache_file)
  }
  res = list()
  pathways = list()

  #Pathway enrichment using target genes
  for (i in 1:length(hub_tfs[[1]])) {
    color = names(hub_tfs[[1]])[i]
    res[[i]] = module_enrich(as.matrix(RNA.tpm), color, hub_tfs, net)
    names(res)[i] = color
    pathways[[i]] = res[[i]]@result[["Description"]]
    names(pathways)[i] = color
  }

  ItemsList <- gplots::venn(pathways, show.plot = FALSE)
  x = attributes(ItemsList)

  #Keep only unique pathways
  for (i in 1:length(hub_tfs[[1]])) {
    color = names(hub_tfs[[1]])[i]
    if(is.null(res[[i]]) == T){
      print(paste0("No enrichment for module ", color))
    }else{
      res[[i]]@result = res[[i]]@result[which(res[[i]]@result[["Description"]]%in%x$intersections[[color]]),] #take unique pathways per module
      if(nrow(res[[i]]@result)==0){
        print(paste0("No enrichment for module ", color))
      }else{
        pdf(paste0("Results/Enrichment by Reactome\nModule ", color))
        print(enrichplot::dotplot(res[[i]],  title=paste0("Enrichment by Reactome\nModule ", color)))
        dev.off()
      }
    }
  }

}

#' Compute modules relationship
#'
#' Performs linear correlation between two matrices (e.g., TF modules and external features such as deconvolution estimates or pathway activities)
#' across the same set of samples, visualizing only significant associations.
#'
#' @param matA A numeric matrix or data frame of features (samples x features).
#' @param matB A numeric matrix or data frame of features to correlate with (samples x features).
#' @param file_name A string indicating the base name (without extension) of the figure to be saved in the "Results/" folder.
#' @param batch Optional vector indicating batch assignment for samples. If provided, partial correlations are computed controlling for batch.
#' @param width An integer indicating the width (in inches) of the output PDF figure. Default is 8.
#' @param height An integer indicating the height (in inches) of the output PDF figure. Default is 8.
#' @param par_mar A numeric vector of length 4 specifying the margin sizes (bottom, left, top, right) for the heatmap. If NULL (default), reasonable defaults are chosen based on plot orientation.
#' @param pval A numeric threshold for the p-value to define significance. Default is 0.05.
#' @param padj Logical; if TRUE, applies Bonferroni correction for multiple testing. Default is FALSE.
#' @param cor_type Type of correlation to compute: "p" (Pearson), "s" (Spearman), or "k" (Kendall). Default is "p".
#' @param return Logical; if TRUE, the function returns a list containing the correlation matrix and a named list of significant feature names per module. Default is FALSE.
#' @param vertical Logical; if TRUE, produces a vertical heatmap (traits on x-axis, modules on y-axis). Otherwise, a horizontal layout is used. Default is FALSE.
#' @param plot Logical; if TRUE, saves the heatmap plot as a PDF. Default is TRUE.
#' @param plot.grid Logical; if TRUE, generates per-pair scatter grid plots for significant associations.
#' @param width.grid Numeric width of the scatter grid output.
#' @param height.grid Numeric height of the scatter grid output.
#' @param ncol.grid Integer number of columns used in scatter grid layout.
#'
#' @return If `return = TRUE`, returns a list with:
#' \itemize{
#'   \item A correlation matrix between modules and external features.
#'   \item A named list with significant features per module (after p-value or adjusted p-value thresholding).
#' }
#' If `return = FALSE`, the function saves a heatmap of significant correlations to "Results/{file_name}.pdf".
#'
#' 
#' @details
#' The function assumes that `matA` and `matB` share the same rownames (i.e., samples in the same order). The correlation is computed using WGCNA's
#' `cor()` and `corPvalueStudent()` functions. Insignificant correlations (based on p-value or adjusted p-value) are excluded from the visualization.
#'
#' @examples
#' data("counts.norm.tuto")
#' data("network.tuto")
#'
#' pathways <- compute.pathway.activity(counts.norm.tuto)
#' compute.modules.relationship(network.tuto[[1]],
#'                              pathways,
#'                              "Pathways_Progeny-TFs_Modules",
#'                              width = 15)
#'
#' data("deconv_subgroups.tuto")
#' corr = compute.modules.relationship(network.tuto[[1]],
#'                                     deconv_subgroups.tuto[[1]],
#'                                     "Deconvolution-TFs_Modules",
#'                                     plot = FALSE,
#'                                     return = TRUE,
#'                                     pval = 0.01)
#'
#' @export
compute.modules.relationship <- function(matA, matB, file_name, batch = NULL, width = 8, height = 8, par_mar = NULL, pval=0.05, padj = F, cor_type = "p", return = F, vertical = F, plot = T, plot.grid = F, width.grid = 12, height.grid = 10, ncol.grid = NULL){

  matA = data.frame(matA)
  matB = data.frame(matB)

  ##check if names from both features are the same
  if(all(rownames(matA)==rownames(matB)) == F){
    stop("No equal names, verify the input objects")
  }

  # ---------- PARTIAL CORRELATION IF BATCH PROVIDED ----------
  if(!is.null(batch)){
    if(length(batch) != nrow(matA)) stop("Length of batch must match number of samples")

    # Convert batch to numeric if it is factor or character
    if(is.factor(batch) || is.character(batch)){
      batch <- as.numeric(as.factor(batch))
    }

    moduleTraitCor <- matrix(NA, nrow = ncol(matA), ncol = ncol(matB))
    moduleTraitPvalue <- matrix(NA, nrow = ncol(matA), ncol = ncol(matB))
    rownames(moduleTraitCor) <- colnames(matA)
    colnames(moduleTraitCor) <- colnames(matB)
    rownames(moduleTraitPvalue) <- colnames(matA)
    colnames(moduleTraitPvalue) <- colnames(matB)

    for(i in 1:ncol(matA)){
      for(j in 1:ncol(matB)){
        pc <- ppcor::pcor.test(matA[,i], matB[,j], batch)
        moduleTraitCor[i,j] <- pc$estimate
        moduleTraitPvalue[i,j] <- pc$p.value
      }
    }
  } else {
    # ---------- STANDARD CORRELATION ----------
    moduleTraitCor = WGCNA::cor(matA, matB, method = cor_type)
    moduleTraitPvalue = WGCNA::corPvalueStudent(moduleTraitCor, nrow(matA))
  }

  if (plot.grid) {
    plot.module.scatter.grid(
      matA = matA,
      matB = matB,
      cor_mat = moduleTraitCor,
      p_mat = moduleTraitPvalue,
      file_name = file_name,
      pval = pval,
      cor_type = cor_type,
      width = width.grid,
      height = height.grid,
      ncol = ncol.grid,
      only_sig = TRUE
    )
  }

  # rev = which(colSums(moduleTraitPvalue > pval)==nrow(moduleTraitPvalue)) #check if there are features no significant with any module
  #
  # if(length(rev)>0){
  #   moduleTraitCor = moduleTraitCor[,-rev]
  #   moduleTraitPvalue = moduleTraitPvalue[,-rev]
  # }

  if(padj == T){
    for (i in 1:ncol(moduleTraitPvalue)) {
      moduleTraitPvalue[,i] = stats::p.adjust(moduleTraitPvalue[,i], method = 'bonferroni')
    }
  }

  ##Plot in vertical
  if(vertical == T){
    if(ncol(data.frame(moduleTraitCor))>1){
      #Extract significant features per trait
      sig = list()
      for (i in 1:nrow(moduleTraitPvalue)) {
        sig[[i]] = names(which(signif(moduleTraitPvalue[i,],2)<=pval))
      }
      names(sig) = rownames(moduleTraitPvalue)

      if(return == T){
        retu = list(moduleTraitCor, sig)
        return(retu)
      }else{
        d <- stats::dist(t(moduleTraitCor), method = "euclidean")
        hc1 <- stats::hclust(d, method = "ward.D2")
        vec = hc1[["order"]]
        textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 2), ")", sep = "")
        dim(textMatrix) = dim(moduleTraitCor)
        idx = which(round(moduleTraitPvalue,2)>pval)
        for (i in idx) {
          textMatrix[i] = NA
        }
        textMatrix = t(textMatrix)
        moduleTraitCor = data.frame(t(moduleTraitCor))
        if(plot){
          if(is.null(par_mar)){
            par_mar = c(3, 25, 5, 3)
          }

          pdf(paste0("Results/",file_name), width = width, height = height)
          par(mar = par_mar)
          WGCNA::labeledHeatmap(Matrix = moduleTraitCor[vec,],
                                xLabels = colnames(moduleTraitCor),
                                yLabels = rownames(moduleTraitCor[vec,]),
                                xLabelsPosition = "top",
                                colors = WGCNA::blueWhiteRed(50),
                                textMatrix = textMatrix[vec,],
                                setStdMargins = F,
                                cex.text = 0.5,
                                zlim = c(-1,1))
          dev.off()
        }
      }}else{
        #Extract significant features per trait
        sig = list()
        for (i in 1:nrow(moduleTraitPvalue)) {
          sig[[i]] = colnames(moduleTraitPvalue)[which(signif(moduleTraitPvalue[i,],2)<=pval)]
        }
        names(sig) = rownames(moduleTraitPvalue)
        if(return == T){
          retu = list(moduleTraitCor, sig)
          return(retu)
        }else{
          textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 2), ")", sep = "")
          idx = which(round(moduleTraitPvalue,2)>pval)
          for (i in idx) {
            textMatrix[i] = NA
          }
          textMatrix = t(textMatrix)
          moduleTraitCor = data.frame(t(moduleTraitCor))
          colnames(moduleTraitCor)[1] = colnames(matB)
          if(plot){
            if(is.null(par_mar)){
              par_mar = c(25, 15, 3, 3)
            }

            pdf(paste0("Results/",file_name), width = width, height = height)
            par(mar = par_mar)
            WGCNA::labeledHeatmap(Matrix = moduleTraitCor,
                                  xLabels = colnames(moduleTraitCor),
                                  yLabels = rownames(moduleTraitCor),
                                  xLabelsPosition = "top",
                                  colors = WGCNA::blueWhiteRed(50),
                                  textMatrix = textMatrix,
                                  setStdMargins = F,
                                  cex.text = 0.5,
                                  zlim = c(-1,1))
            dev.off()
          }
        }}}


  ###Plot in horizontal
  if(vertical == F){
    if(ncol(data.frame(moduleTraitCor))>1){
      #Extract significant features per trait
      sig = list()
      for (i in 1:nrow(moduleTraitPvalue)) {
        sig[[i]] = names(which(signif(moduleTraitPvalue[i,],2)<=pval))
      }
      names(sig) = rownames(moduleTraitPvalue)
      if(return==T){
        retu = list(moduleTraitCor, sig)
        return(retu)
      }else{
        d <- stats::dist(t(moduleTraitCor), method = "euclidean")
        hc1 <- stats::hclust(d, method = "ward.D2")
        vec = hc1[["order"]]
        textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 2), ")", sep = "")
        dim(textMatrix) = dim(moduleTraitCor)
        idx = which(round(moduleTraitPvalue,2)>pval)
        for (i in idx) {
          textMatrix[i] = NA
        }
        moduleTraitCor = data.frame(moduleTraitCor)
        if(plot){
          if(is.null(par_mar)){
            par_mar = c(25, 15, 3, 3)
          }

          pdf(paste0("Results/",file_name), width = width, height = height)
          par(mar = par_mar)
          WGCNA::labeledHeatmap(Matrix = moduleTraitCor[,vec],
                         xLabels = names(moduleTraitCor[,vec]),
                         yLabels = rownames(moduleTraitCor),
                         ySymbols = rownames(moduleTraitCor),
                         colorLabels = FALSE,
                         colors = WGCNA::blueWhiteRed(50),
                         textMatrix = textMatrix[,vec],
                         setStdMargins = FALSE,
                         cex.text = 0.5,
                         zlim = c(-1,1),
                         main = paste("Module-trait relationships"))
          dev.off()
        }
      }}else{
        #Extract significant features per trait
        sig = list()
        for (i in 1:nrow(moduleTraitPvalue)) {
          sig[[i]] = colnames(moduleTraitPvalue)[which(signif(moduleTraitPvalue[i,],2)<=pval)]
        }
        names(sig) = rownames(moduleTraitPvalue)

        if(return == T){
          retu = list(moduleTraitCor, sig)
          return(retu)
        }else{
          textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 2), ")", sep = "")
          idx = which(round(moduleTraitPvalue,2)>pval)
          for (i in idx) {
            textMatrix[i] = NA
          }
          moduleTraitCor = data.frame(moduleTraitCor)
          colnames(moduleTraitCor)[1] = colnames(matB)
          if(plot){
            if(is.null(par_mar)){
              par_mar = c(25, 15, 3, 3)
            }

            pdf(paste0("Results/",file_name), width = width, height = height)
            par(mar = par_mar)
            WGCNA::labeledHeatmap(Matrix = moduleTraitCor,
                           xLabels = names(moduleTraitCor),
                           yLabels = rownames(moduleTraitCor),
                           ySymbols = rownames(moduleTraitCor),
                           colorLabels = FALSE,
                           colors = WGCNA::blueWhiteRed(50),
                           textMatrix = textMatrix,
                           setStdMargins = FALSE,
                           cex.text = 0.5,
                           zlim = c(-1,1),
                           main = paste("Module-trait relationships"))
            dev.off()
          }
        }}}
}

#' Computes TF-modules pathway activities scores
#'
#' This function computes pathway activity scores from normalized gene expression data using a multivariate linear model (MLM) based on the PROGENy resource (Schubert et al., 2018).
#' Optionally, it also performs Gene Set Variation Analysis (GSVA) using hallmark signatures or any user-provided gene sets.
#'
#' @param RNA.tpm A numeric matrix of normalized gene expression values with genes as rows and samples as columns.
#' @param gene_sets A list of gene sets (e.g., hallmark signatures or user-defined sets). If provided, GSVA scores will be computed for these sets. Default is \code{NULL}.
#' @param paths A data frame describing the pathway-gene interactions for use with PROGENy. If \code{NULL}, the human PROGENy resource (top 500 genes) will be used by default.
#' @param return Logical; if TRUE, saves matrices in Results/ folder. Default is TRUE.
#' @param file.name Optional character suffix used when writing output CSV files.
#'
#' @return If \code{gene_sets} is \code{NULL}, a scaled matrix of PROGENy pathway activity scores (samples as rows, pathways as columns).
#' If \code{gene_sets} is provided, a list with two elements:
#' \itemize{
#'   \item \code{sample_acts_progeny}: A scaled matrix of PROGENy pathway activity scores.
#'   \item \code{sample_acts_gsva}: A scaled matrix of GSVA scores based on the provided gene sets.
#' }
#'
#' @export
#'
#' @references
#' Schubert M, Klinger B, Klünemann M, Sieber A, Uhlitz F, Sauer S, Garnett MJ, Blüthgen N, Saez-Rodriguez J.
#' Perturbation-response genes reveal signaling footprints in cancer gene expression. Nature Communications. 2018. \doi{10.1038/s41467-017-02391-6}
#'
#' @examples
#' # Compute only PROGENy activities
#' data("counts.norm.tuto")
#' pathways <- compute.pathway.activity(counts.norm.tuto)
#'
compute.pathway.activity <- function(RNA.tpm, gene_sets = NULL, paths = NULL, return = TRUE, file.name = NULL) {

  RNA.tpm <- as.matrix(RNA.tpm)
  results_list <- list()

  ###### PROGENy
  progeny_cache_file <- "Results/Pathways_collection_PROGENy.csv"
  if (is.null(paths)) {
    if (file.exists(progeny_cache_file)) {
      paths <- utils::read.csv(progeny_cache_file, row.names = 1)
    } else {
      paths <- decoupleR::get_progeny(organism = "human", top = 500)
      utils::write.csv(paths, progeny_cache_file)
    }
  }

  progeny <- decoupleR::run_mlm(
    mat      = RNA.tpm,
    net      = paths,
    .source  = "source",
    .target  = "target",
    .mor     = "weight",
    minsize  = 5
  )

  sample_acts_progeny <- progeny %>%
    tidyr::pivot_wider(id_cols = "condition", names_from = "source", values_from = "score") %>%
    tibble::column_to_rownames("condition") %>%
    as.matrix() %>%
    scale() %>%
    as.data.frame()

  results_list$PROGENy <- sample_acts_progeny

  ###### GSVA (optional)
  if (!is.null(gene_sets)) {
    gsva_results <- GSVA::gsva(
      RNA.tpm,
      gene_sets,
      method  = "gsva",
      kcdf    = "Gaussian",
      min.sz  = 1,
      mx.diff = TRUE,
      verbose = TRUE
    )

    sample_acts_gsva <- t(gsva_results) %>%
      scale() %>%
      as.data.frame()

    results_list$GSVA <- sample_acts_gsva
  }

  ###### Save outputs if requested
  if (return) {
    if (!is.null(results_list$PROGENy)) {
      utils::write.csv(results_list$PROGENy, paste0("Results/Pathway_matrix_PROGENy_", file.name, ".csv"))
    }
    if (!is.null(results_list$GSVA)) {
      utils::write.csv(results_list$GSVA, paste0("Results/Pathway_matrix_GSVA_", file.name,".csv"))
    }
  }

  ###### RETURN FORMAT LOGIC
  # If only one element → return it directly (matrix)
  if (length(results_list) == 1) {
    return(results_list[[1]])
  }

  # If both → return the list
  return(results_list)
}

#' Compute TF Network Classification
#'
#' Classifies transcription factor (TF) modules into clusters based on their correlations with pathway activity values across samples.
#' Uses hierarchical clustering and silhouette width to determine the optimal number of clusters.
#'
#' @param tf.network A list containing TF module eigengenes, typically output from \code{compute.WTCNA()}
#' @param pathways.features A matrix with pathway activities, typically from \code{compute.pathway.activity()}
#' @param return Logical. If TRUE, intermediate plots (e.g. silhouette, dendrogram, PCA) are saved in the \code{Results/} directory. Default is TRUE.
#'
#' @return A named list of TF module clusters.
#' @export
#'
#' @examples
#'
#' data("network.tuto")
#' data("counts.norm.tuto")
#' pathways <- compute.pathway.activity(counts.norm.tuto)
#' tfs.modules.clusters <- compute.TF.network.classification(tf.network = network.tuto,
#'                                                           pathways.features = pathways,
#'                                                           return = FALSE)
#'
compute.TF.network.classification = function(tf.network, pathways.features, return = T){

  tf.network = data.frame(tf.network[[1]])
  pathways.features = data.frame(pathways.features)

  moduleTraitCor = WGCNA::cor(tf.network, pathways.features, method = "p")
  # moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(tf.network))
  # significance_threshold <- 0.05
  #
  # # Replace non-significant correlations with zero
  # moduleTraitCor[moduleTraitPvalue > significance_threshold] <- 0
  #
  # # Identify columns that have variance (i.e., are not constant or all zeros)
  # non_constant_columns <- apply(moduleTraitCor, 2, function(x) var(x) != 0)
  # moduleTraitCor <- moduleTraitCor[, non_constant_columns] # Filter out the columns that are constant or all zeros

  ### Find clusters
  silhouette = factoextra::fviz_nbclust(moduleTraitCor, hcut, method = "silhouette", k.max = nrow(moduleTraitCor)-1)
  k_cluster = as.numeric(silhouette$data$clusters[which.max(silhouette$data$y)])

  if(return){
    pdf(paste0("Results/TFs_modules_Silhouette_scores"))
    print(silhouette)
    dev.off()
  }

  hc_modules = stats::hclust(dist(moduleTraitCor), method = "ward.D2")
  dend_pathways = as.dendrogram(hc_modules)

  if(return){
    pdf(paste0("Results/TFs_modules_clusters"))
    plot(dend_pathways, cex = 0.6)
    stats::rect.hclust(hc_modules, k = k_cluster, border = 2:5)
    dev.off()
  }

  ### Extract clusters
  sub_grp <- dendextend::cutree(hc_modules, k = k_cluster)

  ### Plot PCA and biplot
  p = factoextra::fviz_cluster(list(data = moduleTraitCor, cluster = sub_grp))

  if(return){
    pdf(paste0("Results/PCA_TFs_modules_clusters"))
    print(p)
    dev.off()
  }

  res.pca <- stats::prcomp(moduleTraitCor,  scale = F)
  p = factoextra::fviz_pca_biplot(res.pca, label="all", select.var = list(contrib = 6), addEllipses=TRUE, ellipse.level=0.75)

  if(return){
    pdf(paste0("Results/PCA_Biplot_TFs_modules_clusters"))
    print(p)
    dev.off()
  }

  # Extract the loadings
  loadings <- res.pca$rotation
  contribution <- (loadings^2)*100
  features = contribution %>%
    data.frame() %>%
    tibble::rownames_to_column("Features") %>%
    dplyr::arrange(desc(PC1))

  if(return){
    pdf("Results/Pathways_contribution_pca.pdf", width = 12, height = 8)
    p = ggplot2::ggplot(features, aes(x = reorder(Features, -PC1), y = PC1)) +
      geom_bar(stat = "identity", fill = "skyblue") +
      labs(title = "Contribution of pathways",
           x = "Feature",
           y = "Contribution (%)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size=15))
    print(p)
    dev.off()
  }

  groups = list()
  for (i in 1:length(unique(sub_grp))) {
    groups[[i]] = names(sub_grp)[sub_grp == i]
    names(groups)[i] = paste0(gsub("ME", "", groups[[i]]), collapse = "_")
  }

  return(groups)
}

#' Compute Transcription Factor (TF) activity
#'
#' Infers transcription factor (TF) activity from a gene expression matrix using the VIPER algorithm (Alvarez et al., 2016). The function requires a TF–target gene regulatory network, which can be provided by the user or obtained from OmnipathR resources such as CollecTRI or Dorothea. ARACNE-inferred networks are also supported.
#'
#' @param RNA.counts A gene expression matrix with genes as rows and samples as columns. The matrix should be normalized (e.g., TPM, log2CPM, etc.).
#' @param TF.collection Character. The source of the TF–target network. Options are `"CollecTRI"` (default), `"Dorothea"`, or `"ARACNE"`.
#' - `"CollecTRI"` and `"Dorothea"` use prebuilt collections from OmnipathR.
#' - `"ARACNE"` allows user input of a custom network file in a 3-column format: `regulator`, `target`, and `mutual information`.
#' @param min_targets_size Integer. Minimum number of target genes per regulon required for TF activity inference. Default is 5.
#' @param universe Optional. A user-specified data frame of TF–target interactions. If not provided, the function will fetch the relevant network based on the `TF.collection` argument.
#' @param cores Integer. Number of cores used by VIPER inference. Default is 4.
#' @param return Logical; if TRUE, saves matrix in Results/ folder. Default is TRUE.
#' @param file.name Optional character suffix used when writing the TF activity matrix to disk.
#'
#' @return A data frame of inferred and scaled TF activity scores, with samples as rows and TFs as columns.
#' @export
#'
#' @references
#' Alvarez, M. et al. (2016). Functional characterization of somatic mutations in cancer using network-based inference of protein activity. *Nature Genetics*, 48(8), 838–847. https://doi.org/10.1038/ng.3593
#'
#' Türei, D., Korcsmáros, T., & Saez-Rodriguez, J. (2016). OmniPath: guidelines and gateway for literature-curated signaling pathway resources. *Nature Methods*, 13(12), 966–967. https://doi.org/10.1038/nmeth.4077
#'
#' Garcia-Alonso, L. et al. (2019). Benchmark and integration of resources for the estimation of human transcription factor activities. *Genome Research*. https://doi.org/10.1101/gr.240663.118
#'
#' Lachmann, A. et al. (2016). ARACNe-AP: gene network reverse engineering through adaptive partitioning inference of mutual information. *Bioinformatics*, 32(14), 2233–2235. https://doi.org/10.1093/bioinformatics/btw216
#'
#' Margolin, A.A. et al. (2006). ARACNE: an algorithm for the reconstruction of gene regulatory networks in a mammalian cellular context. *BMC Bioinformatics*, 7(Suppl 1), S7. https://doi.org/10.1186/1471-2105-7-S1-S7
#'
#' @examples
#' data("counts.norm.tuto")
#' tfs_activity <- compute.TFs.activity(counts.norm.tuto, cores = 1)
#'
compute.TFs.activity <- function(RNA.counts, TF.collection = "CollecTRI", min_targets_size = 5, universe = NULL, cores = 4, return = TRUE, file.name = NULL){

  tf_cache_file <- "Results/TF_target_collection.csv"
  if(TF.collection == "CollecTRI"){
    if(is.null(universe)){
      if(file.exists(tf_cache_file)){
        universe = utils::read.csv(tf_cache_file, row.names = 1)
      } else {
        universe = decoupleR::get_collectri(organism = 'human', split_complexes = F)
        utils::write.csv(universe, tf_cache_file)
      }
    }
  } else if(TF.collection == "Dorothea"){
    if(is.null(universe)){
      if(file.exists(tf_cache_file)){
        universe = utils::read.csv(tf_cache_file, row.names = 1)
      } else {
        universe = dplyr::filter(dorothea::dorothea_hs, .data$confidence %in% c("A", "B")) %>%
          dplyr::mutate(source = .data$tf) %>%
          dplyr::select(-tf)
        utils::write.csv(universe, tf_cache_file)
      }
    }
  }

  if(TF.collection == "ARACNE"){
    cat("For ARACNE analysis you need to specify the path of your network file. Remember this file should be a 3 columns text file, with regulator in the first column, target in the second and mutual information in the third column")
    network_file = readline(prompt = "Path for network file from aracne (no quotes): ")
    net_regulons <- viper::aracne2regulon(network_file, as.matrix(RNA.counts), format = "3col")
  }

  sample_acts <- decoupleR::run_viper(mat = RNA.counts, network = universe,
                                      .source = "source", minsize = min_targets_size, eset.filter = FALSE, method = "none",
                                      verbose = FALSE, cores = cores) %>%
    decoupleR::pivot_wider_profile(id_cols = source,
                                   names_from = condition,
                                   values_from = score) %>%
    as.matrix()

  sample_acts = scale(t(sample_acts)) ## Scale NES TFs to account for variability in the number of targets per TF and make them comparable across modules

  if(return){
    utils::write.csv(sample_acts, paste0("Results/TF_matrix_", file.name, ".csv"))
  }

  return(data.frame(sample_acts))

}

#' Compute Weighted TF-coactivity Network Analysis (WTCNA)
#'
#' Construct a weighted signed or unsigned network using TF activity to cluster protein regulators
#' into modules that share similar activity patterns. Each TF module will have a sample-level score
#' represented by the eigenvalue of the module.
#'
#' @param TFs.matrix Matrix of TF activity (samples x TFs).
#' @param batch Logical; if TRUE, performs consensus WGCNA across cohorts provided as a list.
#' @param network.type Network type: "signed", "unsigned", "signed hybrid", or "distance". Default is "signed".
#' @param clustering.method Clustering method for hierarchical clustering. Default is "ward.D2".
#' @param minMod Minimum number of TFs per module. Default is 15.
#' @param corr_mod Correlation threshold (0–1) for merging similar modules. Default is 0.9.
#' @param cor_type Correlation type for adjacency calculation: "p" (Pearson), "s" (Spearman). Default is "p".
#' @param softPower Optional numeric value specifying the soft-thresholding power to be used when constructing
#' @param verbose Boolen value to whether print or no the function messages
#' @param file.name Optional character suffix used when writing WTCNA outputs.
#' @param return Logical, whether to save output plots and module list to "Results/". Default is TRUE.
#'
#' @return A named list with:
#' \itemize{
#'   \item \code{TFs module matrix}: Matrix of module eigengenes (samples x modules).
#'   \item \code{TFs colors}: Vector of module colors assigned to each TF.
#'   \item \code{TFs per module}: List of TF names in each module.
#'   \item \code{Proportion of variance}: Matrix of variance explained per module.
#' }
#'
#' @references
#' Langfelder, P., & Horvath, S. (2008). WGCNA: an R package for weighted correlation network analysis.
#' BMC Bioinformatics, 9, 559. https://doi.org/10.1186/1471-2105-9-559
#'
#' @export
#'
#' @examples
#'
#' data("tfs.tuto")
#' network <- compute.WTCNA(tfs.tuto, corr_mod = 0.9, clustering.method = "ward.D2", return = FALSE)
compute.WTCNA <- function(TFs.matrix, batch = FALSE, network.type = "signed", clustering.method = "ward.D2",
                          minMod = 15, corr_mod = 0.9, cor_type = "p", verbose = F, file.name = NULL,
                          softPower = NULL, return = T) {

  if(verbose){
    cat("Creating weighted TF-coactivity network......................................................................\n\n")
  }

  if(batch){
    if(!is.list(TFs.matrix) || is.data.frame(TFs.matrix)){
      stop("If batch = TRUE, TFs.matrix must be a list of matrices, one per cohort")
    }
    nCohorts <- length(TFs.matrix)
    if(verbose) cat("Running consensus WGCNA across", nCohorts, "cohorts...\n")

    # Step 1: pick soft-threshold per cohort if not provided
    if(is.null(softPower)){
      powers = c(1:10, 12:20)
      softPower <- numeric(nCohorts)
      for(i in 1:nCohorts){
        if(verbose) cat("Picking soft threshold for cohort", i, "...\n")
        sft <- WGCNA::pickSoftThreshold(TFs.matrix[[i]], powerVector = powers, verbose = 0)
        target = 0.9
        diff = abs(-sign(sft$fitIndices[,3]) * sft$fitIndices[,2] - target)
        softPower[i] <- powers[which.min(diff)]
        if(verbose) cat("Cohort", i, "softPower =", softPower[i], "\n")
      }
    } else {
      if(length(softPower) != nCohorts){
        stop("softPower must be a vector with one value per cohort")
      }
    }

    # Step 2: create multiExpr list
    multiExpr <- lapply(TFs.matrix, function(x) list(data = x))

    # Step 3: run blockwiseConsensusModules
    consensusModules <- WGCNA::blockwiseConsensusModules(
      multiExpr,
      power = softPower,                  # vector of powers per cohort
      networkType = network.type,
      minModuleSize = minMod,
      deepSplit = 2,
      pamRespectsDendro = FALSE,
      mergeCutHeight = 0.25,
      corType = ifelse(cor_type == "p", "pearson", "spearman"),
      verbose = ifelse(verbose, 3, 0)
    )

    dynamicColors <- consensusModules$colors

    # Step 4: Extract module eigengenes per cohort and scale within cohort
    MEs_list <- consensusModules$multiMEs
    scaled_MEs_list <- lapply(MEs_list, function(x) {
      ME <- x$data              # extract the data.frame
      ME_scaled <- scale(ME)    # scale MEs scores to account for variabilities across cohorts even though they are already scaled at TF level (MEs can be different in scale)
      colnames(ME_scaled) <- colnames(ME)
      return(ME_scaled)
    })

    # Step 5: concatenate scaled MEs across cohorts
    MEs <- do.call(rbind, scaled_MEs_list)
    colnames(MEs) <- gsub("ME", "", colnames(MEs))

    # Step 6: Map TFs to modules (same as before)
    modtfs <- list()
    modules <- unique(dynamicColors)
    tfs <- colnames(TFs.matrix[[1]])
    for(i in seq_along(modules)){
      modtfs[[i]] <- tfs[dynamicColors == modules[i]]
    }
    names(modtfs) <- modules

    output <- list(
      "TFs module matrix" = MEs,
      "TFs colors" = dynamicColors,
      "TFs per module" = modtfs,
      "TFs_matrix" = TFs.matrix
    )


  } else {


    #####Choose parameter for scale-free network topology
    powers = c(c(1:10), seq(from = 12, to=20, by=1))
    sink(tempfile())
    invisible(sft <- WGCNA::pickSoftThreshold(TFs.matrix, powerVector = powers, verbose = 0, networkType = network.type))
    sink()

    if(return){
      pdf("Results/Soft_Threshold")
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
      diff = abs(-sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2] - target) #Calculate absolute difference
      min_index = which.min(diff) #Identify the index with the minimum difference
      softPower = powers[min_index]

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

    adjacency = WGCNA::adjacency(TFs.matrix, power =softPower, type=network.type, corFnc = "cor", corOptions = list(use = cor_type))
    TOM = WGCNA::TOMsimilarity(adjacency, TOMType = network.type, verbose = 0)
    dissTOM = 1-TOM

    #####Unsupervised hierarchical clustering using dissimilarity matrix
    geneTree = stats::hclust(dist(dissTOM), method = clustering.method)
    dynamicMods = dynamicTreeCut::cutreeDynamic(dendro = geneTree, distM = dissTOM,
                                                deepSplit = 2, pamRespectsDendro = FALSE,
                                                minClusterSize = minMod, verbose = 0);
    dynamicColors = WGCNA::labels2colors(dynamicMods)

    #####Remove variables and clean garbage
    rm(adjacency, TOM, dissTOM)
    gc()

    if(return){
      pdf("Results/Gene_dendrogram_and_module_colors")
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

    if(verbose){
      print(paste0("Merging modules significantly correlated with ", corr_mod, "........"))
    }

    merge = mergeModules(MEs, dynamicColors, corr_mod)
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
      pdf("Results/Gene_dendrogram_and_module_colors_after_merging")
      WGCNA::plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                                 dendroLabels = FALSE, hang = 0.03,
                                 addGuide = TRUE, guideHang = 0.05,
                                 main = "Gene dendrogram and module colors")
      dev.off()
    }

    TFspropVar = WGCNA::propVarExplained(TFs.matrix, dynamicColors, MEs, corFnc = "cor", corOptions = "use = 'p'")

    colnames(MEs) <- gsub("ME", "", colnames(MEs))

    output = list(MEs, dynamicColors, modtfs, TFspropVar, TFs.matrix)
    names(output) = c("TFs module matrix", "TFs colors", "TFs per module", "Proportion of variance", "TFs_matrix")

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
    utils::write.csv(tfs_modules, paste0('Results/TFs_modules_', file.name,'.csv'), row.names = F)
    utils::write.csv(MEs, file = paste0("Results/TF_module_matrix_", file.name, ".csv"), row.names = F)
  }

  return(output)
}

#' Identify hub TFs
#'
#' Identifies hub TFs per module using values of module membership and degree.
#'
#' @param datExpr A matrix of TF activity (TFs as rows and samples as columns).
#' @param TF.network TF network obtained from compute.WTCNA().
#' @param MM_thresh Threshold for module membership (e.g., 0.8).
#' @param degree_thresh Quantile threshold for degree (e.g., 0.9 for top 10%).
#'
#' TFs with high module membership (r > MM_thresh) and among the top percentage
#' of genes by degree (above degree_thresh quantile) are considered hub TFs.
#'
#' @return A list with two elements:
#' \describe{
#'   \item{hubGenes}{Named list of hub TFs per module}
#'   \item{detailedData}{Dataframe with module, degree, and membership info}
#' }
#' @export
#'
#' @examples
#'
#' data("tfs.tuto")
#' data("network.tuto")
#'
#' hub_tfs <- identify_hub_TFs(t(tfs.tuto), network.tuto, MM_thresh = 0.8, degree_thresh = 0.9)
#'
identify_hub_TFs <- function(datExpr, TF.network, MM_thresh = 0.8, degree_thresh = 0.9) {
  moduleEigengenes_df = TF.network[[1]]
  moduleColors = TF.network[[2]]
  # Calculate Module Membership (MM)
  moduleMemberships <- sapply(unique(moduleColors), function(module) {
    genesInModule <- which(moduleColors == module)
    eigengene <- moduleEigengenes_df[, module]
    stats::cor(t(datExpr[genesInModule, ]), eigengene)
  })

  # Calculate adjacency matrices and degrees
  adjacencyList <- lapply(unique(moduleColors), function(module) {
    genesInModule <- which(moduleColors == module)
    moduleData <- datExpr[genesInModule, ]
    adjacency <- stats::cor(t(moduleData))
    adjacency[lower.tri(adjacency)] <- 0
    adjacency
  })

  names(adjacencyList) <- unique(moduleColors)

  moduleDegrees <- lapply(unique(moduleColors), function(module) {
    adjacency <- adjacencyList[[module]]
    degree <- rowSums(adjacency)
    names(degree) <- rownames(datExpr)[which(moduleColors == module)]
    degree
  })
  names(moduleDegrees) <- unique(moduleColors)

  # Identify hub genes
  hubGenesList <- list()
  hubGenesData <- data.frame()  # Initialize an empty data frame

  allDegrees <- numeric()
  allMemberships <- numeric()

  for (module in unique(moduleColors)) {
    genesInModule <- which(moduleColors == module)
    degrees <- moduleDegrees[[module]]
    memberships <- moduleMemberships[[module]]

    # Update the global lists of degrees and memberships
    allDegrees <- c(allDegrees, degrees)
    allMemberships <- c(allMemberships, memberships)

    # Calculate the cutoff for the top 10% by degree
    degreeCutoff <- stats::quantile(degrees, degree_thresh)

    # Plot distribution of Degrees for the current module
    # degree_plot_data <- data.frame(Degree = degrees)
    #
    # degree_plot <- ggplot(degree_plot_data, aes(x = Degree)) +
    #   geom_histogram(binwidth = 5, fill = module, color = "black", alpha = 0.6) +
    #   geom_vline(xintercept = degreeCutoff, linetype = "dashed", color = "red") +
    #   labs(
    #     title = paste("Distribution of Gene Degrees in Module", module),
    #     x = "Degree",
    #     y = "Frequency"
    #   ) +
    #   theme_minimal()
    #
    # print(degree_plot)

    # Plot distribution of Module Membership for the current module
    # membership_plot_data <- data.frame(ModuleMembership = memberships)
    #
    # membership_plot <- ggplot(membership_plot_data, aes(x = ModuleMembership)) +
    #   geom_histogram(binwidth = 0.05, fill = module, color = "black", alpha = 0.6) +
    #   geom_vline(xintercept = MM_thresh, linetype = "dashed", color = "red") +
    #   labs(
    #     title = paste("Distribution of Module Membership in Module", module),
    #     x = "Module Membership",
    #     y = "Frequency"
    #   ) +
    #   theme_minimal()
    #
    # print(membership_plot)

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

  return(list(hubGenes = hubGenesList, detailedData = hubGenesData))
}

#' Identify cell groups
#'
#' Identifies cell dendrograms corresponding to each TF module group based on correlation with deconvolution features.
#'
#' @param features A list with two named elements:
#'   \describe{
#'     \item{correlations}{A matrix of correlations between TF modules and cell type features.}
#'     \item{significant}{A named list of significant features per TF module (e.g., p-value < 0.05).}
#'   }
#' @param clustering.method Clustering method used in hclust. Default: "ward.D2".
#' @param width Width of the saved PDF plots.
#' @param height Height of the saved PDF plots.
#' @param return Logical; whether to save dendrogram plots to the "Results/" folder.
#'
#' @return A named list of dendrograms for each TF module.
#' @export
#'
identify.cell.groups = function(features, clustering.method = "ward.D2", width = 12, height = 18, return = T){

  moduleTraitCor = features[[1]]

  #names(features[[2]]) = paste0("ME", names(features[[2]])) #To match names of columns from corr matrix

  features_vec = features[[2]]

  #Discard modules with no significant features
  zeros = c()
  for (i in 1:length(features_vec)) {
    if(length(features_vec[[i]])<2){
      zeros = c(i, zeros)
    }
  }

  if(length(zeros)!=0){
    moduleTraitCor = moduleTraitCor[-zeros,,drop = F]
    features_vec = features_vec[-zeros]
  }


  if(length(features_vec) == 0){
    return(NULL)
  }
  #Gather significant features across TF modules clusters
  # features_vec = list()
  # for (i in 1:length(tfs.modules.groups)) {
  #   features_vec[[i]] = unique(unlist(unname(features[[2]][tfs.modules.groups[[i]]])))
  #   names(features_vec)[i] = names(tfs.modules.groups)[i]
  # }

  lis.dendrogram = list()

  for (i in 1:length(features_vec)){
    TFmoduleTraitcor = moduleTraitCor[i,colnames(moduleTraitCor)%in%features_vec[[i]], drop = F]
    #TFmoduleTraitcor = TFmoduleTraitcor[rownames(TFmoduleTraitcor)%in%tfs.modules.groups[[i]], , drop=F]
    ###Dendogram by Module

    #d = 1 - TFmoduleTraitcor #Distance based on correlation (close to 0 means similar correlation per module)
    #d = d/sqrt(nrow(TFmoduleTraitcor)) #Adjust/Scale distance matrix for number of features to make dendrograms comparable
    d <- stats::dist(t(TFmoduleTraitcor))
    dendrogram <- stats::hclust(d, method = clustering.method)
    #dendrogram <- hclust(dist(t(d)), method = clustering.method)
    if(return){
      pdf(paste0("Results/Dendogram_cell_types_", names(features_vec)[i]), width = width, height = height)
      par(mar = c(5, 2, 4, 35)) #bottom, left, top, right
      plot(as.dendrogram(dendrogram), horiz= T)
      dev.off()
    }
    lis.dendrogram[[i]] = dendrogram
  }

  names(lis.dendrogram) = names(features_vec)

  return(lis.dendrogram)

}

#' Compute a cell-type composition matrix from deconvolution subgroups
#'
#' Builds a binary presence matrix indicating which original cell types are present in higher-level cell groups.
#'
#' @param deconvolution.subgroupped A list containing "Deconvolution subgroups composition" per model.
#' @param cell.groups A list containing cell group definitions, where the second element holds the groupings.
#' @param cells_extra Optional vector of additional cell identifiers to consider during extraction.
#'
#' @return A binary matrix (data.frame) where rows are cell groups and columns are cell types (1 = present, 0 = absent).
#' @export
#'
compute.composition.matrix = function(deconvolution.subgroupped, cell.groups, cells_extra = NULL){

  # Initialize composition matrix
  cells_types = extract_cells(colnames(deconvolution.subgroupped[["Deconvolution matrix"]]), cells_extra = cells_extra)
  presence_matrix <- data.frame(matrix(data = 0, nrow = length(cell.groups[[2]]), ncol = length(cells_types)))
  colnames(presence_matrix) <- cells_types

  #Extract deconvolution subgroups composition for each ML model
  deconv_subgroups = list()
  contador = 1
  subgroups = deconvolution.subgroupped[["Deconvolution subgroups composition"]]
  for (i in 1:length(subgroups)) {
    if(length(subgroups[[i]])!=0){ #Whether a specific cell type does not contains subgroups
      for (j in 1:length(subgroups[[i]])) {
        deconv_subgroups[[contador]] = subgroups[[i]][[j]]
        names(deconv_subgroups)[contador] = names(subgroups[[i]])[j]
        contador = contador + 1
      }
    }
  }

  # Update presence_matrix for this model
  row_index = 1
  for (j in seq_along(cell.groups[[2]])) {
    features <- cell.groups[[2]][[j]] #Extract cell group j from ML model i
    for (cell_feature in features) { #Iterate over features in cell group j ML model i
      cells <- get_all_cells(cell_feature, deconv_subgroups) #Use recursive function to extract all nested cell features from different subgroup levels
      cells_types = extract_cells(cells, cells_extra = cells_extra)
      presence_matrix[row_index, cells_types] <- 1 #Set 1 if cell feature is present in subgroup
    }
    row_index <- row_index + 1
  }

  rownames(presence_matrix) = names(cell.groups[[2]])

  idy = which(colSums(presence_matrix)==0)

  if(length(idy)!=0){
    warning("Features with no presence in cell groups:", "\n", paste0(names(idy), collapse = ", "))
    presence_matrix = presence_matrix[,-idy]
  }

  # --- MINIMAL ADDITION: duplicate rows ---
  presence_matrix = presence_matrix[rep(1:nrow(presence_matrix), each = 2), ]

  rownames(presence_matrix) = paste0(rep(names(cell.groups[[2]]), each = 2),
                                     c("_pos","_neg"))

  return(presence_matrix)

}

#' Construct cell groups based on TF networks and deconvolution
#'
#' Identifies and projects cell groups using module relationships derived from TF networks and deconvolution outputs.
#' If a binary trait is specified, the function splits the data and constructs cell groups for both classes (supervised analysis).
#'
#' @param network A list containing TF networks for cell types.
#' @param dt A list containing deconvolution subgroup structures.
#' @param batch Optional vector indicating batch assignment for samples.
#' @param pval Numeric. P-value threshold applied both to filter TF module–deconvolution
#'   feature correlations and as the significance cutoff for the CCA permutation test.
#'   Default: 0.05.
#' @param clustering.method Clustering method for hierarchical clustering. Default: "ward.D2".
#' @param n_perm Integer. Number of permutations for the CCA significance test per cell group.
#'   Higher values give more precise p-values but increase runtime. Default: 999.
#'
#' @return A list of 3 elements:
#' \describe{
#'   \item{scores}{A data frame or matrix with the projected cell group scores (samples x groups).}
#'   \item{composition}{A named list where each element is a character vector of original cell types per group.}
#'   \item{loadings}{A list of numeric vectors indicating the loadings (feature contributions) for each group.}
#' }
#' @export
construct_cell_groups = function(network, dt, batch = NULL, pval = 0.05, clustering.method = "ward.D2", n_perm = 999, dendrogram_file = NULL, return_dendrogram = FALSE){

  corr_modules = compute.modules.relationship(network[[1]], dt[[1]], batch = batch, return = T, plot = F, pval = pval)
  cell_dendrograms = identify.cell.groups(corr_modules, clustering.method = clustering.method, height = 20, return = F)
  cell.groups = cell.groups.computation(dt[[1]], cell_dendrograms, network, batch = batch, return = F, pval = pval, n_perm = n_perm, dendrogram_file = dendrogram_file, return_dendrogram = return_dendrogram)

  names(cell.groups) = c("Cell_groups", "Composition", "Weights")
  return(cell.groups)
}

#' Compute composite scores on test set based on previous cell groups
#'
#' This function simulates cell subgroups compositions from deconvolution results in the test deconvolution data,
#' and calculates composite scores using provided cell groups and features in order to replicate cell groups
#' coming from a training set into an independent set.
#' The composite scores represent summarized information from cell subgroup profiles.
#'
#' @param deconv_res A list containing deconvolution results, including
#'   subgroup compositions (list of data frames or matrices).
#' @param cell_groups A list with three elements:
#'   - cell groups scores
#'   - composition: A named list of character vectors where each element represents
#'     cells belonging to a specific group.
#'   - loadings: A corresponding list of numeric vectors (loadings) for each cell group.
#' @param features A character vector of feature names to select relevant cell groups.
#' @param deconvolution_test A data frame or matrix of deconvolution features for the test set,
#'   with cells as columns and samples as rows.
#'
#' @return A data frame where each column corresponds to a composite score
#'   calculated for each feature group in the test set.
#'   If no composite scores can be computed due to zero variance,
#'   returns an empty data frame with a printed message.
#'
#' @details
#' The function first simulates cell subgroups by computing median values across
#' specified iterations and joins them with the original test deconvolution data.
#' Then it extracts the relevant cells for each feature and calculates composite
#' scores.
#'
#' @export
compute.test.set = function(deconv_res, cell_groups, features, deconvolution_test){

  ################################################################################Simulate cell subgroups
  deconv_subgroups = deconv_res[["Deconvolution subgroups composition"]]
  iterations = find.maximum.iteration(deconv_subgroups)

  ## Extract the deconv feature without the cluster type
  features_with_clusters <- colnames(deconv_res[["Deconvolution matrix"]])
  has_clusters <- grepl("_.*(mixed|immunosuppressive|immunoactive)$", features_with_clusters)

  if(any(has_clusters)){

    # Base name = everything before final cluster label
    base_names <- sub("_(mixed|immunosuppressive|immunoactive)$",
                      "",
                      features_with_clusters)

    # Cluster suffix = cluster type
    cluster_suffixes <- sub(".*_(mixed|immunosuppressive|immunoactive)$",
                            "\\1",
                            features_with_clusters)

    map <- data.frame(base = base_names,
                      suffix = cluster_suffixes,
                      stringsAsFactors = FALSE)
  }

  if(is.infinite(iterations) && iterations < 0){
    warning("No subgroups to replicate")
    if(any(has_clusters)){
      ## Paste the corresponding clusters to the deconvolution features
      colnames(deconvolution_test) <- paste0(colnames(deconvolution_test), map$suffix[match(colnames(deconvolution_test), map$base)])
    }

    deconvolution_test = deconvolution_test[,colnames(deconvolution_test) %in% colnames(deconv_res[["Deconvolution matrix"]])]
  }else{
    # Create same groups composition
    for (m in 1:iterations) {
      base_groups = list()
      for (i in 1:length(deconv_subgroups)){
        if(length(deconv_subgroups[[i]])!=0){
          idy = grep(paste0("Iteration.",m), names(deconv_subgroups[[i]]))
          if(length(idy)!=0){
            base_groups = append(base_groups, deconv_subgroups[[i]][idy])
          }
        }
      }

      deconv_subgroups_values = c()
      for (i in 1:length(base_groups)) {

        x = as.matrix(deconvolution_test[, colnames(deconvolution_test) %in% base_groups[[i]], drop = FALSE])

        if(ncol(x) == 0){
          med = rep(0, nrow(deconvolution_test))
        } else {
          med = matrixStats::rowMedians(x)
        }

        deconv_subgroups_values = cbind(deconv_subgroups_values, med) #Compute median using base groups
      }
      colnames(deconv_subgroups_values) = names(base_groups)
      deconvolution_test = cbind(deconv_subgroups_values, deconvolution_test) # Join cell subgroups and deconv features

    }

    if(any(has_clusters)){
      ## Paste the corresponding clusters to the deconvolution features
      colnames(deconvolution_test) <- paste0(colnames(deconvolution_test), "_", map$suffix[match(colnames(deconvolution_test), map$base)])
    }

    deconvolution_test = deconvolution_test[,colnames(deconvolution_test)%in%colnames(deconv_res[[1]])]
  }

  # Compute composite scores
  idx = which(names(cell_groups[[2]]) %in% features)
  cell_dendrogram = c()
  names = c()
  for (i in 1:length(idx)) {
    cells = cell_groups[[2]][[idx[i]]]
    pca_cells = deconvolution_test[,colnames(deconvolution_test) %in% cells, drop = F]
    pca_cells <- pca_cells[, apply(pca_cells, 2, function(x) all(!is.na(x)) && var(x, na.rm = TRUE) != 0), drop = FALSE] #Drop zero-columns or NAs
    name_cell_group = names(cell_groups[[2]][idx[i]])
    color = stringr::str_split(name_cell_group, "_")[[1]][2]
    loadings_cells = cell_groups[[3]][[idx[i]]]
    cell_dendrogram = cbind(cell_dendrogram, compute.test.score(pca_cells, loadings_cells))
    names = c(names, names(cell_groups[[2]])[idx[i]])
  }

  if(is.null(cell_dendrogram)==T){
    print("No composite scores because all features have zero variance.")
  }else{
    colnames(cell_dendrogram) = names
  }

  return(data.frame(cell_dendrogram))

}

#' Extract cells from cell type groups
#'
#' @param groups A character vector of cell type group names, typically from cell.groups.computation()
#' @param cells_extra Optional character vector of additional cell type names to include in the extraction
#'
#' @return A character vector of unique cell types found in the groups, ignoring method and signature suffixes
#' @export
#'
extract_cells = function(groups, cells_extra = NULL){
  names_cells = c("B.cells", "B.naive", "B.memory", "Macrophages.cells", "Macrophages.M0", "Macrophages.M1", "Macrophages.M2", "Monocytes", "Neutrophils", "NK.cells", "NK.activated",
                  "NK.resting", "NKT.cells", "CD4.cells", "CD4.memory.activated", "CD4.memory.resting", "CD4.naive", "CD8.cells", "T.cells.regulatory", "T.cells.non.regulatory","T.cells.helper",
                  "T.cells.gamma.delta", "Dendritic.cells", "Dendritic.activated", "Dendritic.resting", "Cancer", "Endothelial", "Eosinophils", "Plasma.cells", "Myocytes", "Fibroblasts",
                  "Mast.cells", "Mast.activated", "Mast.resting", "CAF")


  if(is.null(cells_extra) == F){
    names_cells = c(names_cells, cells_extra)
  }

  # Create regex to capture base cell type + cluster
  # Matches known cell type, optionally followed by .Iteration.x or .Subgroup.x, then _mixed/_immunosuppressive/_immunoactive
  regex_pattern <- paste0(
    "(",
    paste(names_cells, collapse = "|"),
    ")",
    "((?:[._](?:Subgroup|Iteration)\\.[0-9]+)*)",
    "_(mixed|immunosuppressive|immunoactive)"
  )


  normalized_names <- sapply(groups, function(x) {
    m <- regexpr(regex_pattern, x, perl = TRUE)
    if (m != -1) {
      matched <- regmatches(x, m)
      # remove Iteration / Subgroup artifacts
      gsub("([._](Subgroup|Iteration)\\.[0-9]+)+", "", matched)
    } else {
      NA
    }
  })

  normalized_names <- unique(na.omit(normalized_names))

  return(normalized_names)
}

#' Extract colors
#'
#' Extract TF module colors from cell type group names
#'
#' @param module_colors Character vector of TF module colors, e.g., from compute.WTCNA()
#' @param cell_group_name Character vector or string of cell type group names to search within
#'
#' @return Character vector of matched module colors found in cell_group_name, ordered by their appearance.
#' Returns NA if no matches are found.
#' @export
#'
extract_colors <- function(module_colors, cell_group_name) {
  matches <- c() # For storing the matches
  for (color in module_colors) {
    match <- regexpr(color, cell_group_name) # Find the position of the match
    if (match != -1) {
      matches <- c(matches, regmatches(cell_group_name, match))
    }
  }


  if (length(matches) > 0) {
    order <- sapply(matches, function(m) regexpr(m, cell_group_name)) # Sort matches based on their position in the original string to ensure names are the same
    matches <- matches[order(order)]  # Order the matches based on their position
    return(matches)  # Return the ordered matches
  } else {
    return(NA)  # If no matches are found, return NA
  }

}

#' Create TFs modules
#'
#' This function re-create existing TF modules on a different TF activity matrix.
#'
#' @param TF.matrix TFs activity matrix with samples as rows and TFs as columns (should be the output of compute.TF.activity()).
#' @param network_tfs A TF network object obtained from compute.WTCNA() from which TF modules need to be re-create.
#'
#' @return A matrix of TFs modules scores across samples
#' @export
#'
create_tfs_modules = function(TF.matrix, network_tfs){

  TF.matrix = TF.matrix[, colnames(TF.matrix) %in% colnames(network_tfs$TFs_matrix)]

  tfs.modules = TF.matrix %>%
    t() %>%
    data.frame() %>%
    dplyr::mutate(Module = "na") ## Create column to assign the corresponding module to each TF

  for (i in 1:length(network_tfs[[3]])) {
    tfs.modules$Module[which(rownames(tfs.modules) %in% network_tfs[[3]][[i]])] = names(network_tfs[[3]])[i]
  }

  tfs_colors = tfs.modules %>%
    dplyr::pull(Module)

  MEList = WGCNA::moduleEigengenes(TF.matrix, colors = tfs_colors, scale = F) #Data already scale
  MEs = MEList$eigengenes
  MEs = WGCNA::orderMEs(MEs)

  colnames(MEs) <- gsub("ME", "", colnames(MEs))

  return(MEs)
}

#' Find maximum iteration from cell subgroups
#'
#' Scans a nested list of cell subgroups and returns the highest iteration number
#' found across all elements.
#'
#' @param cells.groups A nested list of cell subgroups. Names of inner elements
#'   must follow the pattern `*.Iteration.<n>` where `<n>` is an integer.
#'
#' @return An integer giving the maximum iteration number across all subgroups.
#'
#' @keywords internal
find.maximum.iteration = function(cells.groups){
  max_iteration = c()
  for (i in 1:length(cells.groups)){
    if(is.null(names(cells.groups[[i]]))==F){
      iterations <- sapply(names(cells.groups[[i]]), function(x) {
        as.numeric(sub(".*\\.Iteration\\.(\\d+)", "\\1", x))
      })
      local_max = max(unlist(iterations))
      max_iteration = c(max_iteration, local_max)
    }
  }

  return(max(max_iteration))
}

#' Merge highly correlated TF modules
#'
#' Identifies pairs of TF modules whose eigengene correlation exceeds \code{corr}
#' and merges them by averaging their columns.
#'
#' @param data A numeric matrix or data frame of TF module eigengenes (samples x modules).
#' @param colors A character vector of module color labels aligned with the columns of \code{data}.
#' @param corr Numeric. Spearman correlation threshold above which two modules are merged. Default 0.9.
#'
#' @return A list of two elements:
#' \itemize{
#'   \item \code{[[1]]}: Data frame of merged module eigengenes (samples x modules).
#'   \item \code{[[2]]}: Updated character vector of module color labels after merging.
#' }
#'
#' @keywords internal
mergeModules = function(data, colors, corr){
  df = correlation(data)
  idx = which(round(df$r,2) > corr)
  if(length(idx)>0){
    for(i in seq(1, length(idx), by=2)){
      module1 = df$measure1[idx[i]]
      module2 = df$measure2[idx[i]]
      if((module1 %in% colnames(data)) && (module2 %in% colnames(data))){
        colors[which(colors%in%c(substring(module1, 3), substring(module2, 3)))] = substring(module1, 3)
        data <- data %>%
          dplyr::mutate(new_column = rowMeans(dplyr::select(., dplyr::all_of(c(module1, module2))))) %>%
          dplyr::rename(!!module1 := new_column) %>%
          dplyr::select(., -dplyr::all_of(module2))
      }
    }
  }

  return(list(data, colors))
}

#' Remove cell groups with duplicate composition
#'
#' Identifies cell groups whose cell-type composition is identical to another group
#' and removes the duplicates, keeping only the first occurrence.
#'
#' @param cell.values A list of numeric vectors of cell group scores.
#' @param cell.composition A list of character vectors describing cell-type membership per group.
#' @param cell.loadings A list of loading vectors corresponding to each cell group.
#'
#' @return A list of three elements:
#' \itemize{
#'   \item \code{[[1]]}: Deduplicated cell group scores.
#'   \item \code{[[2]]}: Deduplicated cell group compositions.
#'   \item \code{[[3]]}: Deduplicated cell group loadings.
#' }
#'
#' @keywords internal
remove_equal = function(cell.values, cell.composition, cell.loadings){

  #Sorted list to avoid no recognizing vectors with equal composition but different order of cells
  for(i in 1:length(cell.composition)){
    for (j in 1:length(cell.composition[[i]])) {
      cell.composition[[i]][[j]] = sort(cell.composition[[i]][[j]])
    }
  }

  #Remove cell groups
  for(i in 1:length(cell.composition)){
    exist = c() #Initialize vector of equalities
    rang = seq(1, length(cell.composition))[-i] #Create sequence to iterate all list elements except the one being analyzed
    for (j in rang){
      idx = which(cell.composition[[i]] %in% cell.composition[[j]] == TRUE) #Map all cell groups which already existed
      exist = c(exist, idx) #Save index cluster
    }
    if(length(exist)!=0){
      cell.composition[[i]] = cell.composition[[i]][-unique(exist)] #Remove cell groups that already exist
      cell.values[[i]] = cell.values[[i]][-unique(exist)] #Remove cell groups that already exist
      cell.loadings[[i]] = cell.loadings[[i]][-unique(exist)] #Remove cell loadings that already exist
    }
  }

  #Remove dendrograms without elements (length equal 0)
  vec = c()
  for(i in 1:length(cell.composition)){
    if(length(cell.composition[[i]]) == 0){
      vec = c(vec, i)
    }
  }

  if(length(vec)>0){
    cell.composition = cell.composition[-vec]
    cell.values = cell.values[-vec]
    cell.loadings = cell.loadings[-vec]
  }

  cell.composition = unlist(cell.composition, recursive = FALSE)
  cell.values = unlist(cell.values, recursive = FALSE)
  cell.loadings = unlist(cell.loadings, recursive = FALSE)

  return(list(cell.values, cell.composition, cell.loadings))
}

#' Remove cell groups composed of a single cell type
#'
#' Filters out cell groups whose composition contains only one cell type, as these
#' groups lack multi-cellular context.
#'
#' @param cell.values A list of numeric vectors of cell group scores.
#' @param cell.composition A list of character vectors describing cell-type membership per group.
#' @param cell.loadings A list of loading vectors corresponding to each cell group.
#'
#' @return A list of three elements (scores, compositions, loadings) with singleton groups
#'   removed, or \code{NULL} if all groups are removed.
#'
#' @keywords internal
remove_single_groups = function(cell.values, cell.composition, cell.loadings){

  vec = c()
  for (i in 1:length(cell.composition)) {
    if(length(cell.composition[[i]])==1){
      vec = c(vec, i)
    }
  }

  if(length(vec)>0){
    cell.composition = cell.composition[-vec]
    cell.values = cell.values[-vec]
    cell.loadings = cell.loadings[-vec]
  }

  if(length(cell.composition)==0){
    return(NULL)
  }else{
    return(list(cell.values, cell.composition, cell.loadings))
  }

}

#' Calculate dendrogram cut heights
#'
#' Computes a sequence of candidate cut heights for each cell-type dendrogram.
#' Heights are distributed between a buffered minimum and maximum derived from
#' the dendrogram's own height distribution, avoiding trivial cuts (single-element
#' or all-in-one clusters).
#'
#' @param cell.group.dendrogram A list of \code{hclust} objects, one per TF module,
#'   as returned by \code{identify.cell.groups()}.
#' @param n_cuts Integer. Number of evenly spaced cut heights to generate per
#'   dendrogram. If \code{NULL} (default), the number is set proportional to
#'   the maximum dendrogram height.
#'
#' @return A list of numeric vectors, one per dendrogram, containing the
#'   candidate cut heights.
#'
#' @keywords internal
calculate_dendrogram_cuts = function(cell.group.dendrogram, deep_split = 4, min_cluster_size = 3){

  cuts = list()
  for (i in 1:length(cell.group.dendrogram)) {
    hc <- as.hclust(as.dendrogram(cell.group.dendrogram[[i]]))

    labels <- dynamicTreeCut::cutreeDynamic(
      dendro         = hc,
      method         = "tree",
      deepSplit      = deep_split,
      minClusterSize = min_cluster_size,
      verbose        = 0
    )

    cuts[[i]] <- labels
  }

  return(cuts)

}

plot_dendrogram_clusters = function(cell.group.dendrogram, cuts_per_dendrogram, file_name = NULL) {
  plots <- list()

  for (i in 1:length(cell.group.dendrogram)) {
    dend   <- as.dendrogram(cell.group.dendrogram[[i]])
    labels <- cuts_per_dendrogram[[i]]
    names(labels) <- labels(dend)

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

    p <- dendextend::as.ggdend(dend) |>
      ggplot2::ggplot() +
      ggplot2::coord_flip() +
      ggplot2::scale_y_reverse() +
      ggplot2::ggtitle(paste("Dendrogram:", names(cell.group.dendrogram)[i])) +
      ggplot2::theme_classic() +
      ggplot2::theme(axis.text   = ggplot2::element_blank(),
                     axis.ticks  = ggplot2::element_blank(),
                     plot.title  = ggplot2::element_text(hjust = 0.5, face = "bold")) +
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

#' Remove highly correlated cell groups
#'
#' Computes pairwise Spearman correlations among cell group score vectors and
#' removes one member of each pair whose absolute correlation exceeds
#' \code{threshold}.
#'
#' @param data A list of three elements:
#'   \describe{
#'     \item{scores}{A numeric data frame or matrix of cell group scores (samples x groups).}
#'     \item{compositions}{A named list of cell-type vectors describing group membership.}
#'     \item{loadings}{A named list of loading vectors corresponding to each cell group.}
#'   }
#' @param threshold Numeric. Correlation threshold above which one of a correlated
#'   pair is removed. Default is 0.95.
#'
#' @return A list of three elements (scores, compositions, loadings) with
#'   redundant cell groups removed.
#'
#' @keywords internal
remove.cell.groups.corr <- function(data, threshold = 0.95) {

  # Compute correlation matrix
  corr_matrix <- stats::cor(data[[1]])
  # Find highly correlated features
  contador = 1
  while(nrow(corr_matrix)>0){
    color_features = c()
    feature = data.frame(corr_matrix[1, , drop = FALSE]) #Extract first row feature
    feature_corr = feature %>%                                #Take only high corr above threshold
      dplyr::mutate_all(~ifelse(. > threshold, ., NA)) %>%
      dplyr::select_if(~all(!is.na(.)))
    feature = feature_corr #In order to save original corr matrix and print names
    corr_matrix = corr_matrix[-which(rownames(corr_matrix)%in%colnames(feature)),-which(colnames(corr_matrix)%in%colnames(feature)), drop = F] #Remove already joined features
    color_features = list()
    if(ncol(feature)>1){
      for (m in 1:ncol(feature)) {
        name_cell_group = colnames(feature)[m]
        color = stringr::str_split(name_cell_group, "_")[[1]][2]
        color_features[[m]] = color
      }

      len = length(unique(sapply(unname(data[[2]][colnames(feature)]), length)))
      if(len != 1){
        feature = feature[,which.min(sapply(unname(data[[2]][colnames(feature)]), length)), drop = F] #Remove elements with higher number of features and keep the min composition that explains the score (high corr > 0.95 of two different compositions implies that extra cell types are just noise)
      }
      new_group_composition = unique(unlist(unname(data[[2]][colnames(feature)])))
      new_group_value = rowMeans(data[[1]][,colnames(data[[1]])%in%colnames(feature),drop=F])

      ################ Merging cell loadings

      # Select the loadings matrices corresponding to the current set of correlated features
      mats <- data[[3]][colnames(data[[1]]) %in% colnames(feature)]

      # Ensure that all matrices have proper column names (set the rownames() as colnames() as they are square)
      mats_fixed <- lapply(mats, function(x) {
        if (is.null(colnames(x))) {
          colnames(x) <- rownames(x)
        }
        x
      })

      # Compute the union of all features across the selected matrices: full set of features that the merged matrix should contain
      all_features <- Reduce(union, lapply(mats_fixed, rownames))

      # Align each matrix to the full set of features
      # - Add missing rows/columns filled with 0 for features not present in the matrix
      # - Reorder rows and columns according to 'all_features'
      mats_aligned <- lapply(mats_fixed, function(x) {
        missing <- setdiff(all_features, rownames(x))  # Identify features missing in this matrix
        if(length(missing) > 0){
          # Add missing rows (initialized to 0)
          x <- rbind(x, matrix(0, nrow = length(missing), ncol = ncol(x),
                               dimnames = list(missing, colnames(x))))
          # Add missing columns (initialized to 0)
          x <- cbind(x, matrix(0, nrow = nrow(x), ncol = length(missing),
                               dimnames = list(rownames(x), missing)))
        }
        # Reorder rows and columns so all matrices have the same order
        x[all_features, all_features, drop = FALSE]
      })

      # # Ensure that all matrices have proper column names (set the rownames() as colnames() as they are square)
      # mats_fixed <- lapply(mats, function(x) {
      #
      #   # Convert vector or scalar -> square matrix
      #   if (is.null(dim(x))) {
      #     # Give names if missing
      #     if (is.null(names(x)) || any(names(x) == "")) {
      #       names(x) <- paste0("Feature_", seq_along(x))
      #     }
      #     x <- matrix(x,
      #                 nrow = length(x),
      #                 ncol = length(x),
      #                 dimnames = list(names(x), names(x)))
      #   }
      #
      #   # Ensure rownames exist
      #   if (is.null(rownames(x)) || any(rownames(x) == "")) {
      #     rownames(x) <- paste0("Feature_", seq_len(nrow(x)))
      #   }
      #
      #   # Ensure colnames exist
      #   if (is.null(colnames(x)) || any(colnames(x) == "")) {
      #     colnames(x) <- rownames(x)
      #   }
      #
      #   x
      # })
      #
      #
      # # Compute the union of all features across the selected matrices: full set of features that the merged matrix should contain
      # all_features <- Reduce(union, lapply(mats_fixed, rownames))
      #
      # # Align each matrix to the full set of features
      # # - Add missing rows/columns filled with 0 for features not present in the matrix
      # # - Reorder rows and columns according to 'all_features'
      # mats_aligned <- lapply(mats_fixed, function(x) {
      #   missing <- setdiff(all_features, rownames(x))  # Identify features missing in this matrix
      #   if(length(missing) > 0){
      #     # Add missing rows (initialized to 0)
      #     x <- rbind(x, matrix(0, nrow = length(missing), ncol = ncol(x),
      #                          dimnames = list(missing, colnames(x))))
      #     # Add missing columns (initialized to 0)
      #     x <- cbind(x, matrix(0, nrow = nrow(x), ncol = length(missing),
      #                          dimnames = list(rownames(x), missing)))
      #   }
      #   # Reorder rows and columns so all matrices have the same order
      #   x[all_features, all_features, drop = FALSE]
      # })

      # Combine all aligned matrices by computing the element-wise average
      new_loadings_value <- Reduce("+", mats_aligned) / length(mats_aligned)

      ################

      if(contador==1){
        #Remove features from original data
        new_data <- data[[1]][, -which(colnames(data[[1]])%in%colnames(feature_corr)), drop = F]
        new_groups = data[[2]][-which(names(data[[2]]) %in% colnames(feature_corr))]
        new_loadings = data[[3]][-which(names(data[[2]]) %in% colnames(feature_corr))] #Using names from data[[2]] cause data[[3]] has the same order of elements
      }else{
        new_loadings = new_loadings[-which(names(new_groups) %in% colnames(feature_corr))]
        new_data <- new_data[, -which(colnames(new_data)%in%colnames(feature_corr)), drop = F]
        new_groups = new_groups[-which(names(new_groups) %in% colnames(feature_corr))]
      }

      #Add new combined features
      class <- unique(stringr::str_extract(colnames(feature_corr), "positive|negative"))
      if(length(class)>1){
        class = paste0(class, collapse = ".")
      }

      ## If corr features are from two different classes don't combine, just discard them (as they don't make any distinction between classes for prediction)
      if(is.na(class) == F){
        new_name = paste0("Dendrogram_",  paste0(unique(unlist(color_features)), collapse = "_"), ".group_combined_", contador, "_", class)
      }else{
        new_name = paste0("Dendrogram_",  paste0(unique(unlist(color_features)), collapse = "_"), ".group_combined_", contador)
      }
      new_data = cbind(new_data, new_group_value)
      colnames(new_data)[length(new_data)] = new_name

      new_groups[[length(new_groups)+1]] = new_group_composition
      names(new_groups)[length(new_groups)] = new_name

      new_loadings[[length(new_loadings)+1]] = new_loadings_value
      names(new_loadings)[length(new_loadings)] = new_name

      contador = contador + 1

    }else{ #Nothing is combined
      if(contador == 1){
        new_data = data[[1]]
        new_groups = data[[2]]
        new_loadings = data[[3]]
      }
    }
  }

  res = list(new_data, new_groups, new_loadings)

  return(res)
}

#' Run Reactome pathway enrichment for a single TF module
#'
#' Given the hub TFs of a module, extracts their target genes from the TF-gene
#' universe, keeps the most variable targets, and runs Reactome over-representation
#' analysis using the full gene expression matrix as the background universe.
#'
#' @param tpm.counts A numeric matrix of normalized expression values (genes x samples).
#' @param module_color Character. Name of the TF module (color label) to analyze.
#' @param hub_genes A list as returned by \code{identify_hub_TFs()}, where the first
#'   element maps module names to vectors of hub TF gene symbols.
#' @param tfs_universe A data frame of TF-target interactions with at minimum columns
#'   \code{source} (TF) and \code{target} (gene).
#'
#' @return A \code{ReactomePA} enrichResult object restricted to pathways with
#'   adjusted p-value < 0.05, or \code{NULL} if no enrichment is found.
#'
#' @keywords internal
module_enrich = function(tpm.counts, module_color, hub_genes, tfs_universe){
  # genes = colnames(TFs.matrix)
  # inModule = is.finite(match(module_colors,module))
  # modGenes = genes[inModule]
  targets = tfs_universe[tfs_universe$source %in% hub_genes[[1]][[module_color]],] #Extract targets from TFs
  targets = unique(targets$target) #Keep only unique targets from TFs

  targets_genes = tpm.counts[rownames(tpm.counts)%in%targets,] #Extract gene expression from targets
  targets_genes = targets_genes[order(matrixStats::rowVars(targets_genes), decreasing = T),][1:round(0.2*nrow(targets_genes)),] #Keep only highly variable targets (20%)

  entrz <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = rownames(targets_genes), columns = "ENTREZID", keytype = "SYMBOL") #Change to EntrezID
  universe = AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = rownames(tpm.counts), columns = "ENTREZID", keytype = "SYMBOL") #Change to EntrezID

  reac <- ReactomePA::enrichPathway(gene    = entrz$ENTREZID,
                                    organism     = 'human',
                                    universe = universe$ENTREZID,
                                    pvalueCutoff = 0.05)

  reac@result = reac@result[reac@result$p.adjust<0.05,]

  if(nrow(reac@result)!=0){
    return(reac)
  }

}

#' Permutation test for the first canonical correlation
#'
#' Estimates the null distribution of the first canonical correlation by
#' row-permuting X and recomputing CCA \code{n_perm} times, then returns the
#' empirical p-value (proportion of null correlations >= observed).
#'
#' @param X Numeric matrix (samples x features), already scaled.
#' @param Y Numeric matrix (samples x features), already scaled.
#' @param n_perm Integer. Number of permutations. Default 999.
#'
#' @return A list with \code{r} (observed first canonical correlation) and
#'   \code{p_value} (empirical one-sided p-value).
#'
#' @keywords internal
permutation_cca_test <- function(X, Y, n_perm = 999) {
  obs_r <- stats::cancor(X, Y)$cor[1]
  null_r <- replicate(n_perm, {
    X_perm <- X[sample(nrow(X)), , drop = FALSE]
    stats::cancor(X_perm, Y)$cor[1]
  })
  list(r = obs_r, p_value = mean(null_r >= obs_r))
}

#' Compute composite score for cell groups
#'
#' Computes a composite score by performing Canonical Correlation Analysis (CCA)
#' between cell group features and corresponding TF module scores.
#'
#' @param cell_group A numeric matrix of cell deconvolution features for a cell group (samples x features).
#' @param module_group A character vector indicating TF module group colors corresponding to the cell group (can be obtained via `extract_colors()`).
#' @param tfs.module.network Output of compute.WTCNA().
#' @param batch Optional vector indicating batch assignment for samples.
#' @param discard Logical; whether to discard cell groups that do not pass the
#'   permutation test for the first canonical correlation (default TRUE).
#' @param pval Numeric. Significance threshold for the permutation test (default 0.05).
#' @param n_perm Integer. Number of permutations used to build the null distribution (default 999).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{selected_components}: Numeric matrix of the first canonical component scores across samples.
#'   \item \code{xcoef}: The canonical weights (coefficients) for the cell group features.
#' }
#' If discarded due to low correlation, returns \code{list("NA", "NA")}.
#'
#' @export
#'
compute_composite_score = function(cell_group, module_group, tfs.module.network, batch = NULL, discard = T, pval = 0.05, n_perm = 999){

  modules = tfs.module.network[["TFs module matrix"]]
  tfs_all = tfs.module.network[["TFs_matrix"]]

  if(!is.null(batch) & is.list(tfs_all)){
    tfs_all <- do.call(rbind, tfs_all)
  }

  ### Only in case there are groups combined
  if(length(grep(".group", module_group))==1){
    module_group = gsub(".group", "", module_group)
  }

  tf_per_module_matrix = tfs_all[, colnames(tfs_all) %in% tfs.module.network[["TFs per module"]][[module_group]]]
  tf_module_matrix = modules[, grep(module_group, colnames(modules)), drop = F]

  # ---- Regress out batch if provided ----
  if(!is.null(batch)){
    if(!is.numeric(batch)) batch <- as.numeric(as.factor(batch))  # Convert to numeric if needed
    cell_group <- apply(cell_group, 2, function(x) residuals(lm(x ~ batch)))
    tf_module_matrix <- apply(tf_module_matrix, 2, function(x) residuals(lm(x ~ batch)))
    tf_per_module_matrix <- apply(tf_per_module_matrix, 2, function(x) residuals(lm(x ~ batch)))
  }

  ### We take the TF module only to verify the overall sign of association, following analysis is done with full TF matrix per color
  signs = sign(stats::cor(cell_group, tf_module_matrix))

  if(signs[1] < 0){ ## Because groups come from hierarchical clustering you expect they will not have multiple signs so you just take the first assuming is the same for the others
    cell_group_matrix = scale(cell_group*-1) #Multiply by -1 because CCA identifies the most positive correlated components so to identify the most negative correlated components we need to invert matrix
  }else{
    cell_group_matrix = scale(cell_group)
  }

  # CCA analysis: permutation test against the TF module eigengene to validate the association
  # before committing to the full per-TF CCA projection.
  if(discard){
    perm <- permutation_cca_test(cell_group_matrix, tf_module_matrix, n_perm = n_perm)
    if(perm$p_value >= pval){
      return(list("NA", "NA"))
    }
  }

  #Once we confirmed the potential association we can further do a CCA to find a linear component that maximizes the linear relationship. If we dont do this, we will invent linear components that maximizes correlation but that were not existing before

  # CCA analysis: Note 'cell_group_matrix' is already scale
  cca_result <- stats::cancor(cell_group_matrix, tf_per_module_matrix) ## Now this is done per tfs_per_module (tfs are already scale check compute.TFs.activity())

  # Apply canonical weights to each matrix
  cell_group <- cell_group[, rownames(cca_result$xcoef), drop = F] #Ensure no features have been discard on the way and if yes, subset the matrix (CCA can discard collinear or zero variance features)

  scaled_obj <- scale(as.matrix(cell_group))
  train_means <- attr(scaled_obj, "scaled:center")
  train_sds   <- attr(scaled_obj, "scaled:scale")

  ### xcoef[,1] corresponds to the most correlated linear component (we scale cause the coef came from the scale matrix)
  weighted_cell_group_matrix <- scaled_obj %*% cca_result$xcoef[, 1] #Multiply by the original matrix even if the coefx came from the inverse matrix because we need to find the inverse relationship

  ### MIGHT BE USEFUL AFTER TO DISCARD VARIABLES THAT DONT HELP TO THE ASSOCIATION AND REDUCE GROUP COMPOSITION

  # weights
  # x_weights <- cca_result$xcoef[, 1]
  #y_weights <- cca_result$ycoef[, 1]

  # # canonical scores (variates)
  # U1 <- as.vector(scale(cell_group_matrix) %*% x_weights)
  # V1 <- as.vector(scale(tf_module_matrix) %*% y_weights)
  #
  #
  # # canonical loadings (= correlations of original vars with their own variate)
  # x_loadings <- as.numeric(cor(as.matrix(cell_group_matrix), U1))   # length = ncol(cell_group_matrix)
  # y_loadings <- as.numeric(cor(as.matrix(tf_module_matrix), V1))   # length = ncol(tf_module_matrix)
  #
  # # cross-loadings (= correlation of X with V1 and Y with U1)
  # x_cross_loadings <- as.numeric(cor(as.matrix(cell_group_matrix), V1))
  # y_cross_loadings <- as.numeric(cor(as.matrix(tf_module_matrix), U1))
  #
  # x_df <- data.frame(variable = colnames(cell_group_matrix),
  #                    weight = x_weights,
  #                    loading = x_loadings,
  #                    cross_loading = x_cross_loadings)
  # x_df <- x_df[order(-abs(x_df$loading)), ]   # order by importance
  projection_params = list(
    xcoef       = cca_result$xcoef[, 1, drop = F],
    train_means = train_means,
    train_sds   = train_sds
  )

  return(list(weighted_cell_group_matrix, projection_params))

}

#' Extract significant features using Wilcoxon test
#'
#' Performs Wilcoxon rank-sum test for each feature comparing groups defined by the trait.
#'
#' @param features A numeric data frame or matrix where columns are features and rows are samples.
#' @param trait A vector or factor defining group labels for each sample.
#'
#' @return A character vector of feature names with significant difference between trait groups (p < 0.05).
#'
#' @export
#'
extract_wilcox_significant = function(features, trait){
  significant_features <- c()
  data = cbind(trait, features)
  colnames(data)[1] = "trait"
  for (feature in colnames(data)[-1]) {  # Exclude the first column (trait)
    test_result <- stats::wilcox.test(data[,feature] ~ data$trait)

    if (test_result$p.value < 0.05) {
      significant_features <- c(significant_features, feature)
    }
  }

  return(significant_features)
}

#' Classify samples by high or low deconvolution values in given cell groups
#'
#' Classifies samples in `coldata` as "High" or "Low" based on whether their deconvolution
#' values in all specified cell groups exceed the median.
#'
#' @param coldata A data frame with sample metadata.
#' @param deconvolution A numeric matrix/data frame with cell deconvolution features (samples x cell groups).
#' @param group A character vector specifying the cell groups (columns) to consider.
#'
#' @return The input `coldata` with an additional factor column `Cells_level` ("High" or "Low").
#'
#' @export
#'
classify.deconvolution = function(coldata, deconvolution, group){
  deconv = deconvolution[,colnames(deconvolution)%in%group]

  #Patients high in group 1 of cells
  vec = c()
  if(is.null(ncol(deconv))==T){
    idx = which(deconv > median(deconv))
    vec = c(vec,idx)
  }else{
    for (i in 1:ncol(deconv)) {
      idx = which(deconv[,i] > stats::median(deconv[,i]))
      vec = c(vec, idx)
    }
  }

  #High in all deconv features from group
  pos = which(table(vec) == length(group))

  coldata = coldata %>%
    dplyr::mutate(Cells_level = "Low")

  coldata$Cells_level[pos] = "High"
  coldata$Cells_level = factor(coldata$Cells_level)

  return(coldata)

}

#' Perform pairwise Spearman correlation across all features
#'
#' @param data A numeric matrix or data frame where columns are features.
#'
#' @return A data frame of pairwise significant correlations (p < 0.05), with
#'   columns \code{measure1}, \code{measure2}, \code{r}, \code{p}, \code{sig_p},
#'   \code{p_if_sig}, \code{r_if_sig}, and \code{AbsR}, ordered by descending \code{r}.
#'
#' @keywords internal
correlation <- function(data) {

  M <- Hmisc::rcorr(as.matrix(data), type = "spearman")
  Mdf <- purrr::map(M[c("r", "P", "n")], ~data.frame(.x))

  corr_df = Mdf %>%
    purrr::map(~tibble::rownames_to_column(.x, var="measure1")) %>%
    purrr::map(~tidyr::pivot_longer(.x, -measure1, names_to = "measure2")) %>%
    dplyr::bind_rows(.id = "id") %>%
    tidyr::pivot_wider(names_from = id, values_from = value) %>%
    dplyr::rename(p = P) %>%
    dplyr::mutate(sig_p = ifelse(p < .05, T, F),
                  p_if_sig = ifelse(sig_p, p, NA),
                  r_if_sig = ifelse(sig_p, r, NA))

  corr_df = stats::na.omit(corr_df)  #remove the ones that are the same features (pval = NA)
  corr_df <- corr_df[which(corr_df$sig_p==T),]  #remove not significant
  corr_df <- corr_df[order(corr_df$r, decreasing = T),]
  corr_df$AbsR =  abs(corr_df$r)

  return(corr_df)

}

#' Recursively retrieve all base cell types for a subgroup
#'
#' Walks a nested subgroup hierarchy and returns the leaf-level cell type names
#' that belong to the requested subgroup.
#'
#' @param subgroup_name Character. Name of the subgroup or base cell type to resolve.
#' @param cell_subgroups A named list mapping subgroup names to vectors of member
#'   names (which can themselves be subgroup keys).
#'
#' @return A character vector of unique base cell type names belonging to
#'   \code{subgroup_name}.
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

#' Compute mean silhouette width for a clustering
#'
#' @param clusters An integer vector of cluster assignments (one per sample).
#' @param distance_matrix A numeric matrix used to compute Euclidean distances
#'   between samples.
#'
#' @return A single numeric value: the mean silhouette width across all samples.
#'
#' @keywords internal
compute_silhouette <- function(clusters, distance_matrix) {
  silhouette_scores <- cluster::silhouette(clusters, stats::dist(distance_matrix))
  mean(silhouette_scores[, 3])  # Return average silhouette width
}

#' Project test-set cell group scores using training CCA parameters
#'
#' Scales a test-set cell group matrix using the mean and standard deviation
#' stored from training, then projects it onto the first canonical component
#' learned during training.
#'
#' @param cell_group A numeric matrix of cell deconvolution features for test
#'   samples (samples x features). Column names must overlap with training features.
#' @param projection_params A list containing:
#'   \describe{
#'     \item{xcoef}{Named numeric matrix of CCA canonical weights (features x 1).}
#'     \item{train_means}{Named numeric vector of training column means.}
#'     \item{train_sds}{Named numeric vector of training column standard deviations.}
#'   }
#'
#' @return A numeric matrix (samples x 1) of projected canonical scores for the
#'   test cohort, or \code{NULL} if no features overlap.
#'
#' @keywords internal
compute.test.score = function(cell_group, projection_params){

  xcoef       = projection_params$xcoef
  train_means = projection_params$train_means
  train_sds   = projection_params$train_sds

  # Subset to common features
  common_features = intersect(colnames(cell_group), rownames(xcoef))
  if(length(common_features) == 0){
    warning("No features overlap between test data and training loadings.")
    return(NULL)
  }

  cell_group  = cell_group[, common_features, drop = F]
  xcoef       = xcoef[common_features, , drop = F]
  train_means = train_means[common_features]
  train_sds   = train_sds[common_features]

  # Scale using training means/SDs (unflipped, matching training)
  cell_group_scaled = scale(cell_group, center = train_means, scale = train_sds)

  return(as.matrix(cell_group_scaled) %*% xcoef)
}

#' Prepare CellTFusion folds for cross-validation with training and test data
#'
#' This function prepares data for k-fold cross-validation using the CellTFusion framework.
#' It supports two modes of operation:
#' \enumerate{
#'   \item If a set of tuned hyperparameters is provided via \code{bestune}, the function will
#'   process the full dataset once with those parameters and return a single processed training set.
#'   \item If no tuned parameters are provided, the function will construct fold-specific training
#'   and test datasets, expanding over a grid of CellTFusion hyperparameters and returning them
#'   for subsequent model training and hyperparameter selection.
#' }
#'
#' For each fold, training and test sets are generated by running the \code{CellTFusion()} pipeline
#' on the training data (gene expression and optional deconvolution features). The trained
#' projection is then applied to the test data to ensure comparability between folds.
#'
#' @param data A data frame containing gene expression data (samples x genes) and a column named \code{target}
#'   indicating class labels.
#' @param folds A list of integer vectors indicating row indices for the training set in each fold.
#'   The test set is implicitly defined as the complement.
#' @param deconv A matrix or data frame of deconvolution features (samples x features). If \code{NULL},
#'   CellTFusion will run without these features.
#' @param universe Optional universe of features for CellTFusion.
#' @param paths Optional list of prior knowledge resources for CellTFusion.
#' @param normalized Logical. Whether the gene expression data is already normalized. Defaults to \code{FALSE}.
#' @param coldata A data frame with metadata (e.g., sample annotations), must match the number and order of samples in \code{data}.
#' @param time_var Optional numeric vector with survival/censoring time values.
#' @param event_var Optional vector with event labels used for survival tasks.
#' @param corr_type Correlation type passed to CellTFusion. Default is \code{"spearman"}.
#' @param ncores Integer. Number of CPU cores to use for parallelization. If \code{NULL}, \code{parallel::detectCores() - 2} is used.
#' @param batch Logical; whether batch correction should be used. Default is \code{FALSE}.
#' @param batch_id Optional vector of batch identifiers, aligned with samples.
#' @param min_targets_size Hyperparameter for CellTFusion: minimum TF regulon size.
#' @param minMod Hyperparameter for CellTFusion: minimum WGCNA module size.
#' @param corr_mod Hyperparameter for CellTFusion: module merging correlation threshold.
#' @param corr Hyperparameter for CellTFusion: deconvolution correlation threshold.
#'   These can be provided as vectors to define a tuning grid when \code{bestune = NULL}.
#' @param bestune Optional. A data frame of tuned hyperparameters. If provided, the function skips fold construction
#'   and directly processes the full dataset using these values.
#'
#' @return A list with two possible structures depending on \code{bestune}:
#' \itemize{
#'   \item If \code{bestune} is provided:
#'     \itemize{
#'       \item \code{train_cell_data_final}: Final processed cell group feature matrix for the full dataset, including the \code{target} column.
#'       \item \code{custom_output}: The full CellTFusion object from training.
#'       \item \code{best_celltfusion_params}: The tuned hyperparameters used.
#'     }
#'   \item If \code{bestune} is \code{NULL}:
#'     \itemize{
#'       \item \code{processed_folds}: A list of folds. Each fold contains:
#'         \itemize{
#'           \item \code{train_data}: Processed training data with cell group features and \code{target}.
#'           \item \code{test_data}: Projected test data in the learned feature space.
#'           \item \code{obs_test}: True class labels for the test set.
#'           \item \code{rowIndex}: Row indices corresponding to the test samples.
#'           \item \code{fold_name}: Fold identifier (if provided).
#'           \item \code{params}: Hyperparameter values used for this fold.
#'         }
#'     }
#' }
#'
#' @details
#' When \code{bestune = NULL}, the function expands all possible hyperparameter combinations
#' via \code{expand.grid()} and runs CellTFusion separately for each fold and parameter set.
#' This can be computationally expensive, so the function supports parallelization via
#' \pkg{foreach}, \pkg{doParallel}, and \pkg{parallel}.
#'
#' The results of this function are designed to be passed directly to \code{compute_custom_k_fold_CV()}
#' for hyperparameter selection and model evaluation.
#'
#' @importFrom dplyr mutate select
#' @importFrom stats setNames
#' @importFrom parallel makeCluster stopCluster detectCores
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach foreach %dopar%
#' @export
#'
prepare_CellTFusion_folds <- function(data, folds = NULL, deconv = NULL, universe = NULL, paths = NULL, normalized = FALSE,
                                      coldata, time_var = NULL, event_var = NULL, corr_type = "spearman",
                                      ncores = NULL, batch = F, batch_id = NULL, min_targets_size, minMod, corr_mod, corr, bestune = NULL){

    if(!is.null(bestune)){
      # Run CellTFusion on the full training set
      obs_train <- NULL

      if ("target" %in% colnames(data)) {
        obs_train <- data$target
        data$target <- NULL
      } else if (!is.null(time_var) && !is.null(event_var)) {
        obs_train <- list(time = time_var, event = as.numeric(event_var == trait.positive))
      } else {
        stop("Either 'target' or both 'time_var' and 'event_var' must be provided.")
      }

      required_cols <- c("min_targets_size", "minMod", "corr_mod", "corr")

      best_celltfusion_params <- if (is.data.frame(bestune)) {
        if (all(required_cols %in% names(bestune))) {
          dplyr::select(bestune, dplyr::all_of(required_cols))
        } else {
          message("No tunable parameters found in bestune; using fixed parameters.")
          tibble::tibble(
            min_targets_size = min_targets_size,
            minMod           = minMod,
            corr_mod         = corr_mod,
            corr             = corr
          )
        }
      } else if (is.list(bestune)) {
        if (all(required_cols %in% names(bestune))) {
          tibble::as_tibble(bestune[required_cols])
        } else {
          message("No tunable parameters found in bestune; using fixed parameters.")
          tibble::tibble(
            min_targets_size = min_targets_size,
            minMod           = minMod,
            corr_mod         = corr_mod,
            corr             = corr
          )
        }
      } else {
        stop("`bestune` must be a data.frame or list.")
      }

      train_processed_final <- CellTFusion(
        t(data),
        deconv = deconv,
        normalized = normalized,
        coldata = coldata,
        universe = universe,
        corr_type = corr_type,
        paths = paths,
        min_targets_size = best_celltfusion_params$min_targets_size,
        minMod           = best_celltfusion_params$minMod,
        corr_mod         = best_celltfusion_params$corr_mod,
        corr             = best_celltfusion_params$corr,
        batch = batch, batch_id = batch_id,
        return = FALSE,
        verbose = FALSE
      )

      # Get cell group features
      train_cell_data_final <- train_processed_final$Latent_spaces$Z %>%
        data.frame()

      # train_cell_data_final <- train_processed_final$Cell_groups[[1]] %>%
      #   data.frame()

      if (is.list(obs_train)) {
        train_cell_data_final <- train_cell_data_final %>%
          dplyr::mutate(
            time  = obs_train$time,
            event = obs_train$event
          )
      } else {
        train_cell_data_final <- train_cell_data_final %>%
          dplyr::mutate(target = obs_train)
      }

      custom_output = train_processed_final

      res = list(train_cell_data_final, custom_output, best_celltfusion_params)

      return(res)
    }else{ ### Means best tune need to be find

      custom_grid <- expand.grid(
        min_targets_size = min_targets_size,
        minMod           = minMod,
        corr_mod         = corr_mod,
        corr             = corr
      )

      # Parallelize over folds
      processed_folds <- lapply(seq_along(folds), function(i) {
        cat("Starting fold", names(folds)[i], "\n")

        train_idx <- folds[[i]]
        test_idx <- setdiff(seq_len(nrow(data)), train_idx)

        ## Subset data
        train_data <- data[train_idx, , drop = FALSE]
        train_deconv <- deconv[train_idx, , drop = FALSE]
        train_coldata <- coldata[train_idx, , drop = FALSE]

        # Determine target or time/event
        if ("target" %in% colnames(train_data)) {
          obs_train <- train_data$target
          train_data$target <- NULL

        } else if (!is.null(time_var) && !is.null(event_var)) {
          obs_train <- list(
            time  = time_var[train_idx],
            event = as.numeric(event_var[train_idx] == trait.positive)
          )
        } else {
          stop("Either 'target' or both 'time_var' and 'event_var' must be provided.")
        }

        fold_results <- lapply(seq_len(nrow(custom_grid)), function(j) {
          message("Running fold ", names(folds)[i], ", grid ", j, " / ", nrow(custom_grid))

          train_processed <- CellTFusion(
            raw.counts = t(train_data),
            deconv = train_deconv,
            normalized = normalized,
            coldata = train_coldata,
            universe = universe,
            corr_type = corr_type,
            paths = paths,
            min_targets_size = custom_grid$min_targets_size[j],
            minMod           = custom_grid$minMod[j],
            corr_mod         = custom_grid$corr_mod[j],
            corr             = custom_grid$corr[j],
            batch = batch, batch_id = batch_id,
            return = FALSE,
            verbose = FALSE
          )

          train_cell_data <- train_processed$Latent_spaces$Z %>%
            data.frame()

          # train_cell_data <- train_processed$Cell_groups[[1]] %>%
          #   data.frame()

          if (is.list(obs_train)) {
            train_cell_data <- train_cell_data %>%
              dplyr::mutate(
                time = obs_train$time,
                event = obs_train$event
              )
          } else {
            train_cell_data <- train_cell_data %>%
              dplyr::mutate(target = obs_train)
          }

          # Prepare test data
          test_deconv <- deconv[test_idx, , drop = FALSE]
          obs_test <- if ("target" %in% colnames(data)) {
            data$target[test_idx]
          } else if (!is.null(time_var) && !is.null(event_var)) {
            list(time = time_var[test_idx], event = as.numeric(event_var[test_idx] == trait.positive))
          } else {
            NULL
          }

          test_data <- project_test_factors(train_processed, test_deconv)

          if (is.list(obs_test)) {
            test_data <- test_data %>%
              dplyr::mutate(
                time  = obs_test$time,
                event = obs_test$event
              )
          }

          list(
            train_data = train_cell_data,
            test_data = test_data,
            obs_test = obs_test,
            rowIndex = test_idx,
            fold_name = names(folds)[i],
            params = custom_grid[j, ]
          )
        })

        if (nrow(custom_grid) == 1) {
          fold_results <- fold_results[[1]]
        }

        filename = file.path("Results", paste0("fold_", names(folds)[i], ".rds"))
        saveRDS(fold_results, file = filename)
      })
    }
}

#' Unregister a parallel backend registered with doParallel
#'
#' Switches the foreach backend back to sequential execution and calls
#' \code{gc()} to release memory held by the parallel workers.
#'
#' @return Called for its side effect; returns \code{NULL} invisibly.
#'
#' @keywords internal
unregister_dopar <- function() {
  if (!is.null(foreach::getDoParRegistered())) {
    # switch back to sequential backend
    foreach::registerDoSEQ()
    gc()
  }
}

#' Student's t-test for cell group comparisons
#'
#' Performs a Student’s t-test comparing cell group scores between two groups of a binary trait.
#' Significant features are plotted as boxplots and saved as PDF files in the "Results/" directory.
#'
#' @param cell.groups A list where the first element is a data frame or matrix of cell group scores,
#'        and the second and third elements contain metadata or identifiers for the cell groups.
#' @param coldata A data frame containing sample-level annotations including the trait to test.
#' @param trait Character. Name of the column in `coldata` used as the grouping variable.
#' @param pval Numeric. P-value threshold for significance (default = 0.05).
#'
#' @return A list containing only significant cell groups after the t-test.
#'         Returns \code{NULL} if no significant groups are found.
#' @export
#'
cell.groups.ttest <- function(cell.groups, coldata, trait, pval = 0.05) {
  sig = c()
  coldata[, trait] = as.factor(coldata[, trait])

  if (length(unique(coldata[, trait])) != 2) {
    stop("T-test requires a binary trait (exactly two groups).")
  }

  for (j in 1:ncol(cell.groups[[1]])) {
    data = data.frame(Value = cell.groups[[1]][, j], Trait = coldata[, trait])

    res.ttest <- stats::t.test(Value ~ Trait, data = data, var.equal = FALSE)

    if (round(res.ttest$p.value, 5) <= pval) {
      cat("Significant p-value after T-test for", colnames(cell.groups[[1]])[j], "\n")

      pdf(paste0("Results/Ttest_", trait, "_", colnames(cell.groups[[1]])[j], ".pdf"),
          width = 12, height = 9)
      print(
        ggplot2::ggplot(data, aes(x = Trait, y = Value, fill = Trait)) +
          geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "gray30") +
          geom_jitter(width = 0.15, size = 2.5, alpha = 0.8, color = "black") +
          stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "white") +
          scale_fill_brewer(palette = "Set2") +
          labs(
            title = paste0("Student's t-test - ", colnames(cell.groups[[1]])[j]),
            subtitle = paste0("P-value: ", round(res.ttest$p.value, 5)),
            x = paste0("Clinical trait: ", trait),
            y = "Cell group score"
          ) +
          theme_minimal(base_size = 16) +
          theme(
            axis.text.x = element_text(size = 20),
            axis.title.x = element_text(size = 20),
            legend.title = element_text(size = 15),
            legend.text = element_text(size = 12),
            axis.title.y = element_text(size = 20, angle = 90),
            plot.title = element_text(size = 20),
            plot.subtitle = element_text(size = 15, face = "bold")
          ) +
          scale_fill_discrete(name = trait)
      )
      dev.off()

      sig = c(sig, j)
    }
  }

  if (length(sig) == 0) {
    message("No significant cell groups (p-value < ", pval, ") after t-test.")
    return(NULL)
  } else {
    cell.groups.sig = list()
    cell.groups.sig[[1]] = cell.groups[[1]][, sig, drop = FALSE]
    cell.groups.sig[[2]] = cell.groups[[2]][sig]
    cell.groups.sig[[3]] = cell.groups[[3]][sig]
    return(cell.groups.sig)
  }
}

#' Kruskal–Wallis test for multi-group comparisons
#'
#' Performs a Kruskal–Wallis test to compare cell group scores across multiple trait levels.
#' Significant results are visualized as annotated boxplots with Dunn post-hoc tests.
#'
#' @param cell.groups A list containing cell group score data, metadata, and identifiers.
#' @param coldata A data frame containing sample annotations, including the grouping trait.
#' @param trait Character. Name of the column in `coldata` used as the grouping variable.
#' @param pval Numeric. P-value threshold for significance (default = 0.05).
#'
#' @return A list containing only significant cell groups after Kruskal–Wallis test.
#'         Returns \code{NULL} if no significant groups are found.
#' @export
#'
cell.groups.kruskal.test <- function(cell.groups, coldata, trait, pval = 0.05) {
  sig = c()
  coldata[, trait] = as.factor(coldata[, trait])

  for (j in 1:ncol(cell.groups[[1]])) {
    data = data.frame(Value = cell.groups[[1]][, j], Trait = coldata[, trait])

    # Perform Kruskal–Wallis test
    res.kruskal <- rstatix::kruskal_test(Value ~ Trait, data = data)

    if (res.kruskal$p < pval) {
      cat("Significant p-value after Kruskal–Wallis test for", colnames(cell.groups[[1]])[j], "\n")

      # Dunn post-hoc test + positions for significance annotation
      pwc <- data %>%
        rstatix::dunn_test(Value ~ Trait, p.adjust.method = "BH") %>%
        rstatix::add_xy_position(x = "Trait")

      pdf(paste0("Results/Kruskal_", trait, "_", colnames(cell.groups[[1]])[j], ".pdf"),
          width = 12, height = 9)
      print(
        ggpubr::ggboxplot(data, x = "Trait", y = "Value", fill = "Trait", add = "jitter") +
          ggpubr::stat_pvalue_manual(pwc, hide.ns = TRUE) +
          labs(
            title = "Kruskal–Wallis Test with Dunn Post-hoc",
            subtitle = rstatix::get_test_label(res.kruskal, detailed = TRUE),
            caption = rstatix::get_pwc_label(pwc),
            x = paste0("Clinical trait: ", trait),
            y = paste0("Values for ", colnames(cell.groups[[1]])[j])
          ) +
          theme(
            axis.text.x = element_text(size = 20),
            axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20, angle = 90),
            plot.title = element_text(size = 20),
            plot.subtitle = element_text(size = 15, face = "bold"),
            legend.position = "none"
          )
      )
      dev.off()

      sig = c(sig, j)
    }
  }

  if (length(sig) == 0) {
    message("No significant cell groups (p-value < ", pval, ") after Kruskal–Wallis test.")
    return(NULL)
  } else {
    cell.groups.sig = list()
    cell.groups.sig[[1]] = cell.groups[[1]][, sig, drop = FALSE]
    cell.groups.sig[[2]] = cell.groups[[2]][sig]
    cell.groups.sig[[3]] = cell.groups[[3]][sig]
    return(cell.groups.sig)
  }
}

#' Wilcoxon rank-sum test for binary traits
#'
#' Performs a Wilcoxon rank-sum (Mann–Whitney U) test comparing cell group scores
#' between two levels of a binary clinical trait.
#' Significant features are plotted as boxplots and saved to the "Results/" folder.
#'
#' @param cell.groups A list containing cell group scores and associated metadata.
#' @param coldata A data frame containing sample annotations and clinical traits.
#' @param trait Character. Name of the column in `coldata` used as the binary grouping variable.
#' @param pval Numeric. P-value threshold for significance (default = 0.05).
#'
#' @return A list containing significant cell groups or \code{NULL} if none are found.
#' @export
#'
cell.groups.wilcox.test <- function(cell.groups, coldata, trait, pval = 0.05) {
  sig = c()
  coldata[, trait] = as.factor(coldata[, trait])

  if (length(unique(coldata[, trait])) != 2) {
    stop("Wilcoxon test requires a binary trait (exactly two groups).")
  }

  for (j in 1:ncol(cell.groups[[1]])) {
    data = data.frame(Value = cell.groups[[1]][, j], Trait = coldata[, trait])

    # Perform Wilcoxon rank-sum test
    res.test <- stats::wilcox.test(Value ~ Trait, data = data, exact = FALSE)

    if (round(res.test$p.value, 5) <= pval) {
      cat("Significant p-value after Wilcoxon test for", colnames(cell.groups[[1]])[j], "\n")

      # Save boxplot
      pdf(paste0("Results/Wilcoxon_", trait, "_", colnames(cell.groups[[1]])[j], ".pdf"),
          width = 12, height = 9)
      print(
        ggplot2::ggplot(data, aes(x = Trait, y = Value, fill = Trait)) +
          geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "gray30") +  # hide boxplot outliers
          geom_jitter(width = 0.15, size = 2.5, alpha = 0.8, color = "black") +  # show all sample points
          stat_summary(fun = median, geom = "point", shape = 23,
                       size = 3, fill = "white") +
          scale_fill_brewer(palette = "Set2") +
          labs(
            title = paste0("Wilcoxon Test - ", colnames(cell.groups[[1]])[j]),
            subtitle = paste0("P-value: ", round(res.test$p.value, 5)),
            x = paste0("Clinical trait: ", trait),
            y = "Cell group score"
          ) +
          theme_minimal(base_size = 16) +
          theme(
            axis.text.x = element_text(size = 20),
            axis.title.x = element_text(size = 20),
            legend.title = element_text(size = 15),
            legend.text = element_text(size = 12),
            axis.title.y = element_text(size = 20, angle = 90),
            plot.title = element_text(size = 20),
            plot.subtitle = element_text(size = 15, face = "bold")
          ) +
          scale_fill_discrete(name = trait)
      )
      dev.off()

      sig = c(sig, j)
    }
  }

  # Collect significant cell groups
  if (length(sig) == 0) {
    message("No significant cell groups (p-value < ", pval, ") after Wilcoxon test.")
    return(invisible(NULL))
  } else {
    cell.groups.sig = list()
    cell.groups.sig[[1]] = cell.groups[[1]][, sig, drop = FALSE]
    cell.groups.sig[[2]] = cell.groups[[2]][sig]
    cell.groups.sig[[3]] = cell.groups[[3]][sig]
    return(cell.groups.sig)
  }
}

#' One-way ANOVA test for multi-group comparisons
#'
#' Performs one-way ANOVA to test for differences in cell group scores across multiple levels of a trait.
#' Tukey post-hoc tests are used to identify pairwise differences and significance is visualized as annotated boxplots.
#'
#' @param cell.groups A list containing cell group score data, metadata, and identifiers.
#' @param coldata A data frame containing sample annotations including the grouping variable.
#' @param trait Character. Name of the column in `coldata` used for the grouping variable.
#' @param pval Numeric. P-value threshold for significance (default = 0.05).
#'
#' @return A list of significant cell groups or \code{NULL} if none are significant.
#' @export
#'
cell.groups.anova.test = function(cell.groups, coldata, trait, pval = 0.05){
  sig = c()
  coldata[, trait] = as.factor(coldata[, trait])

  for (j in 1:ncol(cell.groups[[1]])) {
    data = data.frame(Value = cell.groups[[1]][, j], Trait = coldata[, trait])
    res.aov <- rstatix::anova_test(data = data, dv = Value, between = Trait)

    if (res.aov$p < pval) {
      cat("Significant p-value after ANOVA for", colnames(cell.groups[[1]])[j], "\n")

      # Tukey post-hoc test + position for significance annotation
      pwc <- data %>%
        rstatix::tukey_hsd(Value ~ Trait) %>%
        rstatix::add_xy_position(x = "Trait")

      # Generate annotated boxplot (same style as first function)
      pdf(paste0("Results/ANOVA_", trait, "_", colnames(cell.groups[[1]])[j], ".pdf"),
          width = 12, height = 9)
      print(
        ggpubr::ggboxplot(data, x = "Trait", y = "Value", fill = "Trait", add = "jitter") +
          ggpubr::stat_pvalue_manual(pwc, hide.ns = TRUE) +
          labs(
            title = "One-way ANOVA with Tukey HSD",
            subtitle = rstatix::get_test_label(res.aov, detailed = TRUE),
            caption = rstatix::get_pwc_label(pwc),
            x = paste0("Clinical trait: ", trait),
            y = paste0("Values for ", colnames(cell.groups[[1]])[j])
          ) +
          theme(
            axis.text.x = element_text(size = 20, angle = 0),
            axis.title.x = element_text(size = 20),
            legend.position = "none",
            axis.title.y = element_text(size = 20, angle = 90),
            plot.title = element_text(size = 20),
            plot.subtitle = element_text(size = 15, face = "bold")
          )
      )
      dev.off()

      sig = c(sig, j)
    }
  }

  cell.groups.sig = list()
  if (length(sig) == 0) {
    message("No significant cell groups (p-value < ", pval, ") after ANOVA test")
    return(NULL)
  } else {
    cell.groups.sig[[1]] = cell.groups[[1]][, sig, drop = FALSE]
    cell.groups.sig[[2]] = cell.groups[[2]][sig]
    cell.groups.sig[[3]] = cell.groups[[3]][sig]
    return(cell.groups.sig)
  }
}

#' Fisher test using cell groups scores
#'
#' @param cell.groups A list where the first element is a data frame of cell group scores,
#'        and the second element contains metadata or labels for these groups.
#' @param coldata A data frame containing the clinical or experimental traits.
#' @param trait Character. Name of the column in `coldata` to test with Fisher's exact test.
#' @param pval Numeric. P-value threshold for significance (default 0.05).
#'
#' @return A list containing the significant cell groups after Fisher test. Additionally,
#'         it saves corresponding barplot visualizations in the "Results/" folder.
#' @export
#'
cell.groups.fisher.test = function(cell.groups, coldata, trait, pval = 0.05){

  sig = c()
  coldata[,trait] = as.factor(coldata[,trait])

  for (j in 1:ncol(cell.groups[[1]])) {
    coldata = coldata %>%
      dplyr::mutate(level = cell.groups[[1]][,j],
                    Cells_level = ifelse(level > summary(level)[3], 'High', 'Low'))

    contingency = table(coldata[,"Cells_level"], coldata[,trait])
    test = stats::fisher.test(contingency)

    ##Extract only significant features
    if(round(test$p.value, 5) <= pval){
      cat("Significant pval after doing Fisher test for", colnames(cell.groups[[1]])[j], "\n")
      df = data.frame("Cells_level" = coldata[,"Cells_level"], "Trait" = coldata[,trait])
      pdf(paste0("Results/Fisher_", trait, "_", colnames(cell.groups[[1]])[j], ".pdf"), width = 12, height = 9)
      print(ggstatsplot::ggbarstats(df, Cells_level, Trait, results.subtitle = F,
                                    title= paste0("Dendrogram_", colnames(cell.groups[[1]])[j]),
                                    subtitle = paste0("Fisher's exact test, p-value = ", ifelse(test$p.value < 0.001, "< 0.001", round(test$p.value, 5))))+
              ggplot2::theme(plot.title = ggplot2::element_text(size=15), axis.text = ggplot2::element_text(size=14), legend.title = ggplot2::element_text(size=14)))
      dev.off()

      sig = c(j, sig)
    }
  }

  cell.groups.sig = list()
  cell.groups.sig[[1]] = cell.groups[[1]][,sig]
  cell.groups.sig[[2]] = cell.groups[[2]][sig]
  cell.groups.sig[[3]] = cell.groups[[3]][sig]

  if(length(sig)==0){
    message("No significant cell groups (pvalue < ", pval, ") after Fisher test")
  }else{
    return(cell.groups.sig)
  }

}

#' Perform statistical analysis on cell group scores using a specified test
#'
#' This function provides a unified interface to perform one of several
#' statistical tests (Fisher’s exact, Wilcoxon rank-sum, ANOVA, Kruskal–Wallis,
#' or Student’s t-test) on cell group scores in relation to a given clinical or
#' experimental trait. The choice of test is specified by the `method` argument.
#' Each test identifies significant associations and saves a corresponding
#' visualization for each significant feature in the "Results/" directory.
#'
#' @param cell.groups A list where:
#'   \itemize{
#'     \item The first element is a data frame or matrix containing cell group scores
#'           (rows = samples, columns = cell groups).
#'     \item The second and third elements contain corresponding metadata or annotations
#'           for these groups (e.g., names, features, etc.).
#'   }
#' @param coldata A data frame containing clinical or experimental metadata for samples.
#'   Must include the column specified in `trait`.
#' @param trait Character. The name of the column in `coldata` representing the clinical
#'   or experimental trait to test against (e.g., response, subtype, etc.).
#' @param method Character. Statistical test to perform. One of:
#'   \itemize{
#'     \item `"fisher"` — Fisher’s exact test (for categorical data)
#'     \item `"wilcox"` — Wilcoxon rank-sum test (non-parametric, binary traits)
#'     \item `"anova"` — One-way ANOVA (parametric, >2 groups)
#'     \item `"kruskal"` — Kruskal–Wallis test (non-parametric, >2 groups)
#'     \item `"ttest"` — Student’s t-test (parametric, binary traits)
#'   }
#'   Defaults to all available options, but only one can be used per call.
#' @param pval Numeric. P-value threshold for significance (default: 0.05).
#'
#' @details
#' The function automatically calls the corresponding statistical test function
#' based on the `method` argument:
#' \itemize{
#'   \item `cell.groups.fisher.test()`
#'   \item `cell.groups.wilcox.test()`
#'   \item `cell.groups.anova.test()`
#'   \item `cell.groups.kruskal.test()`
#'   \item `cell.groups.ttest()`
#' }
#'
#' Each test produces both a statistical result and visual outputs (PDF plots)
#' stored in the `"Results/"` folder. These visualizations include the relevant
#' test results (p-values) annotated on the plots.
#'
#' @return A list of significant cell groups, where:
#'   \itemize{
#'     \item The first element contains the subset of the original score matrix
#'           for significant features.
#'     \item The second and third elements contain associated metadata or feature
#'           annotations.
#'   }
#'   Returns `NULL` if no significant features are found.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' sig.groups <- cell.groups.stat.analysis(cell.groups, coldata, trait = "response",
#'                                         method = "kruskal", pval = 0.05)
#' }
#'
#' @export
#'
cell.groups.stat.analysis <- function(cell.groups, coldata, trait,
                                      method = c("fisher", "wilcox", "anova", "kruskal", "ttest"),
                                      pval = 0.05) {
  method <- match.arg(method)

  message("Running ", toupper(method), " test for cell group comparison...\n")

  result <- switch(method,
                   fisher   = cell.groups.fisher.test(cell.groups, coldata, trait, pval),
                   wilcox   = cell.groups.wilcox.test(cell.groups, coldata, trait, pval),
                   anova    = cell.groups.anova.test(cell.groups, coldata, trait, pval),
                   kruskal  = cell.groups.kruskal.test(cell.groups, coldata, trait, pval),
                   ttest    = cell.groups.ttest(cell.groups, coldata, trait, pval))

  if (is.null(result)) {
    message("No significant features found using ", method, " test (p < ", pval, ").")
  } else {
    message("Significant features found: ", ncol(result[[1]]))
  }

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
#'   via elbow on reconstruction MSE across ranks 2:8.
#' @param seed Random seed. Default 123.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{Z}: Sample-level NMF factor scores (samples x rank). Non-negative.
#'   \item \code{W}: Feature weights per factor (2*n_CGs x rank). Non-negative.
#'   \item \code{nmf_input}: The positive-negative split matrix fed to NMF (samples x 2*n_CGs).
#'   \item \code{rank}: The rank used.
#' }
#'
#' @details
#' Signed CCA scores are decomposed as:
#'   score_pos = max(score, 0)  -- patient aligned with TF program
#'   score_neg = max(-score, 0) -- patient anti-aligned with TF program
#' Both are concatenated column-wise before NMF. Column names are suffixed
#' with "_pos" and "_neg" to track direction.
#'
#' @export
compute.latent_factors <- function(X, rank = NULL, seed = 123, file_name = NULL, return = TRUE) {

  if (!requireNamespace("RcppML", quietly = TRUE))
    stop("Install RcppML: install.packages('RcppML')")

  set.seed(seed)
  X <- as.matrix(X)

  # ── Positive-negative split ───────────────────────────────────────────────
  X_pos <- pmax(X, 0)
  X_neg <- pmax(-X, 0)
  colnames(X_pos) <- paste0(colnames(X), "_pos")
  colnames(X_neg) <- paste0(colnames(X), "_neg")
  nmf_input <- cbind(X_pos, X_neg)

  A <- t(nmf_input)   # features x samples (RcppML convention)

  # ── Rank estimation ───────────────────────────────────────────────────────
  if (is.null(rank)) {
    message("Estimating optimal rank via elbow on reconstruction MSE (ranks 2:8)...")
    .mse <- function(A, fit) {
      recon <- fit$w %*% diag(fit$d, nrow = length(fit$d)) %*% fit$h
      mean((A - recon)^2)
    }
    mse_vals <- vapply(2:8, function(k) {
      fit <- RcppML::nmf(A, k = k, seed = seed, verbose = FALSE)
      .mse(A, fit)
    }, numeric(1))
    names(mse_vals) <- as.character(2:8)
    d2   <- diff(diff(mse_vals))
    rank <- as.integer(names(mse_vals)[which.min(d2) + 1L])
    message("MSE by rank: ", paste(sprintf("k%s=%.3g", names(mse_vals), mse_vals), collapse = "  "))
    message("Selected rank (elbow): ", rank)
  }

  # ── NMF fit ───────────────────────────────────────────────────────────────
  message("Running RcppML NMF with rank = ", rank, "...")
  nmf_res <- RcppML::nmf(A, k = rank, seed = seed, verbose = FALSE)
  W <- as.matrix(nmf_res$w)    # features x rank
  Z <- t(as.matrix(nmf_res$h)) # samples x rank

  rownames(Z) <- rownames(X)
  colnames(Z) <- paste0("Factor", seq_len(rank))
  colnames(W) <- paste0("Factor", seq_len(rank))
  rownames(W) <- colnames(nmf_input)

  # ── Patient mixture plot ───────────────────────────────────────────────────

  # normalise Z rows to sum to 1
  Z_norm <- sweep(Z, 1, rowSums(Z), "/")
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
      subtitle = "Each bar = one patient  ·  Height = fraction of each factor",
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
#' For each NMF factor (column of the W basis matrix), selects the features
#' whose weight exceeds the specified quantile threshold.
#'
#' @param latent_factors A list returned by \code{compute.latent_factors()},
#'   containing at minimum a matrix \code{W} (features x factors) of NMF basis weights.
#' @param quantile_cutoff Numeric between 0 and 1. Features with weight above
#'   this quantile of the factor's positive weights are retained. Default 0.7.
#'
#' @return A named list (one element per factor) of named numeric vectors, where
#'   names are feature names and values are their NMF basis weights, sorted
#'   descending.
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
#' For each NMF factor, extracts the top-contributing cell groups (by basis
#' weight), maps them to their constituent cell types, and keeps only cell
#' types enriched relative to the background composition. Saves a star-network
#' PDF per factor and returns the weighted cell-type associations.
#'
#' @param latent_factors A list returned by \code{compute.latent_factors()},
#'   containing \code{W} (features x factors) and the NMF model object.
#' @param dt A named list of deconvolution subgroup results, used to build
#'   the cell-group composition matrix via \code{compute.composition.matrix()}.
#' @param cell.groups A list of cell group definitions as returned by
#'   \code{construct_cell_groups()}.
#' @param enrich_thresh Numeric. Minimum enrichment ratio (foreground / background
#'   frequency) for a cell type to be retained per factor. Default 1.5.
#' @param quantile_cutoff Numeric between 0 and 1. Quantile threshold for
#'   selecting top-contributing cell groups per factor. Default 0.7.
#' @param cells_extra Optional character vector of additional cell-type columns
#'   to include in the composition matrix.
#' @param return Logical. If \code{TRUE}, saves network PDF plots to
#'   \code{Results/}. Default \code{TRUE}.
#' @param file_name Character. Suffix appended to output file names.
#'
#' @return A named list (one element per factor) of named numeric vectors giving
#'   the enriched cell types and their cumulative NMF edge weights, sorted
#'   descending.
#'
#' @keywords internal
compute_cells_niches <- function(latent_factors, dt, cell.groups,
                                 enrich_thresh   = 1.5,
                                 quantile_cutoff = 0.7,
                                 cells_extra     = NULL,
                                 return = TRUE, file_name = NULL) {

  if (!dir.exists("Results")) dir.create("Results", recursive = TRUE)

  # ── extract top contributing cell groups per factor ───────────────────────
  # returns named list: factor → named numeric vector
  # names = cell group names (suffix stripped), values = NMF weights (all > 0)
  top_features_list <- extract_contributing_features(
    latent_factors,
    quantile_cutoff = quantile_cutoff
  )

  # ── binary composition matrix ─────────────────────────────────────────────
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

    # ── edge weights: factor → cell types ────────────────────────────────
    # for each cell type: sum of NMF weights of cell groups containing it
    # all weights are positive — no direction distinction

    # Interpretation of edge
    #   Magnitude = cumulative strength of CGs containing it
    #   Sign = direction of association with factor
    #   Large positive → enriched in positive CGs
    #   Large negative → enriched in negative CGs

    edge_weights <- colSums(
      comp_sub * top_features[rownames(comp_sub)]
    )
    edge_weights <- edge_weights[edge_weights > 0]

    # If a cell type appears mostly in positive CGs → positive edge
    # If appears mostly in negative CGs → negative edge
    # If balanced → near zero

    if (length(edge_weights) == 0) next

    # ── enrichment relative to background ────────────────────────────────
    # is this cell type more common in the top cell groups than globally?
    # enrichment >= enrich_thresh → cell type enriched across cell groups relative to baseline (“Is this cell type more common in my selected subset than expected based on its overall frequency?”)

    fg_freq <- colMeans(comp_sub)
    enrich  <- fg_freq / background_freq
    enrich[!is.finite(enrich)] <- NA

    # Interpretation
    # Suppose CD4.cells_immunosuppressive appears in 70% of the selected top features (fg_freq = 0.7)
    # - Its background frequency = 0.34408602
    # - Enrichment = 0.7 / 0.344 ≈ 2.03 → more than 2× enriched relative to baseline

    keep         <- names(which(enrich >= enrich_thresh))
    edge_weights <- edge_weights[names(edge_weights) %in% keep]

    if (length(edge_weights) == 0) next

    factor_celltype_weights[[factor_name]] <- edge_weights

    # ── network plot ──────────────────────────────────────────────────────
    edges <- data.frame(
      from = factor_name,
      to   = names(edge_weights),
      stringsAsFactors = FALSE
    )

    g <- igraph::graph_from_data_frame(edges, directed = FALSE)

    igraph::V(g)$type <- ifelse(igraph::V(g)$name == factor_name, "factor", "celltype")

    # node size — factor fixed, cell types scaled by edge weight
    ct_names  <- igraph::V(g)$name[igraph::V(g)$type == "celltype"]
    ct_sizes  <- scales::rescale(edge_weights[ct_names], to = c(6, 14))

    igraph::V(g)$size <- ifelse(
      igraph::V(g)$name == factor_name,
      18,
      ct_sizes[igraph::V(g)$name]
    )

    # node color — factor gold, cell types by edge weight magnitude
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
    # Positive if the cell type is mostly in positive-contributing cell groups.
    #
    # Negative if mostly in negative-contributing cell groups.
    #
    # Magnitude = cumulative contribution across all cell groups for this factor.

    lay <- igraph::layout_as_star(g,
                           center = which(igraph::V(g)$name == factor_name))

    if(return){
      pdf(paste0("Results/Network_", factor_name, "_", file_name, ".pdf"))
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

  # ── build TME states output ───────────────────────────────────────────────
  CG_states <- list()

  for (factor_name in names(factor_celltype_weights)) {

    edge_weights <- factor_celltype_weights[[factor_name]]

    if (length(edge_weights) == 0) next

    # single list of cell types per factor — no positive/negative split
    # ordered by weight descending
    CG_states[[factor_name]] <- sort(edge_weights, decreasing = TRUE)
  }

  return(CG_states)
}


#' Project cell group scores onto trained NMF latent factors
#'
#' Applies the positive-negative split to new samples and projects them
#' onto the latent space learned during training using the Moore-Penrose
#' pseudoinverse of the scaled basis matrix W.
#'
#' @param latent_spaces A list returned by \code{compute.latent_factors()},
#'   containing at minimum \code{W} (features x rank), \code{rank}, and
#'   the RcppML model object with \code{nmf_model$d} scaling vector.
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
#' @export
project_factors <- function(latent_spaces, scores_test) {

  scores_test <- as.matrix(scores_test)

  # ── Step 1: apply same pos-neg split as training ───────────────────────────
  test_pos <- pmax(scores_test, 0)
  test_neg <- pmax(-scores_test, 0)
  colnames(test_pos) <- paste0(colnames(scores_test), "_pos")
  colnames(test_neg) <- paste0(colnames(scores_test), "_neg")
  nmf_input_test <- cbind(test_pos, test_neg)   # samples x (2 * n_CGs)

  # ── Step 2: extract W and absorb d scaling ─────────────────────────────────
  # RcppML stores: A ≈ W %*% diag(d) %*% H
  # Must absorb d into W before pseudoinverse: W_scaled = W %*% diag(d)
  W <- latent_spaces$W                          # features x rank
  d <- latent_spaces$nmf_model$d                # length = rank

  W_scaled <- W %*% diag(d, nrow = length(d))  # features x rank, d absorbed

  # ── Step 3: restrict to features shared between training and test ──────────
  common_features <- intersect(rownames(W_scaled), colnames(nmf_input_test))

  if (length(common_features) == 0)
    stop("No overlapping features between training NMF input and test scores. ",
         "Check that column names of scores_test match training cell group names.")

  if (length(common_features) < nrow(W_scaled))
    warning(nrow(W_scaled) - length(common_features),
            " training features missing in test cohort — projection may be degraded.")

  W_sub <- W_scaled[common_features, , drop = FALSE]   # common x rank
  Y_sub <- t(nmf_input_test[, common_features,
                              drop = FALSE])            # common x samples

  # ── Step 4: project via pseudoinverse ─────────────────────────────────────
  # H_test = pinv(W_sub) %*% Y_sub  →  Z_test = t(H_test)
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
#' @param test_deconv A numeric matrix or data frame of deconvolution features
#'   for the test samples (samples x cell types). Column names must match those
#'   used during training.
#'
#' @return A numeric matrix (test samples x NMF factors) of non-negative
#'   projected factor scores.
#'
#' @export
project_test_factors <- function(train_processed, test_deconv) {

  cell_groups_scores <- compute.test.set(
    deconv_res        = train_processed$Processed_deconvolution,
    cell_groups       = train_processed$Cell_groups,
    features          = names(train_processed$Cell_groups[[2]]),
    deconvolution_test = test_deconv
  )

  project_factors(train_processed$Latent_spaces, cell_groups_scores)
}

#' Save a grid of scatter plots for significant module–feature pairs
#'
#' For each pair of columns from \code{matA} and \code{matB} whose p-value in
#' \code{p_mat} is below \code{pval}, draws a scatter plot with a regression
#' line and annotates it with the correlation coefficient. The full grid is
#' saved as an SVG file in \code{Results/}.
#'
#' @param matA A numeric matrix (samples x modules) for the y-axis of each panel.
#' @param matB A numeric matrix (samples x features) for the x-axis of each panel.
#' @param cor_mat A numeric matrix of correlation coefficients between columns of
#'   \code{matA} (rows) and \code{matB} (columns).
#' @param p_mat A numeric matrix of p-values corresponding to \code{cor_mat}.
#' @param file_name Character. Base name for the output SVG file (saved to
#'   \code{Results/<file_name>_scatter_grid.svg}).
#' @param pval Numeric. P-value cutoff for displaying a pair. Default 0.05.
#' @param cor_type Character. Label used in axis text (e.g., \code{"p"} for
#'   Pearson). Default \code{"p"}.
#' @param width Numeric. Width of the SVG output in inches. Default same as
#'   \code{height}.
#' @param height Numeric. Height of the SVG output in inches. Default same as
#'   \code{width}.
#' @param only_sig Logical. If \code{TRUE} (default), only pairs with
#'   \code{p <= pval} are plotted.
#' @param ncol Integer or \code{NULL}. Number of columns in the plot grid.
#'   If \code{NULL}, set to \code{ceiling(sqrt(n_pairs))}.
#'
#' @return Called for its side effect (saves SVG); returns \code{NULL} invisibly.
#'
#' @keywords internal
plot.module.scatter.grid <- function(matA, matB, cor_mat, p_mat,
                                     file_name,
                                     pval = 0.05,
                                     cor_type = "p",
                                     width = width,
                                     height = width,
                                     only_sig = TRUE,
                                     ncol = NULL) {   # allow NULL for auto

  pairs <- expand.grid(A = colnames(matA),
                       B = colnames(matB),
                       stringsAsFactors = FALSE)

  pairs$keep <- mapply(function(a, b) {
    if (only_sig) p_mat[a, b] <= pval else TRUE
  }, pairs$A, pairs$B)

  pairs <- pairs[pairs$keep, ]

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
      width = width, height = height)

  # Increase margins a bit if needed
  par(
    mfrow = c(nrow_grid, ncol),
    mar = c(5, 6, 5, 3),        # bottom, left, top, right
    cex.lab = 1.5,               # axis labels
    cex.axis = 1.3,              # axis tick labels
    cex.main = 1.8               # plot title
  )

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
         pch = 16, cex = 1.5,        # increase point size
         xlab = "",
         ylab = "",
         main = paste0(i, " vs ", j))

    # Add bigger axis titles
    title(
      xlab = j, ylab = i,
      cex.lab = 2        # increases both x and y labels
    )

    abline(lm(y ~ x, df), col = "steelblue", lwd = 2.5)  # thicker regression line

    mtext(
      paste0("R = ", signif(r, 3),
             " | p = ", format.pval(pv, digits = 2)),
      side = 3,
      line = 0,                  # move a bit higher
      cex = 1.3                     # increase text size
    )
  }

  dev.off()
}


#' Summarize TF module–trait associations as ANOVA boxplot grids
#'
#' For each categorical trait in \code{coldata}, runs a one-way ANOVA for each
#' TF module, performs Tukey HSD post-hoc tests for significant modules, and
#' saves a multi-panel SVG boxplot grid to \code{Results/}.
#'
#' @param tfs.modules A numeric matrix or data frame of TF module scores
#'   (samples x modules), typically from \code{compute.WTCNA()}.
#' @param coldata A data frame of sample metadata. Only character and factor
#'   columns are used as traits.
#' @param pval Numeric. ANOVA p-value threshold for significance. Default 0.05.
#' @param file.name Character. Base name appended to output SVG file names.
#' @param ncol Integer. Number of columns in the boxplot facet grid. Default 5.
#' @param y_min Numeric. Lower y-axis limit for boxplots. Default 0.
#' @param y_max Numeric. Upper y-axis limit for boxplots. Default 0.5.
#' @param width Numeric. Width of the SVG output in inches. Default 18.
#' @param height Numeric. Height of the SVG output in inches. Default 10.
#'
#' @return Called for its side effect (saves SVG files); returns \code{NULL}
#'   invisibly when no significant traits are found.
#'
#' @keywords internal
compute.metadata.association.boxplot_summary <- function(
    tfs.modules, coldata, pval = 0.05, file.name,
    ncol = 5, y_min = 0, y_max = 0.5, width = 18, height = 10
){

  coldata_cat <- coldata %>%
    dplyr::select(where(is.character) | where(is.factor))

  if (ncol(coldata_cat) == 0) {
    message("No categorical traits found.")
    return(invisible(NULL))
  }

  # ---- LOOP OVER EACH TRAIT ----
  for(tr in colnames(coldata_cat)){

    df_long <- cbind(tfs.modules, group = coldata_cat[[tr]]) %>%
      as.data.frame() %>%
      tidyr::pivot_longer(
        cols = colnames(tfs.modules),
        names_to = "module",
        values_to = "value"
      ) %>%
      dplyr::mutate(
        trait = tr,
        group = as.factor(group)
      )

    # ---- ANOVA ----
    anova_df <- df_long %>%
      dplyr::group_by(module) %>%
      dplyr::do({
        res <- rstatix::anova_test(data = ., dv = value, between = group)
        tibble(F = res$F, p = res$p)
      }) %>%
      dplyr::ungroup()

    sig_pairs <- anova_df %>% dplyr::filter(p < pval)
    if(nrow(sig_pairs) == 0){
      cat("No significant pairs found in trait:", tr)
      next
    }

    feature_labels <- sig_pairs %>%
      dplyr::mutate(
        feature_label = paste0(module, "\nF=", signif(F,3), ", p=", signif(p,3))
      ) %>%
      dplyr::distinct(module, feature_label) %>%
      tibble::deframe()

    df_long <- df_long %>%
      dplyr::semi_join(sig_pairs, by = "module")


    # ---- Tukey ----
    tukey_df <- df_long %>%
      dplyr::group_by(module) %>%
      dplyr::filter(n_distinct(group) > 1) %>%
      dplyr::group_modify(~ {
        tuk <- tukey_hsd(.x, value ~ group)
        if(nrow(tuk) == 0) return(NULL)
        add_xy_position(tuk, x = "group")
      }) %>%
      dplyr::ungroup()

    # ---- SAVE SVG PER TRAIT ----
    svg(paste0("Results/ANOVA_boxplot_summary_", file.name, "_", tr, ".svg"),
        width = width, height = height)

    p <- ggplot2::ggplot(df_long, aes(x = group, y = value, fill = group)) +
      geom_boxplot(width = 0.6, outlier.size = 0.4, alpha = 0.85) +
      geom_jitter(width = 0.2, size = 1, alpha = 0.7, color = "black") +
      facet_wrap(~ module, ncol = ncol,
                 labeller = labeller(module = feature_labels)) +
      coord_cartesian(ylim = c(y_min, y_max)) +
      labs(
        y = "Feature value",
        fill = "Group",
        title = paste0("Features-trait associations: ", tr),
        subtitle = paste0("One-way ANOVA with Tukey HSD | p-value < ", pval)
      ) +
      theme_bw(base_size = 12) +
      theme(
        strip.text = element_text(size = 14),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold"),
        plot.title = element_text(size = 20, face = "bold"),
        plot.subtitle = element_text(size = 16),
        legend.position = "top"
      )

    if(nrow(tukey_df) > 0){
      p <- p + stat_pvalue_manual(
        tukey_df,
        hide.ns = TRUE,
        size = 6,
        tip.length = 0.02,
        bracket.size = 0.6,
        inherit.aes = FALSE
      )
    }

    print(p)
    dev.off()
  }
}

#' Run multivariate feature-based GSEA using limma and Hallmark gene sets
#'
#' This function fits a multivariate linear model for each gene using all features
#' in \code{features_df} as continuous covariates. For each feature, it extracts
#' the moderated t-statistics and p-values, ranks genes, and performs GSEA using
#' the Hallmark gene sets from MSigDB. Optional dotplots for the top enriched
#' pathways can be saved as PDFs.
#'
#' @param RNA.tpm A numeric matrix or data frame of gene expression values
#'   (genes in rows, samples in columns).
#' @param features_df A data frame of continuous features (samples in rows,
#'   features in columns) to be modeled as covariates.
#' @param plot_dot Logical; if TRUE, generates and saves dotplots of top
#'   enriched Hallmark pathways for each feature. Default is TRUE.
#' @param top_n Integer; number of top pathways to display in the dotplot. Default is 10.
#' @param file_name Character; optional suffix for saved PDF files. Default is NULL.
#' @param width Numeric; width of the PDF plot in inches. Default is 8.
#' @param height Numeric; height of the PDF plot in inches. Default is 10.
#'
#' @return A list containing:
#' \describe{
#'   \item{DE_results}{A named list of \code{topTable} results for each feature,
#'     including logFC, moderated t-statistics, p-values, and adjusted p-values.}
#'   \item{GSEA_results}{A named list of \code{GSEA} results from
#'     \code{clusterProfiler} for each feature.}
#' }
#'
#' @details
#' The function works as follows:
#' \enumerate{
#'   \item Hallmark gene sets are retrieved from MSigDB using \code{msigdbr}.
#'   \item A multivariate linear model is fitted for each gene using \code{limma::lmFit}.
#'   \item Empirical Bayes moderation is applied via \code{limma::eBayes}.
#'   \item For each feature:
#'     \enumerate{
#'       \item Differential expression results are extracted using \code{topTable} for
#'         the coefficient of that feature.
#'       \item Genes are ranked by moderated t-statistics.
#'       \item Hallmark GSEA is performed using the ranked gene list.
#'       \item Optionally, a dotplot of the top enriched pathways is generated.
#'     }
#' }
#'
#' @import limma
#' @import msigdbr
#' @import clusterProfiler
#' @import dplyr
#' @import enrichplot
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
                             height = 10) {


  features_df = data.frame(features_df)
  # -----------------------------------------
  # Retrieve Hallmark gene sets
  # -----------------------------------------
  hallmark_df <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
  gene_sets <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)

  term2gene <- data.frame(
    term = rep(names(gene_sets), lengths(gene_sets)),
    gene = unlist(gene_sets)
  )

  # -----------------------------------------
  # Fit multivariate limma model (ALL features)
  # -----------------------------------------
  design <- stats::model.matrix(~ ., data = features_df)

  fit <- limma::lmFit(RNA.tpm, design)
  fit <- limma::eBayes(fit)

  deg_list <- list()
  gsea_results_list <- list()

  # -----------------------------------------
  # Loop over features (extract coefficients)
  # -----------------------------------------
  for (feature_name in colnames(features_df)) {

    cat("Processing:", feature_name, "\n")

    # Extract DE results for this feature
    deg <- limma::topTable(fit,
                    coef = feature_name,
                    number = Inf,
                    adjust.method = "BH",
                    sort.by = "P")

    deg_list[[feature_name]] <- deg

    # Rank genes by moderated t-statistic
    tstats <- fit$t[, feature_name]
    names(tstats) <- rownames(RNA.tpm)
    tstats <- sort(tstats, decreasing = TRUE)

    # Run Hallmark GSEA
    gsea_res <- clusterProfiler::GSEA(
      geneList = tstats,
      TERM2GENE = term2gene,
      verbose = FALSE,
      eps = 0,
      nPermSimple = 10000
    )

    gsea_results_list[[feature_name]] <- gsea_res


    # Optional plotting
    if (plot_dot && !is.null(gsea_res) && nrow(gsea_res@result) > 0) {

      if (!dir.exists("Results")) dir.create("Results")

      pdf(paste0("Results/GSEA_", feature_name, "_", file_name, ".pdf"),
          width = width,
          height = height)

      print(enrichplot::dotplot(gsea_res,
                              showCategory = top_n,
                              title = paste("Hallmark GSEA -", feature_name)))

      dev.off()
    }
  }

  return(list(
    DE_results = deg_list,
    GSEA_results = gsea_results_list
  ))
}


#' Run differential expression analysis with edgeR/limma-voom
#'
#' Filters low-expression genes, applies TMM normalization, runs voom
#' transformation, fits a linear model, and returns the top differentially
#' expressed genes via \code{limma::topTable}.
#'
#' @param counts A raw count matrix (genes x samples).
#' @param coldata A data frame of sample metadata whose row names match
#'   the column names of \code{counts}.
#' @param group_col Character. Name of the column in \code{coldata} used as
#'   the grouping factor for differential expression.
#' @param ref_level Character or \code{NULL}. Reference level for the group
#'   factor. If \code{NULL}, the default factor ordering is used.
#'
#' @return A data frame of differentially expressed genes (p.adj < 0.05) as
#'   returned by \code{limma::topTable}, with columns \code{logFC},
#'   \code{AveExpr}, \code{t}, \code{P.Value}, \code{adj.P.Val}, and \code{B}.
#'
#' @keywords internal
run_deg_analysis <- function(counts, coldata, group_col, ref_level = NULL) {
  # Prepare counts
  counts_mat <- as.matrix(counts)
  mode(counts_mat) <- "numeric"
  counts_mat <- counts_mat[, rownames(coldata)]

  # Create group factor
  group <- factor(coldata[[group_col]])

  # Set reference level if provided
  if (!is.null(ref_level)) {
    group <- stats::relevel(group, ref = ref_level)
  }

  # Create DGE object and filter
  dge <- edgeR::DGEList(counts = counts_mat, group = group)
  keep <- edgeR::filterByExpr(dge)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)

  # Design matrix
  design <- stats::model.matrix(~ group)

  # voom transformation
  v <- limma::voom(dge, design)

  # Fit model
  fit <- limma::lmFit(v, design)
  fit <- limma::eBayes(fit)
  # Extract coefficient name (second column of design)
  coef_name <- colnames(design)[2]

  # Get results
  res <- limma::topTable(fit, coef = coef_name, p.value = 0.05, number = Inf)

  return(res)
}

#' Derive TME meta-programs by clustering Hallmarks across NMF factors
#'
#' Hierarchically clusters Hallmark gene sets by their NES profile across NMF
#' factors to identify recurrent transcriptional programs in the TME.
#'
#' @param nes_mat A numeric matrix of NES values (Hallmarks x factors) as
#'   returned by \code{build_nes_matrix()}.
#' @param k Integer. Number of meta-programs (clusters) to extract. If
#'   \code{NULL} (default), estimated automatically from the dendrogram.
#' @param nes_thresh Numeric. Minimum absolute NES for a Hallmark to be
#'   considered active in a factor. Default 1.0.
#' @param plot Logical. If \code{TRUE} (default), saves a clustering heatmap.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{meta_programs}: Named list mapping meta-program labels to
#'     character vectors of Hallmark names.
#'   \item \code{hallmark_clusters}: Data frame with columns \code{Hallmark},
#'     \code{meta_program}, and \code{mean_NES}.
#'   \item \code{heatmap}: The \code{pheatmap} object.
#'   \item \code{k}: The number of meta-programs used.
#' }
#'
#' @export
derive_meta_programs <- function(nes_mat,
                                  k          = NULL,
                                  nes_thresh = 1.0, file_name = NULL,
                                  plot       = TRUE) {

  # ── filter: keep Hallmarks active in at least one factor ─────────────────
  active <- apply(nes_mat, 1, function(x) any(abs(x) >= nes_thresh))
  nes_filtered <- nes_mat[active, , drop = FALSE]

  if (nrow(nes_filtered) < 3) {
    stop("Too few active Hallmarks — lower nes_thresh")
  }

  # ── cluster Hallmarks by NES profile across factors ───────────────────────
  # distance between Hallmarks = dissimilarity of their NES profiles
  d  <- dist(nes_filtered, method = "euclidean")
  hc <- hclust(d, method = "ward.D2")

  # ── determine k ───────────────────────────────────────────────────────────
  if (is.null(k)) {
    # use elbow on within-cluster sum of squares
    max_k <- min(10, nrow(nes_filtered) - 1)
    wss   <- sapply(2:max_k, function(ki) {
      cuts <- cutree(hc, k = ki)
      sum(sapply(unique(cuts), function(cl) {
        members <- nes_filtered[cuts == cl, , drop = FALSE]
        if (nrow(members) == 1) return(0)
        sum(apply(members, 2, var)) * (nrow(members) - 1)
      }))
    })
    # second derivative elbow
    d2 <- diff(diff(wss))
    k  <- which.min(d2) + 2L
    message("Selected k = ", k, " meta-programs")
  }

  cluster_assignments <- cutree(hc, k = k)

  # ── build meta-program list ───────────────────────────────────────────────
  meta_programs <- lapply(seq_len(k), function(cl) {
    names(cluster_assignments)[cluster_assignments == cl]
  })

  # ── name each meta-program by its dominant Hallmarks ─────────────────────
  meta_names <- paste0("MP", seq_len(k))
  names(meta_programs) <- meta_names

  # ── summary data frame ────────────────────────────────────────────────────
  hallmark_df <- data.frame(
    hallmark      = names(cluster_assignments),
    meta_program  = meta_names[cluster_assignments],
    cluster_id    = cluster_assignments,
    mean_NES      = rowMeans(nes_filtered)[names(cluster_assignments)],
    stringsAsFactors = FALSE
  ) %>%
    arrange(cluster_id, desc(abs(mean_NES)))

  # ── heatmap ───────────────────────────────────────────────────────────────
  p <- NULL
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

    pdf(paste0("Results/TCGA_meta_programs_", file_name, ".pdf"),
        width = 8, height = 10)
    p <- pheatmap::pheatmap(
      nes_filtered,
      cluster_rows     = hc,
      cluster_cols     = TRUE,
      cutree_rows      = k,
      annotation_row   = ann_row,
      annotation_colors = list(Meta_program = mp_colors),
      color            = colorRampPalette(
                           c("#2471A3", "white", "#C0392B"))(100),
      breaks           = seq(-3, 3, length.out = 101),
      fontsize_row     = 8,
      fontsize_col     = 9,
      main             = paste0("Hallmark meta-programs across NMF factors\n",
                                "(TCGA reference, k = ", k, ")")
    )
    dev.off()

  }

  list(
    meta_programs     = meta_programs,
    hallmark_clusters = hallmark_df,
    nes_matrix        = nes_filtered,
    heatmap           = p,
    k                 = k
  )
}

#' Build a Hallmarks x factors NES matrix from GSEA results
#'
#' Combines the per-factor GSEA outputs from \code{compute_factor_gsea()} into
#' a single matrix. Hallmarks not significant in a given factor are filled with 0.
#'
#' @param gsea_results Output list from \code{compute_factor_gsea()}, containing
#'   a \code{GSEA_results} element (named list of \code{enrichResult} objects,
#'   one per factor).
#'
#' @return A numeric matrix of NES values with Hallmarks as rows and NMF factors
#'   as columns. Missing Hallmark–factor combinations are set to 0.
#'
#' @export
build_nes_matrix <- function(gsea_results) {

  gsea_list <- gsea_results$GSEA_results

  # get all Hallmarks seen across any factor
  all_hallmarks <- unique(unlist(
    lapply(gsea_list, function(res) {
      if (is.null(res) || nrow(res@result) == 0) return(character(0))
      res@result$ID
    })
  ))

  # build NES matrix: hallmarks x factors
  nes_mat <- sapply(names(gsea_list), function(fac) {

    res <- gsea_list[[fac]]

    # empty factor
    if (is.null(res) || nrow(res@result) == 0) {
      return(setNames(rep(0, length(all_hallmarks)), all_hallmarks))
    }

    # NES per hallmark — fill missing with 0
    nes <- setNames(res@result$NES, res@result$ID)
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
#' by computing the mean NES of the meta-program's Hallmarks in that factor.
#' The meta-program with the highest mean NES is the best match.
#'
#' @param gsea_study Output from compute_factor_gsea() on study cohort.
#' @param cancer_type Character. TCGA cancer type abbreviation (e.g., \code{"blca"}, \code{"brca"},
#'   \code{"cesc"}, \code{"chol"}, \code{"coad"}, \code{"skcm"}) identifying which pre-built
#'   TCGA meta-program file to load from \code{inst/extdata/}.
#' @param nes_thresh Minimum absolute NES to report a Hallmark as active.
#'                   Default 1.0.
#'
#' @return Data frame with one row per study factor:
#'   factor, best_MP, best_score, all_scores, active_hallmarks
#'
#' @export
map_to_metaprograms <- function(gsea_study,
                                 cancer_type,
                                 nes_thresh = 1.0) {

  cancer_type <- tolower(cancer_type)
  mp_file <- system.file("extdata",
                         paste0("TCGA_meta_programs_", cancer_type, ".RData"),
                         package = "CellTFusion")
  if (mp_file == "") {
    stop("No meta-program file found for cancer type '", cancer_type,
         "'. Available types: blca, brca, cesc, chol, coad, skcm.")
  }
  load(mp_file)

  # ── build study NES matrix: all 50 Hallmarks x study factors ──────────────
  all_hallmarks <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H") %>%
    dplyr::pull(gs_name) %>% unique() %>% sort()

  nes_study <- sapply(names(gsea_study$GSEA_results), function(fac) {
    scores <- setNames(rep(0, length(all_hallmarks)), all_hallmarks)
    res    <- gsea_study$GSEA_results[[fac]]
    if (!is.null(res) && nrow(res@result) > 0) {
      nes <- setNames(res@result$NES, res@result$ID)
      scores[names(nes)] <- nes
    }
    scores
  })
  # nes_study: 50 hallmarks x study factors

  # ── score each study factor against each meta-program ─────────────────────
  # score = mean NES of the meta-program Hallmarks in this study factor
  results <- lapply(colnames(nes_study), function(fac) {

    fac_nes <- nes_study[, fac]
    
    # score per meta-program
    mp_scores <- sapply(names(meta_programs$meta_programs), function(mp) {
      h <- intersect(meta_programs$meta_programs[[mp]], names(fac_nes))
      if (length(h) == 0) return(NA_real_)
      mean(fac_nes[h])
    })

    best_mp    <- names(which.max(mp_scores))
    best_score <- mp_scores[best_mp]

    # active Hallmarks for reporting
    active_h <- names(fac_nes)[abs(fac_nes) >= nes_thresh]
    active_h <- stringr::str_remove(active_h, "^HALLMARK_")

    data.frame(
      factor           = fac,
      best_MP          = best_mp,
      best_score       = round(best_score, 3),
      all_scores       = paste(
        sprintf("%s:%.3f", names(mp_scores), mp_scores),
        collapse = ", "
      ),
      active_hallmarks = paste(active_h, collapse = ", "),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, results)
}
