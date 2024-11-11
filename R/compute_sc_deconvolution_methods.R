#' Compute second-generation deconvolution methods
#'
#' @param raw_counts A matrix with raw counts (samples as columns and genes symbols as rows)
#' @param sc_object A matrix with the counts from scRNAseq object (genes as rows and cells as columns)
#' @param cell_annotations A character vector with the cell labels (need to be of the same order as in the sc_object)
#' @param samples_ids A character vector with the samples labels (need to be of the same order as in the sc_object)
#' @param name_object Signature name to use in the generated single cell signature for deconvolving the bulk RNAseq data
#' @param normalized Boolean value to specify if raw_counts need to be normalized (If no raw_counts are available and argument corresponds to already normalized counts this arguments needs to be set to False)
#' @param n_cores Number of cores to use for paralellization. If no number is set, detectCores() - 1 will be set as the number.
#' @param cbsx_name Cibersortx credential mail if CBSX will be run
#' @param cbsx_token Cibersortx credential token if CBSX will be run
#'
#' @return A matrix of deconvolution features across samples from your bulk counts based on the second generation methods
#' @export
#'
#' @examples
#'
#' deconv_sc = compute_sc_deconvolution_methods(raw.counts, sc_matrix, cell_annotations, cell_samples, name_sc_signature, normalized = T, n_cores = 4, cbsx_name = "XXX", cbsx_token = "XXX")
#'
compute_sc_deconvolution_methods = function(raw_counts, sc_object, cell_annotations, samples_ids, name_object, normalized = T,
                                            n_cores = NULL, cbsx_name = NULL, cbsx_token = NULL){

  if(normalized){
    bulk_counts = ADImpute::NormalizeTPM(raw_counts)
  }else{
    bulk_counts = raw_counts
  }


  if(is.null(n_cores)){
    n_cores = parallel::detectCores() - 1
    message("\nUsing ",n_cores, " cores available for running.........................................................\n")
  }

  message("\nRunning AutogeneS...............................................................\n")
  autogenes = deconvolute_autogenes(bulk_gene_expression = bulk_counts, single_cell_object = sc_object,
                                    cell_type_annotations = as.character(cell_annotations), verbose = T)$proportions

  message("\nRunning BayesPrism...............................................................\n")
  bayesprism = deconvolute_bayesprism(bulk_gene_expression = raw_counts, single_cell_object = sc_object,
                                      cell_type_annotations = cell_annotations, n_cores = n_cores)$theta

  message("\nRunning Bisque...............................................................\n")
  bisque = deconvolute_bisque(bulk_gene_expression = as.matrix(raw_counts), single_cell_object = sc_object,
                              cell_type_annotations = cell_annotations, batch_ids = samples_ids, verbose = T)$bulk_props


  message("\nRunning CibersortX...............................................................\n")
  if(is.null(cbsx_name)&is.null(cbsx_token)){
    message("CibersortX credentials not given, method will not be used.\n")
    cibersortx = NULL
  }else{
    set_cibersortx_credentials(cbsx_name, cbsx_token)
    model_cbsx <- omnideconv::build_model(sc_object, as.character(cell_annotations),
                                          batch_ids = samples_ids, method = "cibersortx")
    cibersortx = deconvolute_cibersortx(bulk_gene_expression = as.matrix(bulk_counts),
                                        signature = model_cbsx, verbose = T)
  }

  message("\nRunning CPM...............................................................\n")
  cpm = deconvolute_cpm(bulk_gene_expression = data.frame(raw_counts), single_cell_object = sc_object, no_cores = 4, #raw counts
                        cell_type_annotations = as.character(cell_annotations), verbose = T)$cellTypePredictions

  message("\nRunning DWLS...............................................................\n")
  model_dwls <- omnideconv::build_model_dwls(sc_object, as.character(cell_annotations),
                                             batch_ids = samples_ids, dwls_method = "mast_optimized", ncores = n_cores)
  deconvolution_dwls = deconvolute_dwls(bulk_gene_expression = bulk_counts, signature = model_dwls,
                                        dwls_submethod = "DampenedWLS", verbose = T)

  message("\nRunning MOMF...............................................................\n")
  model_momf <- omnideconv::build_model(bulk_gene_expression = raw_counts, sc_object,
                                        as.character(cell_annotations), batch_ids = samples_ids, method = "momf") #raw counts
  momf = deconvolute_momf(bulk_gene_expression = as.matrix(bulk_counts), single_cell_object = sc_object,
                          signature = model_momf, method = "KL", verbose = T)

  message("\nRunning MuSiC...............................................................\n")
  music = deconvolute_music(bulk_gene_expression = as.matrix(bulk_counts), single_cell_object = sc_object,
                            cell_type_annotations = cell_annotations,  batch_ids = samples_ids, verbose = T)$Est.prop.weighted

  message("\nRunning SCDC...............................................................\n")
  scdc = deconvolute_scdc(bulk_gene_expression = as.matrix(bulk_counts), single_cell_object = sc_object,
                          cell_type_annotations = cell_annotations, batch_ids = samples_ids, verbose = T)$prop.est.mvw


  results = list(AutogeneS = autogenes, BayesPrism = bayesprism, CBSX = cibersortx, CPM = cpm,
                 DWLS = deconvolution_dwls, MOMF = momf, MuSic = music, SCDC = scdc)


  results <- lapply(names(results), function(method) {
    deconv_method <- results[[method]]

    colnames(deconv_method) <- paste0(method, "_", name_object, "_", colnames(deconv_method))
    colnames(deconv_method) <- str_replace_all(colnames(deconv_method), " ", "_")

    return(data_element)
  })

  results = do.call(cbind, results)

  return(results)

  ### METHODS NOT YET IMPLEMENTED

  # 1. BSeq-sc: Need CIBERSORT source code
  # message("\nRunning BSeq-sc...............................................................\n")
  # bseqsc = deconvolute_bseqsc(bulk_gene_expression = as.matrix(bulk.data), single_cell_object = sc_object,
  #                             cell_type_annotations = cell_annotations, batch_ids = samples_ids, verbose = T)

  # 2. CDSeq: Crash machine
  #message("\nRunning CDSeq...............................................................\n")
  # cdseq = deconvolute(bulk_gene_expression = raw_counts, single_cell_object = counts.matrix, no_cores = 6, method = "cdseq",
  #                     cell_type_annotations = cell_annotations, batch_ids = samples_ids, verbose = T)

  # 3. SCADEN: Takes a lot of time
  # message("\nRunning Scaden...............................................................\n")
  # model_scaden <- omnideconv::build_model(bulk_gene_expression = bulk.data, counts.matrix, as.character(cell_annotations),
  #                                         batch_ids = samples_ids, method = "scaden")
  # scaden = deconvolute_scaden(bulk_gene_expression = bulk_counts,
  #                             signature = model_scaden, verbose = T)

}
