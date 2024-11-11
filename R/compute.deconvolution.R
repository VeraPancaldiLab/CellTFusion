compute.deconvolution <- function(raw.counts, methods = c("Quantiseq", "MCP", "xCell", "CBSX", "Epidish", "DeconRNASeq", "DWLS"), signatures_exclude = NULL, normalized = T, doParallel = F, workers = NULL, return = T,
                                  credentials.mail = NULL, credentials.token = NULL, sc_deconv = F, ncores = NULL, sc_matrix = NULL, cell_annotations = NULL, cell_samples = NULL, name_sc_signature = NULL, file_name = NULL){

  path_signatures = 'src/signatures'

  if(normalized == T){
    cat("Performing TPM normalization log transformed...............................................................\n\n")
    TPM_matrix = data.frame(ADImpute::NormalizeTPM(raw.counts))
  }else{ #If no raw counts are available
    TPM_matrix = data.frame(raw.counts)
  }

  cat("Running deconvolution using the following methods...............................................................\n\n")
  for (method in methods) {
    cat("* ", method, "\n", sep = "")
  }

  if("Quantiseq" %in% methods){
    cat("\nRunning Quantiseq...............................................................\n")
    quantiseq = computeQuantiseq(TPM_matrix)}
  if("MCP" %in% methods){
    cat("\nRunning MCPCounter...............................................................\n")
    mcp = computeMCP(TPM_matrix, path_signatures)}
  if("xCell" %in% methods){
    xcell = computeXCell(TPM_matrix)
    cat("\nRunning XCell...............................................................\n")}

  default_sig = c("Quantiseq", "MCP", "xCell")
  methods = methods[!(methods %in% default_sig)]
  if(length(methods) == 0){
    methods = NULL
  }
  deconv_sig = compute_methods_variable_signature(TPM_matrix, signatures = path_signatures, method = methods, exclude = signatures_exclude, cbsx.name = credentials.mail, cbsx.token = credentials.token, doParallel, workers)

  deconv_default <- NULL
  if (exists("quantiseq")) {
    deconv_default <- quantiseq
  }
  if (exists("mcp")) {
    if (is.null(deconv_default)) {
      deconv_default <- mcp
    } else {
      deconv_default <- cbind(deconv_default, mcp)
    }
  }
  if (exists("xcell")) {
    if (is.null(deconv_default)) {
      deconv_default <- xcell
    } else {
      deconv_default <- cbind(deconv_default, xcell)
    }
  }

  if(is.null(deconv_sig)){
    all_deconvolution_table = deconv_default
  }else{
    all_deconvolution_table = cbind(deconv_default, deconv_sig)
  }

  if(sc_deconv){
    message("Running second generation cell-type deconvolution methods using scRNAseq\n")
    if(is.null(sc_object)==T){
      stop("No single cell object has been provided for deconvolution.")
    }else{
      deconv_sc = compute_sc_deconvolution_methods(raw.counts, sc_matrix, cell_annotations, cell_samples, name_sc_signature, normalized = normalized,
                                                   n_cores = ncores, cbsx_name = credentials.mail, cbsx_token = credentials.token)
      all_deconvolution_table = cbind(data.frame(all_deconvolution_table), deconv_sc)
    }
  }

  deconvolution = compute.deconvolution.preprocessing(data.frame(all_deconvolution_table))

  if(return){
    write.csv(deconvolution, paste0("Results/Deconvolution_", file_name, ".csv"))
  }

  return(deconvolution)

}
