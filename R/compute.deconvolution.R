#' Compute deconvolution
#'
#'The function calculates cell abundance based on cell type signatures using different methods and signatures. Methods available are Quantiseq, MCP, XCell, CibersortX, EpiDISH, DWLS and DeconRNASeq. Provided signatures included signatures based on bulk and methylation data (7 methods and 10 signature in total). Signatures are present in the src/signatures directory, user can add its own signatures by adding the .txt files in this same folder. Second generation methods to perform deconvolution based on single cell data are also available if scRNAseq object is provided.
#'
#' @param raw.counts A matrix with the raw counts (samples as columns and genes symbols as rows)
#' @param methods A character vector with the deconvolution methods to run. Default are "Quantiseq", "MCP", "xCell", "CBSX", "Epidish", "DeconRNASeq", "DWLS"
#' @param signatures_exclude A character vector with the signatures to exclude from the src/signatures folder.
#' @param normalized If raw.counts are not available, user can input its normalized counts. In that case this argument need to be set to False.
#' @param doParallel Whether to do or not parallelization. Only CBSX and DWLS methods will run in parallel.
#' @param workers Number of processes available to run on parallel. If no number is set, this will correspond to detectCores() - 1
#' @param return Whether to save or not the csv file with the deconvolution features
#' @param credentials.mail (Optional) Credential email for running CibersortX. If not provided, cibersortX method will not be run.
#' @param credentials.token (Optional) Credential token for running CibersortX. If not provided, cibersortX method will not be run.
#' @param sc_deconv Whether to run or not deconvolution methods based on single cell.
#' @param ncores If sc_deconv = T, number of cores to use for running the second-generation methods.
#' @param sc_matrix If sc_deconv = T, the matrix of counts across cells from the scRNAseq object is provided.
#' @param cell_annotations If sc_deconv = T, a character vector indicating the cell labels (same order as the count matrix)
#' @param cell_samples If sc_deconv = T, a character vector indicating the cell samples IDs (same order as the count matrix)
#' @param name_sc_signature If sc_deconv = T, the name you want to give to the signature generated
#' @param file_name File name for the csv files and plots saved in the Results/ directory
#'
#' @return
#'
#' A matrix of cell type deconvolution features across samples
#'
#' @export
#'
#' @examples
#'
#' deconv = compute.deconvolution(raw.counts, normalized = T, credentials.mail = "xxxx", credentials.token = "xxxxxx", file_name = "Tutorial")
#' deconv = compute.deconvolution(raw.counts, normalized = T, credentials.mail = "xxxx", credentials.token = "xxxxxx", methods = c("Quantiseq", "MCP", "XCell", "DWLS"), file_name = "Test")
#' deconv = compute.deconvolution(raw.counts, normalized = T, credentials.mail = "xxxx", credentials.token = "xxxxxx", signatures_exclude = "BPRNACan", file_name = "Tutorial")
#' deconv = compute.deconvolution(raw.counts, normalized = T, credentials.mail = "xxxx", credentials.token = "xxxxxx", sc_deconv = T, sc_matrix = sc.object, cell_annotations = cell_labels, cell_samples = bath_ids, name_sc_signature = "Signature_test", file_name = "Test")
#'
#'
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
