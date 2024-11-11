#' Compute deconvolution methods with variable signatures
#'
#' @param TPM_matrix A matrix with TPM normalized counts (samples as columns and genes symbols as rows)
#' @param signatures A path with a directory where signatures are located
#' @param methods A character vector with the methods to compute (Default methods are CBSX, Epidish, DeconRNASeq and DWLS)
#' @param exclude (Optional) A character vector with the signature to exclude
#' @param cbsx.name Cibersortx credential mail if CBSX will be run
#' @param cbsx.token Cibersortx credential token if CBSX will be run
#' @param doParallel Boolean value to specify if DWLS and CBSX should run in parallel (default is False)
#' @param workers Number of worker process to run during parallelization (default is NULL)
#'
#' @return A matrix with the deconvolution features corresponding to all combinations of methods-signatures specified
#' @export
#'
#' @examples
#'
#' deconvolution = compute_methods_variable_signature(TPM_matrix, "~/src/signatures", methods = c("CBSX", "Epidish", "DeconRNASeq", "DWLS"), exclude = "BPRNACan", cbsx.name = "XXX", cbsx.token = "XXX", doParallel = T, workers = 4)
#'
compute_methods_variable_signature = function(TPM_matrix, signatures, methods = c("CBSX", "Epidish", "DeconRNASeq", "DWLS"), exclude = NULL, cbsx.name, cbsx.token, doParallel = F, workers = NULL){

  db=list.files(signatures, full.names = T, pattern = "\\.txt$")
  name_exclude = c()
  if(is.null(methods)==F){
    cat("\nThe following method-signature combinations are going to be calculated...............................................................\n")

    cat("\nMethods\n")
    for (method in methods) {
      cat("* ", method, "\n", sep = "")
    }
    cat("\nSignatures\n")
    for (i in 1:length(db)) {
      name = str_split(basename(db[[i]]), "\\.")[[1]][1]
      cat("* ", name, "\n", sep = "")
      if(is.null(exclude)==F && name %in% exclude){
        name_exclude = c(name_exclude, name)
      }
    }

    if(length(name_exclude)>0){
      cat("\nExcluding signatures: ", paste0(name_exclude, collapse = ", "), "\n")
    }

    deconvolution = list()

    if("CBSX" %in% methods){
      if(is.null(cbsx.name)==T || is.null(cbsx.token)==T){
        cat("\nYou select to run CBSX but no credentials were found")
        cat("\nPlease set your credentials in the function for running CibersortX")
        stop()
      }
    }

    for (i in 1:length(db)) {
      signature <- read.delim(db[[i]], row.names=1)
      signature_name = str_split(basename(db[[i]]), "\\.")[[1]][1]
      if(!is.null(exclude) && signature_name %in% exclude) {
        next
      }else{
        if("DeconRNASeq"%in%methods){
          cat("\nRunning DeconRNASeq...............................................................\n\n")
          deconrnaseq <- computeDeconRNASeq(TPM_matrix, signature, signature_name)}
        if("Epidish"%in%methods){
          cat("\nRunning Epidish...............................................................\n\n")
          epidish <- computeEpiDISH(TPM_matrix, signature, signature_name)}
        if("DWLS"%in%methods){
          if(doParallel == F){
            cat("\nRunning DWLS...............................................................\n\n")
            dwls <- computeDWLS(TPM_matrix, signature, signature_name)}}
        if("CBSX"%in%methods){
          if(doParallel == F){
            cat("\nRunning CBSX...............................................................\n\n")
            cbsx <- computeCBSX(TPM_matrix, signature, cbsx.name, cbsx.token, signature_name)}}
        combined_data <- NULL
        if (exists("deconrnaseq")) {
          combined_data <- deconrnaseq
        }
        if (exists("epidish")) {
          if (is.null(combined_data)) {
            combined_data <- epidish
          } else {
            combined_data <- cbind(combined_data, epidish)
          }
        }
        if (exists("cbsx")) {
          if (is.null(combined_data)) {
            combined_data <- cbsx
          } else {
            combined_data <- cbind(combined_data, cbsx)
          }
        }
        if (exists("dwls")) {
          if (is.null(combined_data)) {
            combined_data <- dwls
          } else {
            combined_data <- cbind(combined_data, dwls)
          }
        }

        deconvolution[[i]] <- combined_data
      }
    }

    deconv = do.call(cbind, deconvolution)

    if("DWLS"%in%methods & doParallel == T){
      cat("\nRunning DWLS in parallel using", workers,"workers...............................................................\n\n")
      dwls <- computeDWLS_parallel(TPM_matrix, db, workers)
      deconv = cbind(deconv, dwls)
    }

    if("CBSX"%in%methods & doParallel == T){
      cat("\nRunning CBSX in parallel using", workers,"workers...............................................................\n\n")
      cbsx <- computeCBSX_parallel(TPM_matrix, db, cbsx.name, cbsx.token, workers)
      deconv = cbind(deconv, cbsx)
    }

    return(deconv)
  }else{
    cat("\nNo methods to be calculated using variable signatures.")
    return(NULL)
  }

}
