#' Computes EpiDISH
#'
#' @param TPM_matrix A matrix with TPM normalized counts (genes symbols as rows and samples as columns).
#' @param signature_file The signature file to use.
#' @param name_signature Signature name to set for the deconvolution results.
#'
#' @return A matrix with cell abundance deconvolve with EpiDISH
#' @export
#'
#' @examples
#'
#' epidish <- computeEpiDISH(TPM_matrix, signature, signature_name)
#'
computeEpiDISH = function(TPM_matrix, signature_file, name_signature){
  require(EpiDISH)
  epi <- epidish(TPM_matrix, as.matrix(signature_file), method = "RPC", maxit = 200)
  epidish = epi$estF

  colnames(epidish) = paste0("Epidish_", name_signature, "_", colnames(epidish))
  colnames(epidish) <- colnames(epidish) %>%
    str_replace_all(., " ", "_")

  return(epidish)
}
