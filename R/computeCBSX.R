#' Computes CibersortX (CBSX) using one signature
#'
#' @param TPM_matrix A matrix with TPM normalized counts (genes symbols as rows and samples as columns).
#' @param signature_file The signature file to use.
#' @param name Credential email for running CibersortX.
#' @param password Credential token for running CibersortX.
#' @param name_signature Signature name to set for the deconvolution results.
#'
#' @return A matrix with cell abundance deconvolve with CBSX
#' @export
#'
#' @examples
#'
#' cbsx <- computeCBSX(TPM_matrix, signature, cbsx.name, cbsx.token, signature_name)
#'
computeCBSX = function(TPM_matrix, signature_file, name, password, name_signature){
  set_cibersortx_credentials(name, password)
  cbsx = omnideconv::deconvolute_cibersortx(TPM_matrix, signature_file)

  colnames(cbsx) = paste0("CBSX_", name_signature, "_", colnames(cbsx))
  colnames(cbsx) <- colnames(cbsx) %>%
    str_replace_all(., " ", "_")

  return(cbsx)
}
