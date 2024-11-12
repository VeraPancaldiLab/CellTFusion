#' Computes DeconRNASeq
#'
#' @param TPM_matrix A matrix with TPM normalized counts (genes symbols as rows and samples as columns).
#' @param signature_file The signature file to use.
#' @param name_signature Signature name to set for the deconvolution results.
#'
#' @return A matrix with cell abundance deconvolve with DeconRNASeq
#' @export
#'
#' @examples
#'
#' deconrnaseq <- computeDeconRNASeq(TPM_matrix, signature, signature_name)
#'
computeDeconRNASeq = function(TPM_matrix, signature_file, name_signature){
  require(DeconRNASeq)
  decon <- DeconRNASeq(TPM_matrix, data.frame(signature_file))
  deconRNAseq = decon$out.all
  rownames(deconRNAseq) = colnames(TPM_matrix)

  colnames(deconRNAseq) = paste0("DeconRNASeq_", name_signature, "_", colnames(deconRNAseq))
  colnames(deconRNAseq) <- colnames(deconRNAseq) %>%
    str_replace_all(., " ", "_")

  return(deconRNAseq)
}
