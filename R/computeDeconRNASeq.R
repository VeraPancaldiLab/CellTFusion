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
