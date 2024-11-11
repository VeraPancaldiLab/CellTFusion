computeDWLS = function(TPM_matrix, signature_file, name_signature){
  require(DWLS)
  genes = rownames(signature_file)

  signature_file <- signature_file %>%
    apply(., 2, as.numeric) %>%
    data.frame() %>%
    mutate("Genes" = genes) %>%
    column_to_rownames("Genes") %>%
    as.matrix()

  dwls = omnideconv::deconvolute_dwls(TPM_matrix, signature_file, dwls_submethod = "SVR", verbose = T)

  colnames(dwls) = paste0("DWLS_", name_signature, "_", colnames(dwls))
  colnames(dwls) <- colnames(dwls) %>%
    str_replace_all(., " ", "_")

  return(dwls)
}
