computeCBSX = function(TPM_matrix, signature_file, name, password, name_signature){
  set_cibersortx_credentials(name, password)
  cbsx = omnideconv::deconvolute_cibersortx(TPM_matrix, signature_file)

  colnames(cbsx) = paste0("CBSX_", name_signature, "_", colnames(cbsx))
  colnames(cbsx) <- colnames(cbsx) %>%
    str_replace_all(., " ", "_")

  return(cbsx)
}
