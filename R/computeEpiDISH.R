computeEpiDISH = function(TPM_matrix, signature_file, name_signature){
  require(EpiDISH)
  epi <- epidish(TPM_matrix, as.matrix(signature_file), method = "RPC", maxit = 200)
  epidish = epi$estF

  colnames(epidish) = paste0("Epidish_", name_signature, "_", colnames(epidish))
  colnames(epidish) <- colnames(epidish) %>%
    str_replace_all(., " ", "_")

  return(epidish)
}
