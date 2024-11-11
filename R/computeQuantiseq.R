computeQuantiseq <- function(TPM_matrix) {
  require(immunedeconv)
  TPM_matrix = TPM_matrix[rownames(TPM_matrix)%in%rownames(immunedeconv::dataset_racle$expr_mat),] #To avoid problems regarding gene names (quantiseq error)

  quantiseq = immunedeconv::deconvolute(TPM_matrix, "quantiseq", tumor = T) %>%
    column_to_rownames("cell_type") %>%
    t()

  colnames(quantiseq) = paste0("Quantiseq_", colnames(quantiseq))
  colnames(quantiseq) <- colnames(quantiseq) %>%
    str_replace_all(., " ", "_")

  return(quantiseq)
}
