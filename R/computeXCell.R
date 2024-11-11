computeXCell <- function(TPM_matrix) {
  require(immunedeconv)

  xcell = immunedeconv::deconvolute(TPM_matrix, "xcell") %>%
    column_to_rownames("cell_type") %>%
    t()

  colnames(xcell) = paste0("XCell_", colnames(xcell))
  colnames(xcell) <- colnames(xcell) %>%
    str_replace_all(., " ", "_")

  return(xcell)
}
