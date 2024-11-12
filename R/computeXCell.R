#' Computes XCell
#'
#' @param TPM_matrix A matrix with TPM normalized counts (genes symbols as rows and samples as columns).
#'
#' @return A matrix with cell enrichment scores from XCell.
#' @export
#'
#' @examples
#'
#' xcell = computeXCell(TPM_matrix)
#'
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
