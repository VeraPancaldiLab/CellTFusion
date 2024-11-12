#' Computes MCPcounter
#'
#' @param TPM_matrix A matrix with TPM normalized counts (genes symbols as rows and samples as columns).
#' @param genes_path Path containing the MCP genes
#'
#' @return A matrix with cell enrichment scores from MCP
#' @export
#'
#' @examples
#'
#' mcp = computeMCP(TPM_matrix, path_signatures)
#'
computeMCP <- function(TPM_matrix, genes_path) {
  require(MCPcounter)
  genes <- read.table(paste0(genes_path, "/MCPcounter/MCPcounter-genes.txt"), sep = "\t", stringsAsFactors = FALSE, header = TRUE, colClasses = "character", check.names = FALSE)
  mcp <- MCPcounter.estimate(TPM_matrix, genes = genes, featuresType = "HUGO_symbols", probesets = NULL) %>%
    t()

  colnames(mcp) = paste0("MCP_", colnames(mcp))
  colnames(mcp) <- colnames(mcp) %>%
    str_replace_all(., " ", "_")

  return(mcp)
}
