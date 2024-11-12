#' Compute CibersortX (CBSX) in parallel across multiple signatures
#'
#' @param TPM_matrix A matrix with TPM normalized counts (genes symbols as rows and samples as columns).
#' @param signatures Path where signatures files are located
#' @param name Credential email for running CibersortX.
#' @param password Credential token for running CibersortX.
#' @param workers Number of processes available to run on parallel.
#'
#' @return A matrix with cell abundance deconvolve with CBSX
#' @export
#'
#' @examples
#'
#' cbsx <- computeCBSX_parallel(TPM_matrix, db, cbsx.name, cbsx.token, workers)
#'
computeCBSX_parallel = function(TPM_matrix, signatures, name, password, workers){
  cl = parallel::makeCluster(workers)
  registerDoParallel(cl)

  cbsx = foreach (i=1:length(signatures), .combine=cbind) %dopar% {
    source("src/environment_set.R")
    signature <- read.delim(signatures[[i]], row.names=1)
    signature_name = str_split(basename(signatures[[i]]), "\\.")[[1]][1]
    computeCBSX(TPM_matrix, signature, name, password, signature_name)
  }

  parallel::stopCluster(cl)
  unregister_dopar()

  return(cbsx)
}
