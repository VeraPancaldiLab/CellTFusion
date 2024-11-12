#' Compute DWLS in parallel across multiple signatures
#'
#' @param TPM_matrix A matrix with TPM normalized counts (genes symbols as rows and samples as columns).
#' @param signatures Path where signatures files are located
#' @param workers Number of processes available to run on parallel.
#'
#' @return A matrix with cell abundance deconvolve with DWLS
#'
#' @export
#'
#' @examples
#'
#' dwls <- computeDWLS_parallel(TPM_matrix, 'src/signatures', workers)
#'
computeDWLS_parallel = function(TPM_matrix, signatures, workers){
  cl = parallel::makeCluster(workers)
  registerDoParallel(cl)

  dwls = foreach (i=1:length(signatures), .combine=cbind) %dopar% {
    source("src/environment_set.R")
    signature <- read.delim(signatures[[i]], row.names=1)
    signature_name = str_split(basename(signatures[[i]]), "\\.")[[1]][1]
    computeDWLS(TPM_matrix, signature, signature_name)
  }

  parallel::stopCluster(cl)
  unregister_dopar()

  return(dwls)
}
