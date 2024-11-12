#' Create TFs modules
#'
#' This function re-create existing TF modules on a different TF activity matrix.
#'
#' @param TF.matrix TFs activity matrix with samples as rows and TFs as columns.
#' @param network_tfs A TF network object obtained from compute.WTCNA() from which TF modules need to be re-create.
#'
#' @return A matrix of TFs modules scores across samples
#' @export
#'
#' @examples
#'
#' TF.matrix_simulated = create_tfs_modules(TFs_test, network_res)
#'
create_tfs_modules = function(TF.matrix, network_tfs){

  tfs.modules = TF.matrix %>%
    t() %>%
    data.frame() %>%
    mutate(Module = "na")

  for (i in 1:length(network_tfs[[3]])) {
    tfs.modules$Module[which(rownames(tfs.modules) %in% network_tfs[[3]][[i]])] = names(network_tfs[[3]])[i]
  }

  tfs_colors = tfs.modules %>%
    pull(Module)

  MEList = moduleEigengenes(TF.matrix, colors = tfs_colors, scale = F) #Data already scale
  MEs = MEList$eigengenes
  MEs = orderMEs(MEs)

  return(MEs)
}
