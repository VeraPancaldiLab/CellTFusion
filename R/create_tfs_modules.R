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
