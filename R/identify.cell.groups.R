identify.cell.groups = function(features, tfs.modules.groups, cor_type = "p", clustering.method = "ward.D2", distance.method = "euclidean", width = 12, height = 18, return = T){

  moduleTraitCor = features[[1]]

  names(features[[2]]) = paste0("ME", names(features[[2]])) #To match names of columns from corr matrix

  #Gather significant features across TF modules clusters
  features_vec = list()
  for (i in 1:length(tfs.modules.groups)) {
    features_vec[[i]] = unique(unlist(unname(features[[2]][tfs.modules.groups[[i]]])))
    names(features_vec)[i] = names(tfs.modules.groups)[i]
  }

  lis.dendrogram = list()

  for (i in 1:length(features_vec)){
    TFmoduleTraitcor = moduleTraitCor[,colnames(moduleTraitCor)%in%features_vec[[i]]]
    TFmoduleTraitcor = TFmoduleTraitcor[rownames(TFmoduleTraitcor)%in%tfs.modules.groups[[i]], , drop=F]
    data_scaled = scale(t(TFmoduleTraitcor)) #Scale per TF module
    ###Dendogram by Module
    d <- dist(data_scaled, method = distance.method)
    d = d/sqrt(ncol(data_scaled)) #Adjust/Scale distance matrix for number of features to make dendrograms comparable
    dendrogram <- hclust(d, method = clustering.method)
    if(return){
      pdf(paste0("Results/Dendogram_cell_types_", names(features_vec)[[i]]), width = width, height = height)
      par(mar = c(5, 2, 4, 35)) #bottom, left, top, right
      plot(as.dendrogram(dendrogram), horiz= T)
      dev.off()
    }
    lis.dendrogram[[i]] = dendrogram
  }

  names(lis.dendrogram) = names(features_vec)

  # for (i in 1:length(features_vec)){
  #   moduleTraitCor[i,!(colnames(moduleTraitCor)%in%features_vec[[i]])] = 0 ###Set to 0 not significant correlations
  # }

  #Add dendrogram "all" considering all TFs modules

  # data_scaled = scale(t(moduleTraitCor)) #Scale per TF module
  # d <- dist(data_scaled, method = "euclidean")
  # d = d/sqrt(ncol(data_scaled))
  # dendrogram_all <- hclust(d, method = clustering.method)
  # if(return){
  #   pdf("Results/Dendogram_cell_types_all", width = width, height = height)
  #   par(mar = c(5, 2, 4, 35)) #bottom, left, top, right
  #   plot(as.dendrogram(dendrogram_all), horiz= T)
  #   dev.off()
  # }
  #
  # lis.dendrogram[[length(lis.dendrogram)+1]] = dendrogram_all
  # names(lis.dendrogram)[[length(lis.dendrogram)]] = "all"

  return(lis.dendrogram)

}
