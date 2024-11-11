compute.WTCNA <- function(TFs.matrix, network.type = "unsigned", clustering.method = "ward.D2", minMod = 15, corr_mod = 0.9, cor_type = "p", return = T){

  cat("Creating weighted TF-coactivity network......................................................................\n\n")
  #####Choose parameter for scale-free network topology
  powers = c(c(1:10), seq(from = 12, to=20, by=1))
  sft = pickSoftThreshold(TFs.matrix, powerVector = powers, verbose = 0, networkType = network.type)

  if(return){
    pdf("Results/Soft Threshold")
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
         main = paste("Scale independence"))
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=0.9,col="red");
    abline(h=0.90,col="red")
    dev.off()
  }

  #Automatic choosing of soft-threshold
  target = 0.9 #Target SFT.R.sq value
  diff = abs(-sign(sft$fitIndices[,3])*sft$fitIndices[,2] - target) #Calculate absolute difference
  min_index = which.min(diff) #Identify the index with the minimum difference
  softPower = powers[min_index]
  cat("Choosing", softPower, "as soft-threshold......................................................................\n\n")

  #####Co-expression matrix using nodes adjacency and topological overlapping nodes
  cat("Calculating nodes adjacency and topological overlapping nodes.................................................\n\n")
  adjacency = adjacency(TFs.matrix, power =softPower, type=network.type, corFnc = "cor", corOptions = list(use = cor_type))
  TOM = TOMsimilarity(adjacency, TOMType = network.type)
  dissTOM = 1-TOM

  #####Unsupervised hierarchical clustering using dissimilarity matrix
  geneTree = hclust(dist(dissTOM), method = clustering.method)
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = minMod);
  dynamicColors = labels2colors(dynamicMods)

  #####Remove variables and clean garbage
  rm(adjacency, TOM, dissTOM)
  gc()

  if(return){
    pdf("Results/Gene dendrogram and module colors")
    plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05,
                        main = "Gene dendrogram and module colors")
    dev.off()
  }

  #####Calculate eigenvectors from modules
  cat("Calculating eigenvectors from modules.................................................\n\n")
  MEList = moduleEigengenes(TFs.matrix, colors = dynamicColors, scale = F) #Data already scale
  MEs = MEList$eigengenes
  MEs = orderMEs(MEs)

  print(paste0("Merging modules significantly correlated with ", corr_mod, "........"))
  merge = mergeModules(MEs, dynamicColors, corr_mod)
  MEs = merge[[1]]
  dynamicColors = merge[[2]]
  modtfs = list()
  for(i in 1:ncol(MEs)){
    tfs = colnames(TFs.matrix)
    modules = c(substring(names(MEs), 3))
    inModule = is.finite(match(dynamicColors,modules[i]))
    modtfs[[i]] = tfs[inModule]
  }
  names(modtfs) = modules

  if(return){
    pdf("Results/Gene dendrogram and module colors after merging")
    plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05,
                        main = "Gene dendrogram and module colors")
    dev.off()
  }

  TFspropVar = propVarExplained(TFs.matrix, dynamicColors, MEs, corFnc = "cor", corOptions = "use = 'p'")

  output = list(MEs, dynamicColors, modtfs, TFspropVar)
  names(output) = c("TFs module matrix", "TFs colors", "TFs per module", "Proportion of variance")

  contador = 1
  tfs_modules = data.frame(matrix(nrow = length(modtfs), ncol = 2))
  colnames(tfs_modules) = c("TFs module", "Composition")
  for (i in 1:length(modtfs)) {
    tfs_modules[contador,1] = names(modtfs)[i]
    tfs_modules[contador,2] = paste(modtfs[[i]], collapse = ",")
    contador = contador + 1
  }

  if(return){
    write.csv(tfs_modules, 'Results/TFs_modules.csv', row.names = F)
  }

  return(output)

}
