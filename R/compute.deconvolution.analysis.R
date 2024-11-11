compute.deconvolution.analysis <- function(deconvolution, corr, zero = 0.9, high_corr = 0.9, seed = NULL, return = T){
  deconvolution.mat = deconvolution

  #####Unsupervised filtering

  #Remove NA (this need to be check -- not possible to have NAs values in deconv)
  deconvolution.mat <- deconvolution.mat %>%
    mutate(across(everything(), ~ replace_na(.x, 0)))

  #Remove high zero number features
  cat(paste0("Removing features with high zero number ", round(zero*100,2), "%...............................................................\n\n"))
  deconvolution.mat = deconvolution.mat[, colSums(deconvolution.mat == 0, na.rm=TRUE) < round(zero*nrow(deconvolution.mat)) , drop=FALSE]
  diff_colnames <- setdiff(colnames(deconvolution), colnames(deconvolution.mat))
  zero_features <- deconvolution[, diff_colnames]

  #Remove low_variance features
  variance = remove_low_variance(deconvolution.mat, plot = return)
  deconvolution.mat = variance[[1]]
  low_variance_features = variance[[2]]

  #Scale deconvolution features by columns for making them comparable between cell types (0-1).
  cat("Scaling deconvolution features for comparison between cell types...............................................................\n\n")
  for (i in 1:ncol(deconvolution.mat)) {
    deconvolution.mat[,i] = deconvolution.mat[,i]/max(deconvolution.mat[,i])
  }

  #####Cell types split
  cat("Splitting deconvolution features per cell type...............................................................\n\n")
  cells_types = compute.cell.types(deconvolution.mat)
  cells = cells_types[[1]]
  cells_discarded = cells_types[[2]]

  ######Pairwise correlation filtering (Highly correlated variables >0.9) within cell types
  cat("Finding group of features with high correlation between each other...............................................................\n\n")
  features_high_corr = list()
  j = 1
  for (i in 1:length(cells)) {
    data = cells[[i]]
    if(is.null(ncol(data))==T){
      cells[[i]] = data
    }else if(ncol(data)>1){
      data = removeCorrelatedFeatures(data, high_corr, names(cells)[i], seed)
      cells[[i]] = data[[1]]
      if(length(data[[2]])>0 && is.null(data[[3]])==F){
        features_high_corr[[j]] = data[[2]]
        names(features_high_corr)[j] = data[[3]]
        j = j+1
      }
    }
  }

  #####Subgrouping of deconvolution features
  res = list()
  groups = list()
  groups_similarity = list()
  groups_discard = list()
  for (i in 1:length(cells)) {
    x = compute_subgroups(cells[[i]], file_name = names(cells)[i], thres_corr = corr)
    res = c(res, x[1])
    groups = c(groups, x[2])
    groups_similarity = c(groups_similarity, x[3])
    groups_discard = c(groups_discard, x[4])
  }
  names_cells = c("B.cells", "B.naive", "B.memory", "Macrophages.cells", "Macrophages.M0", "Macrophages.M1", "Macrophages.M2", "Monocytes", "Neutrophils", "NK.cells", "NK.activated", "NK.resting", "NKT.cells", "CD4.cells", "CD4.memory.activated",
                  "CD4.memory.resting", "CD4.naive", "CD8.cells", "T.cells.regulatory", "T.cells.non.regulatory","T.cells.helper", "T.cells.gamma.delta", "Dendritic.cells", "Dendritic.activated", "Dendritic.resting", "Cancer", "Endothelial",
                  "Eosinophils", "Plasma.cells", "Myocytes", "Fibroblasts", "Mast.cells", "Mast.activated", "Mast.resting", "CAF")

  names(res) = names_cells
  names(groups) = names_cells
  names(groups_similarity) = names_cells
  names(groups_discard) = names_cells

  #####Preparing output
  dt = c()
  for (i in 1:length(res)) {
    dt = c(dt, res[[i]])
  }
  dt = data.frame(dt)
  rownames(dt) = rownames(deconvolution.mat)

  #####Create and export table with subgroups

  #Count number of subgroups - Linear-based
  idx = c()
  for (i in 1:length(groups)){
    if(length(groups[[i]])>0){
      for (j in 1:length(groups[[i]])){
        idx = c(idx, names(groups[[i]])[[j]])
      }
    }
  }
  data.groups = data.frame(matrix(nrow = length(idx), ncol = 2)) #Create table
  colnames(data.groups) = c("Cell_subgroups", "Methods-signatures")
  data.groups$Cell_subgroups = idx #Assign subgroups

  #Save methods corresponding to each subgroup
  contador = 1
  for (i in 1:length(groups)){
    if(length(groups[[i]])>0){
      for (j in 1:length(groups[[i]])){
        data.groups[contador,2] = paste(groups[[i]][[j]], collapse ="\n")
        contador = contador + 1
      }
    }
  }


  #Count number of subgroups - Proportionality-based
  idy = c()
  for (i in 1:length(groups_similarity)){
    if(length(groups_similarity[[i]])>0){
      for (j in 1:length(groups_similarity[[i]])){
        idy = c(idy, names(groups_similarity[[i]])[[j]])
      }
    }
  }
  data.groups.similarity = data.frame(matrix(nrow = length(idy), ncol = 2)) #Create table
  colnames(data.groups.similarity) = c("Cell_subgroups", "Methods-signatures")
  data.groups.similarity$Cell_subgroups = idy #Assign subgroups

  #Save methods corresponding to each subgroup
  contador = 1
  for (i in 1:length(groups_similarity)){
    if(length(groups_similarity[[i]])>0){
      for (j in 1:length(groups_similarity[[i]])){
        data.groups.similarity[contador,2] = paste(groups_similarity[[i]][[j]], collapse ="\n")
        contador = contador + 1
      }
    }
  }

  #Save data to export
  if(return){
    data.output = rbind(data.groups.similarity, data.groups)
    write.csv(dt, 'Results/Deconvolution_after_subgrouping.csv')
    write.csv(data.output, 'Results/Cell_subgroups.csv', row.names = F)
  }

  message("Deconvolution features subgroupped")

  results = list(dt, res, groups, groups_similarity, groups_discard, zero_features, low_variance_features, cells_discarded, features_high_corr)
  names(results) = c("Deconvolution matrix", "Deconvolution groups per cell types", "Deconvolution groups - Linear-based correlation", "Deconvolution groups - Proportionality-based correlation",
                     "Discarded groups with equal method", "Discarded features with high number of zeros", "Discarded features with low variance", "Discarded cell types",
                     "High correlated deconvolution groups (>0.9) per cell type")
  return(results)

}
