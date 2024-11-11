cell.groups.computation = function(deconvolution, tfs.module.network, cell.dendrograms, return = T){

  cuts = calculate_dendrogram_cuts(cell.dendrograms) #Identify dendrograms cuts
  module_colors = unique(tfs.module.network[[2]]) #Identify module colors

  cell.groups = list()
  for(k in 1:length(cell.dendrograms)){
    groups_cut = list()
    for (i in 1:length(cuts[[k]])){
      clusters <- dendextend::cutree(cell.dendrograms[[k]], h = cuts[[k]][i], order_clusters_as_data=FALSE)
      y = list() #Store cell groups composition
      x = list() #Store cell groups scores
      for (j in 1:length(table(clusters))) {
        cells = names(clusters)[clusters==j]
        y[[j]] = cells
        ###################################################Compute score for each cell group
        if(length(cells)>1){
          pca_group = deconvolution[,colnames(deconvolution) %in% cells, drop = F]
          color = extract_colors(module_colors, names(cell.dendrograms)[k])
          x[[j]] <- compute_composite_score(pca_group, color, tfs.module.network[[1]])
        }else{
          x[[j]] = deconvolution[,colnames(deconvolution) %in% cells]
        }
      }

      ###################################################Remove groups with one feature
      groups_cut[[i]] = remove_single_groups(x, y)
    }

    cell.groups.values <- lapply(groups_cut, function(x) x[[1]])
    cell.groups.names <- lapply(groups_cut, function(x) x[[2]])

    ###################################################Remove groups with equal composition (IMPORTANT TO MAKE IT WITHIN DENDROGRAM - Cell groups can have same composition but belong to a different TF module
    cell.groups[[k]] = remove_equal(cell.groups.values, cell.groups.names)
  }

  ###################################################Naming of cell groups
  for (i in 1:length(cell.groups)) {
    for (j in 1:length(cell.groups[[i]][[1]])) {
      names(cell.groups[[i]][[1]])[j] = paste0("Dendrogram_", names(cell.dendrograms)[i], "_group_", j)
      names(cell.groups[[i]][[2]])[j] = paste0("Dendrogram_", names(cell.dendrograms)[i], "_group_", j)
    }
  }

  cell.groups.values <- lapply(cell.groups, function(x) x[[1]])
  cell.groups.names <- lapply(cell.groups, function(x) x[[2]])

  cell.groups = list(unlist(cell.groups.values, recursive = F), unlist(cell.groups.names, recursive = F))
  cell.groups[[1]] = cell.groups[[1]] %>%
    data.frame() %>%
    mutate("Samples" = rownames(deconvolution)) %>%
    column_to_rownames("Samples")

  cell.groups = remove.cell.groups.corr(cell.groups, module_colors, threshold = 0.9)

  cat("Removing low variance features (if present)........................\n")
  zero = nearZeroVar(cell.groups[[1]], saveMetrics = TRUE)
  cell.groups[[1]] = cell.groups[[1]][, !zero$nzv]
  cell.groups[[2]] = cell.groups[[2]][!zero$nzv]

  ###################################################Export clusters in a table and save
  clusters = data.frame(matrix(nrow = length(cell.groups[[2]]), ncol = 2))
  for (i in 1:length(cell.groups[[2]])) {
    clusters[i,1] = names(cell.groups[[2]])[[i]]
    clusters[i,2] = paste(cell.groups[[2]][[i]], collapse ="\n")
  }
  colnames(clusters) = c("Cell groups", "Methods-signatures")

  if(return){
    write.csv(clusters, "Results/Cell.groups.composition.csv", row.names = F)
    write.csv(cell.groups[[1]], "Results/Cell.groups.scores.csv")
  }

  return(cell.groups)
}
