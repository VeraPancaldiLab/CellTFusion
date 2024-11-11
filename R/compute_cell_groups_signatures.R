compute_cell_groups_signatures = function(deconv_res, network_res, cell_groups, features, deconvolution_test, TFs_test){

  #Remove colors indicatives from deconvolution features to be able to project them in the raw deconvolution results
  #pattern_colors <- paste0("_(", paste(unique(network_res[[2]]), collapse = "|"), ")$")

  # Use the regular expression to remove the suffix colors from the deconvolution features
  # for (i in 1:length(cell_groups[[2]])) {
  #   cell_groups[[2]][[i]] = gsub(pattern_colors, "", cell_groups[[2]][[i]])
  # }

  ################################################################################Simulate TFs module scores
  TF.matrix_simulated = create_tfs_modules(TFs_test, network_res)
  module_colors = unique(network_res[[2]])

  ################################################################################Simulate cell subgroups
  #Scale deconvolution features by columns for making them comparable between cell types (0-1).
  for (i in 1:ncol(deconvolution_test)) {
    deconvolution_test[,i] = deconvolution_test[,i]/max(deconvolution_test[,i])
  }

  deconv_subgroups <- mapply(c, deconv_res[[3]], deconv_res[[4]], SIMPLIFY = FALSE) #Join cell groups
  iterations = find.maximum.iteration(deconv_subgroups)

  # Create same groups composition
  for (m in 1:iterations) {
    base_groups = list()
    for (i in 1:length(deconv_subgroups)){
      if(length(deconv_subgroups[[i]])!=0){
        idy = grep(paste0("Iteration.",m), names(deconv_subgroups[[i]]))
        if(length(idy)!=0){
          base_groups = append(base_groups, deconv_subgroups[[i]][idy])
        }
      }
    }

    deconv_subgroups_values = c()
    for (i in 1:length(base_groups)) {
      deconv_subgroups_values = cbind(deconv_subgroups_values, rowMedians(as.matrix(deconvolution_test[,base_groups[[i]]]))) #Compute median using base groups
    }
    colnames(deconv_subgroups_values) = names(base_groups)
    deconvolution_test = cbind(deconv_subgroups_values, deconvolution_test) # Join cell subgroups and deconv features

  }

  deconvolution_test = deconvolution_test[,colnames(deconvolution_test)%in%colnames(deconv_res[[1]])]

  # Compute composite scores
  idx = which(names(cell_groups[[2]]) %in% features)
  cell_dendrogram = c()
  names = c()
  for (i in 1:length(idx)) {
    pca_cells = deconvolution_test[,cell_groups[[2]][[idx[i]]]]
    pca_cells <- pca_cells[, apply(pca_cells, 2, function(x) all(!is.na(x)) && var(x, na.rm = TRUE) != 0), drop = FALSE] #Drop zero-columns or NAs
    name_cell_group = names(cell_groups[[2]][idx[i]])
    color = extract_colors(module_colors, name_cell_group)

    if(ncol(pca_cells) > 1){ #Check if there are more than 2 columns
      cell_dendrogram = cbind(cell_dendrogram, compute_composite_score(pca_cells, color, TF.matrix_simulated))
      names = c(names, names(cell_groups[[2]])[idx[i]])
    }
  }

  if(is.null(cell_dendrogram)==T){
    print("No composite scores because all features have zero variance.")
  }else{
    colnames(cell_dendrogram) = names
    rownames(cell_dendrogram) = rownames(TF.matrix_simulated)
  }

  return(data.frame(scale(cell_dendrogram)))
}
