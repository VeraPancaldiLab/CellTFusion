compute_composite_score = function(cell_group, color_group, tfs.module.matrix, prop_var = 0.7){

  module_group = paste0("ME", color_group) #To match with columns of TFs modules
  pca_group <- cbind(scale(tfs.module.matrix[,module_group, drop=F]), scale(cell_group)) #Combined TF module corresponding to each cell group
  svd_result <- svd(pca_group, nu = min(nrow(pca_group), ncol(pca_group)), nv = min(nrow(pca_group), ncol(pca_group))) #Performs Singular Value Decomposition (SVD)
  singular_values <- svd_result$d #Extract singular values
  nPCs <- sum(cumsum(singular_values^2) / sum(singular_values^2) < prop_var) + 1 #Extract number of PCs necessary to cover the prop_var
  variance_explained <- (singular_values^2/sum(singular_values^2)) #Extract the % variance explained by each PC
  weights <- variance_explained[1:nPCs]
  PCs <- svd_result$u #Take components
  aveg_group = rowMeans(pca_group) #Consider average expression to aligned direction of components

  for (i in 1:ncol(PCs)) {
    cor_group = cor(aveg_group, PCs[,i], use = "p")

    # Check if the correlation is finite; if not, set it to 0
    if (!is.finite(cor_group)) {
      cor_group = 0
    }

    if(cor_group<0){
      PCs[,i] = PCs[,i]*-1
    }
  }

  #selected_components = minMax(PCs[,1:nPCs]) # avoiding negatives (refers only to the direction)
  selected_components = PCs[,1] # Taking eigenvalue

  # if(nPCs != 1){
  #   composite_score = rowSums(selected_components * weights)
  # }else{
  #   composite_score = selected_components * weights
  # }

  #composite_score = scale(composite_score)
  composite_score = scale(selected_components)

  return(composite_score)
}
