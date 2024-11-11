remove.cell.groups.corr <- function(data, colors, threshold = 0.9) {

  features_high_corr = c()
  # Compute correlation matrix
  corr_matrix <- cor(data[[1]])
  # Find highly correlated features
  contador = 1
  while(nrow(corr_matrix)>0){
    color_features = c()
    feature = data.frame(corr_matrix[1, , drop = FALSE]) #Extract first row feature
    feature = feature %>%                                #Take only high corr above threshold
      mutate_all(~ifelse(. > threshold, ., NA)) %>%
      select_if(~all(!is.na(.)))

    corr_matrix = corr_matrix[-which(rownames(corr_matrix)%in%colnames(feature)),-which(colnames(corr_matrix)%in%colnames(feature)), drop = F] #Remove already joined features
    color_features = list()
    if(ncol(feature)>1){
      for (m in 1:ncol(feature)) {
        color_group = extract_colors(colors, colnames(feature)[m])
        color_features[[m]] = color_group
      }

      all_equal <- all(sapply(color_features, function(x) identical(x, color_features[[1]]))) #Check whether highly correlated features belong to the same group of TFs modules

      if(all_equal == TRUE){
        new_group_composition = unique(unlist(unname(data[[2]][colnames(feature)])))
        new_group_value = rowMeans(data[[1]][,colnames(data[[1]])%in%colnames(feature)])

        message("Highly correlated features (r>", threshold,"): ", paste(colnames(feature), collapse = ', '), ". Combining.")

        if(contador==1){
          #Remove features from original data
          new_data <- data[[1]][, -which(colnames(data[[1]])%in%colnames(feature)), drop = F]
          new_groups = data[[2]][-which(names(data[[2]]) %in% colnames(feature))]
        }else{
          new_data <- new_data[, -which(colnames(new_data)%in%colnames(feature)), drop = F]
          new_groups = new_groups[-which(names(new_groups) %in% colnames(feature))]
        }

        #Add new combined features
        new_name = paste0("Dendrogram_",  paste0(color_group, collapse = "_"), ".group_combined_", contador)
        new_data = cbind(new_data, new_group_value)
        colnames(new_data)[length(new_data)] = new_name

        new_groups[[length(new_groups)+1]] = new_group_composition
        names(new_groups)[length(new_groups)] = new_name

        contador = contador + 1
      }else{
        message("Highly correlated features between different clusters of TFs modules (r>", threshold,"): ", paste(colnames(feature), collapse = ', '), ". Not joining.")
      }
    }else{
      if(contador == 1){
        new_data = data[[1]]
        new_groups = data[[2]]
      }
    }
  }

  res = list(new_data, new_groups)

  return(res)
}
