removeCorrelatedFeatures <- function(data, threshold, name, n_seed) {

  features_high_corr = c()
  cell_name = c()
  # Compute correlation matrix
  corr_matrix <- cor(data)
  # Find highly correlated features
  contador = 1
  while(nrow(corr_matrix)>0){
    set.seed(n_seed)
    feature = data.frame(corr_matrix[1, , drop = FALSE]) #Extract first row feature
    feature = feature %>%                                #Take only high corr above threshold
      mutate_all(~ifelse(. > threshold, ., NA)) %>%
      select_if(~all(!is.na(.)))

    corr_matrix = corr_matrix[-which(rownames(corr_matrix)%in%colnames(feature)),-which(colnames(corr_matrix)%in%colnames(feature)), drop = F] #Remove already joined features

    if(ncol(feature)>1){
      keep = colnames(feature)[sample(ncol(feature), size = 1)] #From high corr group, keep only one feature

      print(paste0("Highly correlated features (r>", threshold,"): ", paste(colnames(feature), collapse = ', ')))
      cat(paste0("Keeping only feature: ", keep, "\n\n"))

      if(length(features_high_corr)>0){
        features_high_corr = c(features_high_corr, colnames(feature))
      }else{
        features_high_corr = colnames(feature)
      }

      feature = feature[,-which(colnames(feature)%in%keep), drop = F]

      if(contador==1){
        new_data <- data[, -which(colnames(data)%in%colnames(feature)), drop = F] #Remove rest of the features from original data
      }else{
        new_data <- new_data[, -which(colnames(new_data)%in%colnames(feature)), drop = F] #Remove rest of the features from original data
      }
      contador = contador + 1
      cell_name = name
    }else{
      if(contador == 1){
        new_data = data
      }else{ #If it already started the loop
        new_data = new_data
      }
    }
  }

  if(length(cell_name)==0){
    cell_name = NULL
  }

  return(list(new_data, features_high_corr, cell_name))
}
