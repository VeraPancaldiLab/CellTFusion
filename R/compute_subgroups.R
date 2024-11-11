compute_subgroups = function(deconvolution, thres_corr, file_name){
  data = data.frame(deconvolution)
  cell_subgroups = list()
  cell_groups_similarity = list()
  cell_groups_discard = list()
  if (ncol(data) < 2) {
    warning("Deconvolution features with less than two columns for subgrouping (skipping)\n")
    return(list(data, cell_subgroups, cell_groups_similarity, cell_groups_discard))
  }else{
    #################### Proportionality-based correlation
    is_similar <- function(value1, value2, threshold) {return(abs(value1 - value2) <= threshold)}
    similarity_matrix <- matrix(FALSE, nrow = ncol(data), ncol = ncol(data), dimnames = list(names(data), names(data)))
    for (col1 in names(data)) {
      for (col2 in names(data)) {
        similarity <- all(mapply(is_similar, data[[col1]], data[[col2]], MoreArgs = list(0.05))) #similarity threshold = 0.05
        similarity_matrix[col1, col2] <- similarity
      }
    }
    get_upper_tri <- function(cormat){
      cormat[lower.tri(cormat, diag = T)]<- NA
      return(cormat)
    }
    upper_tri <- get_upper_tri(similarity_matrix)
    x <- melt(upper_tri) %>%
      na.omit() %>%
      mutate_all(as.character)
    indice = 1
    subgroup = list()
    vec = unique(x$Var1)
    while(length(vec)>0){
      sub = x[which(x$Var1%in%vec[1]),]
      sub = sub[which(sub$value==T),]
      if(nrow(sub)!=0){
        subgroup[[indice]] = c(vec[1], sub$Var2)
        x = x[-which(x$Var1%in%subgroup[[indice]]),] #Variable 1
        x = x[-which(x$Var2%in%subgroup[[indice]]),] #Variable 2
        vec = vec[-which(vec%in%subgroup[[indice]])]
        indice = indice + 1
      }else{
        indice = indice
        vec = vec[-1]
      }
    }

    if(length(subgroup)!=0){
      for (i in 1:length(subgroup)){
        names(subgroup)[i] = paste0(file_name, "_Subgroup.Similarity.", i) #Name subgroups
      }
      lis = remove_subgroups(subgroup) #Map subgroups with same method
      if(length(lis)>0){
        cell_groups_discard = subgroup[lis]
        subgroup = subgroup[-lis] #Remove subgroups if all subgroupped features belong to the same method
      }

      if(length(subgroup)!=0){  #check if after removal of subgroups with equal method, you still have subgroups
        cell_groups_similarity = subgroup
        data_sub = c()
        for(i in 1:length(cell_groups_similarity)){ #Create data frame with features subgroupped
          sub = data.frame(data[,colnames(data)%in%cell_groups_similarity[[i]]]) #Map features that are inside each subgroup from input (deconvolution)
          sub$median = rowMedians(as.matrix(sub), useNames = FALSE) #Compute median of subgroups across patients
          data_sub = data.frame(cbind(data_sub, sub$median)) #Save median in a new data frame
          colnames(data_sub)[i] = names(cell_groups_similarity)[i]
          name = colnames(data)[which(!(colnames(data)%in%cell_groups_similarity[[i]]))]
          data = data[,-which(colnames(data)%in%cell_groups_similarity[[i]])] #Remove from deconvolution features that are subgrouped
          if(ncol(data.frame(data))==1){
            data = as.data.frame(data)
            colnames(data)[1] = name
          }
        }

        rownames(data_sub) = rownames(data) #List of patients
        data_sub = data.frame(data_sub[,colnames(data_sub)%in%names(cell_groups_similarity)])
        colnames(data_sub) = names(cell_groups_similarity)

        data = cbind(data, data_sub) #Join subgroups in deconvolution file
      }else{
        cell_groups_similarity = list()
      }
      k = 2
    }else{
      k = 3
    }
    if(ncol(data) == 1){ #everything is already subgroupped
      return(list(data, cell_subgroups, cell_groups_similarity, cell_groups_discard))
    }

    #################### Linear-based correlation
    if(k==2 | k==3){
      terminate = FALSE
      iteration = 1
      while (terminate == FALSE) {
        corr_df <- correlation(data.matrix(data))
        vec = colnames(data)
        indice = 1
        subgroup = list()
        data_sub = c()
        while(length(vec)>0){ #Keep running until no features are left
          if(vec[1] %in% corr_df$measure1){ #Check if feature still no-grouped
            tab = corr_df[corr_df$measure1 == vec[1],] #Take one feature against the others
            tab = tab[tab$r>thres_corr,] #Select features corr above the threshold
            if(nrow(tab)!=0){ #If algorithm found features above corr
              subgroup[[indice]] = c(vec[1], tab$measure2) #Save features as subgroup
              idx = which(corr_df$measure1 %in% subgroup[[indice]])
              if(length(idx)>0){corr_df = corr_df[-idx,]} #Remove features already subgroupped
              idy = which(corr_df$measure2 %in% subgroup[[indice]])
              if(length(idy)>0){corr_df = corr_df[-idy,]} #Remove features already subgroupped
              vec = vec[-which(vec%in%subgroup[[indice]])] #Remove feature already subgroupped from vector
              indice = indice + 1
            }else{ #Condition when there is no correlation above the threshold (features no subgroupped)
              corr_df = corr_df[-which(corr_df$measure1 == vec[1]),] #Remove variable from corr matrix to keep subgrouping the others
              if(length(which(corr_df$measure2==vec[1]))>0){corr_df = corr_df[-which(corr_df$measure2 == vec[1]),]}
              vec = vec[-1] #Remove variable from vector to keep analyzing the others
              indice = indice #Not increase index cause no subgroup appeared
            }
          }else{ #If feature is not in corr matrix it means that there is no any significant correlation against it and other features
            vec = vec[-1] #Remove variable from vector to keep analyzing the others
            indice = indice  #Not increase index cause no subgroup appeared
          }
        }

        if(length(subgroup)!=0){
          for (i in 1:length(subgroup)){ #Name subgroups
            names(subgroup)[i] = paste0(file_name, "_Subgroup.", i, ".Iteration.", iteration)
          }
          ###Check whenever some subgroups belong to the same method
          if(iteration == 1){
            idx = remove_subgroups(subgroup) #Map subgroups with same method
            if(length(idx)>0){
              if(length(cell_groups_discard)>0){
                cell_groups_discard = c(cell_groups_discard, subgroup[idx])
                duplica = which(duplicated(cell_groups_discard)) #check if there are subgroups duplicated discarded
                if(length(duplica)>0){
                  cell_groups_discard = cell_groups_discard[-duplica]
                }
              }
              else{
                cell_groups_discard = subgroup[idx]
              }
              subgroup = subgroup[-idx] #Remove subgroups if all subgroupped features belong to the same method
            }
          }

          if(length(subgroup)!=0){ #check if after removal of subgroups with equal method, you still have subgroups (when iteration == 1)
            #Take median expression of subgroups
            for(i in 1:length(subgroup)){ #Create data frame with features subgroupped
              sub = data.frame(data[,colnames(data)%in%subgroup[[i]]]) #Map features that are inside each subgroup from input (deconvolution)
              sub$median = rowMedians(as.matrix(sub), useNames = FALSE) #Compute median of subgroup across patients
              data_sub = data.frame(cbind(data_sub, sub$median)) #Save median in a new data frame
              colnames(data_sub)[i] = names(subgroup)[i]
              name = colnames(data)[which(!(colnames(data)%in%subgroup[[i]]))]
              data = data.frame(data[,-which(colnames(data)%in%subgroup[[i]])]) #Remove from deconvolution features that are subgrouped
              if(ncol(data.frame(data))==1){
                data = as.data.frame(data)
                colnames(data)[1] = name
              }
            }

            rownames(data_sub) = rownames(data) #List of patients

            if(iteration == 1){ #Save what is inside the first subgroups
              cell_subgroups = subgroup
              data_sub = data.frame(data_sub[,colnames(data_sub)%in%names(cell_subgroups)])
              colnames(data_sub) = names(cell_subgroups)
            }else{
              for (i in 1:length(subgroup)) {
                cell_subgroups[[length(cell_subgroups)+1]] = subgroup[[i]]
                names(cell_subgroups)[length(cell_subgroups)] = names(subgroup)[i]
              }
            }

            if(ncol(data)!=0){
              data = cbind(data, data_sub)
            }else{
              data = data_sub #data will be 0 if all deconvolution features have been subgroupped
              terminate = TRUE
            }
            iteration = iteration + 1
          }else{
            terminate = TRUE #when the only subgroup that keep grouping is composed from the same method
          }

        }else{
          terminate = TRUE
        }
      }

      if(is.null(data_sub)==FALSE){
        data = cbind(data, data_sub)
      }

      idx = which(duplicated(t(data)))
      if(length(idx)>0){
        names = colnames(data)[idx]
        data = data[,-idx, drop = F]
        colnames(data) = names
      }
    }

    return(list(data, cell_subgroups, cell_groups_similarity, cell_groups_discard))
  }

}
