mergeModules = function(data, colors, corr){
  df = correlation(data)
  idx = which(round(df$r,2) > corr)
  if(length(idx)>0){
    for(i in seq(1, length(idx), by=2)){
      module1 = df$measure1[idx[i]]
      module2 = df$measure2[idx[i]]
      if((module1 %in% colnames(data)) && (module2 %in% colnames(data))){
        colors[which(colors%in%c(substring(module1, 3), substring(module2, 3)))] = substring(module1, 3)
        data <- data %>%
          mutate(new_column = rowMeans(dplyr::select(., module1, module2))) %>%
          dplyr::rename(module1 = new_column) %>%
          dplyr::select(., -module1, -module2)
      }
    }
  }

  return(list(data, colors))
}
