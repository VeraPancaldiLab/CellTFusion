#' Choose three base models for stacking based either in Accuracy or AUC scores
#'
#' @param models List with trained machine learning models
#' @param metric Metric to choose the top base models (either Accuracy or AUC). Default is Accuracy
#'
#' @return A character vector with the base models
#' @export
#'
#' @examples
#'
#' models <- list(BAG = fit.treebag,RF = fit.rf, C50 = fit.c50, GLM = fit.glm, LDA = fit.lda, KNN = fit.knn, CART = fit.cart)
#' base_models = choose_base_models(models, metric = "Accuracy")
#'
choose_base_models = function(models, metric = "Accuracy"){

  #Bind accuracy values from each model
  resample_df = list()
  for (i in 1:length(models)){
    resample_df[[i]] = models[[i]]$resample %>%
      mutate(model = names(models)[i])
    names(resample_df)[i] = names(models)[i]
  }
  resample_df = do.call(rbind, resample_df)

  if(metric == "Accuracy"){
    #Prepare data frame for ploting
    resample_df <- resample_df %>%
      group_by(model) %>%
      summarise(Accuracy = mean(Accuracy))
  }else if(metric == "AUC"){
    #Prepare data frame for ploting
    resample_df <- resample_df %>%
      group_by(model) %>%
      summarise(AUC = mean(AUC))
  }

  resample_df <- resample_df %>%
    mutate(Category = case_when(
      model %in% c("BAG", "C50", "CART", "RF") ~ "Tree-based Methods",
      model %in% c("GLM", "LDA", "GLMNET", "LASSO", "RIDGE") ~ "Linear Models",
      model %in% c("KNN", "SVM_linear", "SVM_radial") ~ "Instance-based Methods",
      TRUE ~ "Other"  # In case there are models not in the above lists
    ))

  if(metric == "Accuracy"){
    groupped_df <- resample_df %>%
      group_by(Category) %>%
      filter(Accuracy == max(Accuracy)) %>%
      ungroup()
  }else if(metric == "AUC"){
    groupped_df <- resample_df %>%
      group_by(Category) %>%
      filter(AUC == max(AUC)) %>%
      ungroup()
  }else{
    stop("No metric defined")
  }

  #Retrieve top model based on accuracy/auc
  base_models <- groupped_df %>%
    pull(model)

  return(base_models)
}
