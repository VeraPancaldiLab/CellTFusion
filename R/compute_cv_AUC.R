#' Compute AUC values for each machine learning model
#'
#' @param models List of trained models
#' @param file_name (Optional) File name for plot
#' @param base_models Boolean value to specify if base models for stacking need to be chosen or not
#' @param return Boolean value to specify if plot the AUC values across models should be saved or not
#'
#' @return A list containing
#'
#' - AUC values for each model
#' - Top model with best AUC
#' - If base_models = T, it returns a character vector with the chosen base models (see choose_base_models())
#'
#' @export
#'
#' @examples
#'
#' res = compute_cv_AUC(ml_models, base_models = T, file_name = "Test", return = T)
#'
compute_cv_AUC = function(models, file_name = NULL, base_models = F, return = T){

  #Bind AUC values from each model
  auc = list()
  for (i in 1:length(models)){
    auc[[i]] = models[[i]]$resample %>% #we use the resample matrix and not directly the results matrix as some have hyperparameters so we will need to define best on the tuned parameter (=more code) - resample matrix is made based on the best tuning
      mutate(model = names(models)[i])
    names(auc)[i] = names(models)[i]
  }
  auc_data = do.call(rbind, auc)

  #Retrieve top model based on accuracy
  res_auc <- auc_data %>%
    group_by(model) %>%
    summarise(AUC = mean(AUC))  %>%
    arrange(desc(AUC))

  top_model = res_auc %>%
    slice(1) %>%
    pull(model)

  if(return){
    pdf(paste0("Results/AUC_CV_methods_", file_name, ".pdf"), width = 10)
    plot(ggplot(auc_data, aes(x = model, y = AUC, fill = model)) +
           geom_boxplot() +
           labs(title = "Distribution of AUC scores by Model",
                x = "Model",
                y = "AUC") +
           theme_minimal() +
           theme(legend.position = "none"))
    dev.off()
  }

  if(base_models == T){
    cat("Choosing base models for stacking.......................................\n\n")
    base_models = choose_base_models(models, metric = "AUC")
    cat("Models chosen are:", paste0(base_models, collapse = ", "), "\n\n")
    return(list("AUC" = res_auc, "Top_model" = top_model, "Base_models" = base_models))
  }else{
    return(list("AUC" = res_auc, "Top_model" = top_model))
  }

}
