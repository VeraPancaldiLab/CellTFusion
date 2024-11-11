#' Compute weighted feature importance from base models and meta-learner for stacking models
#'
#' @param base_importance A matrix of feature importance from base models
#' @param base_models A character vector with chosen base models
#' @param meta_learner A caret object with the meta-learner model trained using base models
#'
#' @return A matrix of feature importance weighted from the base models and the meta-learner
#' @export
#'
#' @examples
#'
#' var_importance = calculate_feature_importance_stacking(variable_importance, base_models, meta-learner)
#'
calculate_feature_importance_stacking = function(base_importance, base_models, meta_learner){

  #Extract features importance values within each base model for the meta-learner
  base_importance_list = list()
  for (i in 1:length(base_models)) {
    check = ncol(base_importance[[base_models[i]]][["importance"]])
    if(check > 1){ #Means importance is given for each class
      base_importance_list[[i]] = base_importance[[base_models[i]]][["importance"]] %>%
        rownames_to_column("features") %>%
        dplyr::select(features, yes) %>% #Take only importance for positive class
        dplyr::rename(importance = yes)
    }else{
      base_importance_list[[i]] = base_importance[[base_models[i]]][["importance"]] %>%
        rownames_to_column("features") %>%
        dplyr::rename(importance = Overall)
    }
    names(base_importance_list)[i] = base_models[i]
  }

  #Combine all base model importances in one data frame and add the model name
  combined_importance <- bind_rows(
    lapply(names(base_importance_list), function(model) {
      base_importance_list[[model]] %>%
        data.frame() %>%
        mutate(model = model)
    })
  )

  #Calculate base-models importance for the meta-learner
  meta_importance = varImp(meta_learner, scale = F)$importance %>%
    rownames_to_column("model")

  #Normalize the meta-learner's importance scores so they sum to 1
  meta_importance$Overall <- meta_importance$Overall / sum(meta_importance$Overall)

  #Combine features importance within base models with the overall importance for meta-learner
  combined_importance <- combined_importance %>%
    left_join(meta_importance, by = "model") %>%
    mutate(weighted_importance = importance * Overall) # importance is from base, Overall is from meta

  #Sum the weighted importance by feature across all models
  final_importance <- combined_importance %>%
    group_by(features) %>%
    summarise(final_importance = sum(weighted_importance, na.rm = TRUE)) %>%
    arrange(desc(final_importance))

  return(final_importance)

}
