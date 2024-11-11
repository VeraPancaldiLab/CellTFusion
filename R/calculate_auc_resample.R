#' Calculates AUC values for each resampling during cross-validation (CV)
#'
#' @param obs A character vector with observed values (grountruth)
#' @param pred A character vector with predicted values
#'
#' @return A numeric value corresponding to the AUC score
#' @export
#'
#' @examples
#'
#' obs = c("yes", "no", "yes", "no")
#' yes = c("yes", "yes", "yes", "no")
#' AUC = calculate_auc_resample(obs, yes)
#'
calculate_auc_resample = function(obs, pred){

  prob_obs = data.frame("yes" = pred, "obs" = obs)

  prob_obs = prob_obs %>%
    arrange(desc(pred)) %>% #need to be arrange for apply cumulative sum
    mutate(is_yes = (obs == "yes"),
           tp = cumsum(is_yes), #true positive above the threshold - cumulative sum to refer to the threshold
           fp = cumsum(!is_yes), #false positive above the threshold - cumulative sum to refer to the threshold
           fpr = fp/sum(obs == 'no'),
           tpr = tp/sum(obs == 'yes'))

  auc_value = calculate_auc(prob_obs$fpr, prob_obs$tpr)

  return(auc_value)
}
