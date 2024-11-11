get_sensitivity_specificity = function(predictions, observed, ml.model){
  prob_obs = bind_cols(predictions, observed = observed)

  prob_obs = prob_obs %>%
    arrange(desc(yes)) %>% #need to be arrange for apply cumulative sum
    mutate(is_yes = (observed == "yes"),
           tp = cumsum(is_yes), #true positive above the threshold - cumulative sum to refer to the threshold
           fp = cumsum(!is_yes), #false positive above the threshold - cumulative sum to refer to the threshold
           sensitivity = tp/sum(observed == 'yes'),
           fpr = fp/sum(observed == 'no'),
           specificity = 1 - fpr) %>%
    select(sensitivity, specificity, fpr) %>%
    mutate(model = ml.model)

  # starts_at_zero <- any(prob_obs$sensitivity == 0 & prob_obs$fpr == 0)

  # ##Add dummy row if it doesnt start at 0
  # if(!starts_at_zero){
  #   dummy_row <- data.frame(
  #     sensitivity = 0,
  #     specificity = 1,
  #     fpr = 0,
  #     model = ml.model
  #   )
  #
  #   prob_obs = rbind(dummy_row, prob_obs)
  # }

  prob_obs = prob_obs %>%
    mutate(Accuracy = calculate_accuracy(., observed),
           precision = calculate_precision(., observed),
           recall = calculate_recall(., observed))


  return(prob_obs)

}
