#' Cakculate recall values
#'
#' @param metrics A Dataframe with metrics obtained using get_sensitivity_specificity()
#' @param target A character vector containing the true values from the target variable
#'
#' @return A numeric vector with recall values
#' @export
#'
#' @examples
#'
#' observed = c("yes", "no", "yes", "no", "no")
#' metrics = get_sensitivity_specificity(predict, target, model)
#' recall = calculate_recall(metrics, observed)
#'
calculate_recall <- function(metrics, target) {
  confusion_values <- calculate_confusion_values(metrics, target)
  TP <- confusion_values$TP
  FN <- confusion_values$FN

  # Calculate recall (sensitivity)
  recall <- TP / (TP + FN)

  return(recall)
}
