#' Calculate confusion values
#'
#' @param metrics A Dataframe with metrics obtained using get_sensitivity_specificity()
#' @param target A character vector containing the true values from the target variable
#'
#' @return A list containing the True positives (TP), False negatives (FN), True negatives (TN) and False positives (FP)
#' @export
#'
#' @examples
#'
#' target = c("yes", "no", "yes", "no", "no")
#' metrics = get_sensitivity_specificity(predict, target, model)
#' confusion_values = calculate_confusion_values(metrics, target)
#'
calculate_confusion_values <- function(metrics, target) {
  sensitivity <- metrics[,"sensitivity"]
  specificity <- metrics[,"specificity"]

  # Count total positives and negatives
  total_positives <- sum(target == "yes")
  total_negatives <- sum(target == "no")

  # Calculate confusion matrix values
  TP <- sensitivity * total_positives
  FN <- total_positives - TP
  TN <- specificity * total_negatives
  FP <- total_negatives - TN

  return(list(TP = TP, FN = FN, TN = TN, FP = FP))
}
