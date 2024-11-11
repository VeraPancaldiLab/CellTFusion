#' Calculates accuracy values from prediction
#'
#' @param metrics A Dataframe with metrics obtained using get_sensitivity_specificity() function
#' @param target A character vector containing the true values from the target variable
#'
#' @return A numeric vector with the accuracy values
#' @export
#'
#' @examples
#'
#'observed_values = c("yes", "yes", "no", "yes")
#'accuracy = calculate_accuracy(metrics, observed_values)
#'
calculate_accuracy <- function(metrics, target) {
  sensitivity = metrics[,"sensitivity"]
  specificity = metrics[,"specificity"]
  total_positives = sum(target == "yes")
  total_negatives = sum(target == "no")
  TP <- sensitivity * total_positives
  FN <- total_positives - TP
  TN <- specificity * total_negatives
  FP <- total_negatives - TN

  accuracy <- (TP + TN) / (TP + TN + FP + FN)

  return(accuracy)
}
