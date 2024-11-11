#' Calculate precision values
#'
#' @param metrics A Dataframe with metrics obtained using get_sensitivity_specificity()
#' @param target A character vector containing the true values from the target variable
#'
#' @return A numeric vector with precision values
#' @export
#'
#' @examples
#'
#' observed = c("yes", "no", "yes", "no", "no")
#' metrics = get_sensitivity_specificity(predict, target, model)
#' precision = calculate_precision(metrics, observed)
#'
calculate_precision <- function(metrics, target) {
  confusion_values <- calculate_confusion_values(metrics, target)
  TP <- confusion_values$TP
  FP <- confusion_values$FP

  precision <- TP / (TP + FP)

  return(precision)
}
