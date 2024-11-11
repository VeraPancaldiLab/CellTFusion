#' Calculates AUC from precision-recall curve
#'
#' @param recall A numeric column with the recall values obtained from get_sensitivity_specificity()
#' @param precision A numeric column with the precision values obtained from get_sensitivity_specificity()
#'
#' @return A numeric value corresponding to the AUC from the precision-recall curve
#' @export
#'
#' @examples
#'
#' sens_spec = get_sensitivity_specificity(predict, target, model)
#' auprc = calculate_auprc(sens_spec$recall, sens_spec$precision)
#'
calculate_auprc <- function(recall, precision) {
  # Sort by Recall to ensure trapezoidal rule is correctly applied
  ordered <- order(recall)
  recall <- recall[ordered]
  precision <- precision[ordered]

  auprc <- 0
  for (i in 1:(length(recall) - 1)) { # -1 to avoid NA from the last terms
    # Trapezoidal rule: (precision[i] + precision[i+1]) / 2 * (recall[i+1] - recall[i])
    auprc <- auprc + ((precision[i+1] + precision[i]) / 2) * (recall[i+1] - recall[i])
  }
  return(auprc)
}
