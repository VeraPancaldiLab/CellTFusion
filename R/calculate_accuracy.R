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
