#'
#' Get performance curves
#'
#' Get ROC and precision-recall curves
#'
#' @param data A matrix with the prediction metrics
#' @param spec Name of column with the specificity values
#' @param sens Name of column with the sensitivity values
#' @param reca Name of column with the recall values
#' @param prec Name of column with the precision values
#' @param color Name of column with the cohort names. Each cohort will have a color. If several cohorts are present, different curves will be plot.
#' @param auc_roc AUC-ROC value
#' @param auc_prc AUC-PRC value
#' @param plot_title Title for the plots
#'
#' @return ROC and precision-recall curves saved in Results/ directory.
#' @export
#'
#' @examples
#'
#' get_curves(metrics, "specificity", "sensitivity", "recall", "precision", "model", auc_roc_score, auc_prc_score, "Test")
#'
get_curves = function(data, spec, sens, reca, prec, color, auc_roc, auc_prc, plot_title){

  data = data %>%
    mutate(specificity = data[,spec],
           sensitivity = data[,sens],
           recall = data[,reca],
           precision = data[,prec],
           color = data[,color])

  #Add AUC scores to data frame
  data$color.roc <- paste(data$color, "\n(AUC-ROC =", round(auc_roc, 2), ")\n")

  # Plot the ROC curves
  roc = ggplot(data = data, aes(x = 1- specificity, y = sensitivity, color = color.roc)) +
    geom_line() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
    labs(title = "ROC Curve", subtitle = plot_title, x = "1 - Specificity", y = "Sensitivity") +
    theme_minimal() +
    theme(legend.title = element_blank())

  #Add AUC scores to data frame
  data$color.prc <- paste(data$color, "\n(AUC-PRC =", round(auc_prc, 2), ")\n")

  # Plot recall curves
  recall = ggplot(data = data, aes(x = recall, y = precision, color = color.prc)) +
    geom_line() +
    labs(title = "Precision-Recall Curve", subtitle = plot_title, x = "Recall", y = "Precision") +
    ylim(0, 1) +
    theme_minimal() +
    theme(legend.title = element_blank())

  pdf(paste0("Results/ROC_curve_", file.name, ".pdf"))
  print(roc)
  dev.off()

  pdf(paste0("Results/Recall_curve_", file.name, ".pdf"))
  print(recall)
  dev.off()

}
