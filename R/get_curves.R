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
