get_pooled_roc_curves = function(file.name){

  # Get a list of all RDS files in the folder
  folder_path <- "Results/ML_models"
  res <- list.files(folder_path, pattern = "\\.rds$", full.names = TRUE)
  res <- lapply(res, readRDS) #Read each RDS file
  iterations = length(res)

  #########Boxplot
  auc_roc_list <- lapply(res, function(sublist) {
    sublist[["result"]][["AUC"]][["AUROC"]]
  })

  auc_prc_list <- lapply(res, function(sublist) {
    sublist[["result"]][["AUC"]][["AUPRC"]]
  })

  auc_roc_list = as.numeric(unlist(auc_roc_list))
  auc_prc_list = as.numeric(unlist(auc_prc_list))

  data = data.frame("AUC_roc" = auc_roc_list,
                    "AUC_prc" = auc_prc_list,
                    "Cohort" = file.name)

  mean_auc_roc = data %>%
    group_by(Cohort) %>%
    dplyr::summarize(meanAUROC = mean(AUC_roc))

  mean_auc_prc = data %>%
    group_by(Cohort) %>%
    dplyr::summarize(meanAUPRC = mean(AUC_prc))

  # Plot boxplot with mean AUC annotations
  plot_roc = ggplot(data, aes(x = Cohort, y = AUC_roc, fill = Cohort)) +
    geom_boxplot() +
    labs(title = paste0("Distribution of AUROC values across ", iterations, " splits"),
         x = "Model",
         y = "AUROC") +
    theme_minimal() +
    theme(legend.position = "right") +
    geom_text(data = mean_auc_roc, aes(x = Cohort, y = max(data$AUC_roc),
                                       label = paste("Mean AUC:", round(meanAUROC, 3))),
              size = 4, color = "black", vjust = -0.5)

  pdf(paste0("Results/Boxplot_AUPRC_performance_", file.name, ".pdf"))
  print(plot_roc)
  dev.off()

  plot_prc = ggplot(data, aes(x = Cohort, y = AUC_prc, fill = Cohort)) +
    geom_boxplot() +
    labs(title = paste0("Distribution of AUPRC values across ", iterations, " splits"),
         x = "Model",
         y = "AUPRC") +
    theme_minimal() +
    theme(legend.position = "right") +
    geom_text(data = mean_auc_prc, aes(x = Cohort, y = max(data$AUC_prc),
                                       label = paste("Mean AUC:", round(meanAUPRC, 3))),
              size = 4, color = "black", vjust = -0.5)

  pdf(paste0("Results/Boxplot_AUROC_performance_", file.name, ".pdf"))
  print(plot_prc)
  dev.off()



}
