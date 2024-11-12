#' Merge boruta results
#'
#' Merge boruta results after performing several iterations of the algorithm.
#'
#' @param importance_values A matrix with the importance values along each iteration.
#' @param decisions A matrix with the decision labels along each iteration.
#' @param file_name File name for plots.
#' @param iterations Number of iterations performed.
#' @param threshold A numeric value with the threshold for considering final labels (e.g. features labeled as 'confirmed' more than > threshold are final considered as 'confirmed')
#' @param return Whether to save or not the plots in Results/ directory
#'
#' @return A list containing:
#'
#' - Confirmed features
#' - Tentative features
#' - Matrix of feature importance
#'
#' @export
#'
#' @examples
#'
#' res = merge_boruta_results(matrix_of_importance, features_labels, file_name = "Test", iterations = 100, threshold = 0.8, return = T)
#'
merge_boruta_results = function(importance_values, decisions, file_name, iterations, threshold, return = T){

  ### Construct matrix of importance
  combined_importance <- do.call(rbind, importance_values)
  combined_results_long <- combined_importance %>% #Matrix for plotting
    pivot_longer(cols = meanImp, names_to = "Measure", values_to = "Value")

  median_df <- combined_importance %>% #Calculate the median for each column, grouped by the variable name
    group_by(Variable) %>%
    dplyr::summarize(across(everything(), \(x) median(x, na.rm = TRUE)))

  ### Retrieve important and tentatives variables

  combined_results <- do.call(cbind, decisions)
  rownames(combined_results) = median_df$Variable
  decisions_summary <- apply(combined_results, 1, function(x) {
    table(factor(x, levels = c("Confirmed", "Tentative", "Rejected")))
  })
  confirmed_vars <- names(which(decisions_summary["Confirmed",] >= round(threshold*iterations)))
  tentative_vars <- names(which(decisions_summary["Tentative",] >= round(threshold*iterations)))

  # For plotting
  combined_results_long$Decision = "Rejected"
  combined_results_long$Decision[which(combined_results_long$Variable %in% confirmed_vars)] = "Confirmed"
  combined_results_long$Decision[which(combined_results_long$Variable %in% tentative_vars)] = "Tentative"

  mean_order <- median_df %>% #Extract the order of variables for plotting
    arrange(meanImp) %>%
    pull(Variable)

  # For result
  median_df$Decision = "Rejected"
  median_df$Decision[which(median_df$Variable %in% confirmed_vars)] = "Confirmed"
  median_df$Decision[which(median_df$Variable %in% tentative_vars)] = "Tentative"

  # Plot variable importance boxplots
  if(return){
    pdf(paste0("Results/Boruta_variable_importance_", file_name, ".pdf"), width = 8, height = 12)
    print(ggplot(combined_results_long, aes(x = factor(Variable, levels = mean_order), y = Value, fill = Decision)) +
            geom_bar(stat = "identity", position = "dodge") +
            theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
            coord_flip() +
            labs(x = "Features", y = "Importance", title = paste0("Variable Importance by Boruta after ", iterations, " bootstraps\n", file_name)) +
            scale_fill_manual(values = c("Confirmed" = "green", "Tentative" = "yellow", "Rejected" = "red")) +
            facet_wrap(~ Measure, scales = "free_y"))
    dev.off()
  }

  return(list(Confirmed = confirmed_vars, Tentative = tentative_vars, Matrix_Importance = median_df))
}
