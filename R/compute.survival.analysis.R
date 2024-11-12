#' Compute survival analysis
#'
#' Computes univariate and multivariate cox proportional hazards (coxPH) models and Kaplan-Meier curves to evaluate across cell type groups. After ﬁtting the CoxPH models to different cell type groups combinations, patient are stratiﬁed based on the linear predictors of the model (risk scores) from which we deﬁne as ‘high’ the patients with risk scores above the median value of the cox model’s linear predictors and as ‘low’ the patients below it. It then performes a Kaplan Meier analysis and plotted the survival curves for each risk group. Finally both survival curves are assessed via a log rank test.
#'
#' @param features A matrix with the cell groups scores
#' @param survival.data Survival information containing progression free survival information (named as 'PFS') and the binary value indicating death/progression/recurrence (named as 'DRP_st')
#' @param time_unit Time unit for variables (either days/months/years)
#' @param p.value pvalue for significance. Default is 0.05
#' @param thres Threshold for classifying the risk scores. Default is 0.5 (median).
#' @param max_factors Maximum number of covariables in the lineal model to consider. Default is Inf, meaning it will test all possible combinations and keep only the significant ones.
#'
#' @return A list with the significant combinations formulas along with their survival curves saved in the Results/ directory.
#'
#' @export
#'
#' @examples
#'
#' survival_groups = compute.survival.analysis(features = cell_groups, survival_data, time_unit = "days", p.value = 0.01, max_factors = 3)
#'
compute.survival.analysis = function(features, survival.data, time_unit, p.value = 0.05, thres = 0.5, max_factors = Inf) {
  n_features <- ncol(features)
  significant_combinations <- list() # To store significant feature combinations

  # Generate all possible combinations of the features
  contador = 1
  for (n in 1:min(n_features, max_factors)) {
    combinations <- combn(1:n_features, n, simplify = FALSE)

    for (comb in combinations) {
      # Create a formula dynamically based on the combination
      formula <- as.formula(paste("Surv(time, status) ~", paste(colnames(features)[comb], collapse = " + ")))

      # Prepare the data for survival analysis
      data_for_model <- data.frame("time" = survival.data$PFS,
                                   "status" = survival.data$DRP_st)

      data_for_model = cbind(data_for_model, features[,comb, drop=F])

      # Fit the Cox PH model with the combination of features (cox PH take into account covariates and measure the impact of each variable in the survival time)
      cox <- cph(formula, data = data_for_model)
      data_for_model$CoxPredictors <- cox$linear.predictors #linear predictors is the risk score for each individual in the dataset

      # Check that the model is significant as a predictor (maybe not useful?, it gives the same linear.predictos - to be check)
      cphmodel <- coxph(Surv(time, status) ~ CoxPredictors, data = data_for_model)
      data_for_model$CoxPredictors <- cphmodel$linear.predictors

      quantiles <- quantile(data_for_model$CoxPredictors, thres)

      # Binarize the Cox model output to draw two KM lines (linear predictors are used to stratify between high-risk and low-risk groups)
      data_for_model$coxHL <- ifelse(cphmodel$linear.predictors >= quantiles, 'High', "Low")

      # Perform Kaplan-Meier based on coxHL
      km_fit <- survfit(Surv(time, status) ~ coxHL, data = data_for_model)

      pval <- surv_pvalue(km_fit, data = data_for_model)$pval #Performs log-rank test to see whether both survival curves are significantly different

      if (!is.na(pval) && pval < p.value) {
        significant_combinations[[contador]] <- formula
        names(significant_combinations)[contador] = paste0("Formula_", contador)

        pdf(paste0("Results/SurvPlot_", names(significant_combinations)[contador]), width = 10, height = 5, onefile = FALSE)
        print(ggsurvplot(km_fit,
                         data = data_for_model,
                         size = 1,
                         palette = c("#E7B800", "#2E9FDF"),
                         conf.int.style = "step",
                         pval = TRUE,
                         risk.table = TRUE,
                         risk.table.col = "strata",
                         legend.labs = c("High", "Low"),
                         risk.table.height = 0.3,
                         ggtheme = theme_grey(),
                         title = paste0("Cox PH for ", names(significant_combinations)[contador]),
                         xlab = paste0("Time to death/recurrence/progression (", time_unit, ")")
        ))
        dev.off()
        contador = contador + 1
      }
    }
  }

  if (length(significant_combinations) == 0) {
    print("No significant combinations found.")
  } else {
    return(significant_combinations)
  }

}
