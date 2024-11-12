#' Compute prediction using stacking approach
#'
#' @param super.learner Meta-learner model
#' @param test_data A matrix with the testing dataset
#' @param target Character vector with the true values
#' @param ml.models Machine learning models
#' @param base.models A character vector with the base models for the meta-learner
#'
#' @return A list containing
#'
#' - Prediction metrics
#' - AUROC and AUPRC
#' - Prediction values
#'
#' @export
#'
#' @examples
#'
#' prediction = compute.prediction.stacked(model, testing_set, target, ML_models, base_models)
#'
compute.prediction.stacked = function(super.learner, test_data, target, ml.models, base.models){

  #Learning from simple meta-learner
  base_predictions = list()
  for (i in 1:length(base.models)) {
    base_predictions[[i]] = predict(ml.models[[base.models[i]]], test_data, type = "prob")$yes
    names(base_predictions)[i] = base.models[i]
  }

  base_predictions = do.call(cbind, base_predictions)

  prediction_simple = data.frame(predict(super.learner, base_predictions, type = "prob"))

  # #Learning from simple meta-learner
  # all_predictions = list()
  # for (i in 1:length(ml.models)) {
  #   all_predictions[[i]] = predict(ml.models[[i]], test_data, type = "prob")$yes
  #   names(all_predictions)[i] = names(ml.models)[i]
  # }
  #
  # all_predictions = do.call(cbind, all_predictions)
  #
  # prediction_all = data.frame(predict(super.learner[["all"]], all_predictions, type = "prob"))
  #
  #Metrics

  #Meta-learner simple
  sens_spec_simple = get_sensitivity_specificity(prediction_simple, target, "Meta-learner_simple")
  auroc_simple = calculate_auroc(sens_spec_simple$fpr, sens_spec_simple$sensitivity)
  auprc_simple = calculate_auprc(sens_spec_simple$recall, sens_spec_simple$precision)
  #Meta-learner all
  # sens_spec_all = get_sensitivity_specificity(prediction_all, target, "Meta-learner_all")
  # auc_all = calculate_auc(sens_spec_all$fpr, sens_spec_all$sensitivity)

  #Not returning all (discarded)

  return(list(Metrics = sens_spec_simple, AUC = list("AUROC" = auroc_simple, "AUPRC" = auprc_simple), Predictions = prediction_simple))

}
