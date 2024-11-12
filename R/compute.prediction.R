#' Compute prediction
#'
#' @param model Trained machine learning model
#' @param test_data A matrix with the testing dataset
#' @param target Character vector with the true values
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
#' prediction = compute.prediction(model, testing_set, target)
#'
compute.prediction = function(model, test_data, target){

  cat("Predicting target variable using provided ML model")

  features <- colnames(test_data)
  are_equal = setequal(model[["coefnames"]], features)
  if(are_equal == T){
    #Predict target variable
    predict <- data.frame(predict(model, test_data, type = "prob"))
    #Get metrics
    sens_spec = get_sensitivity_specificity(predict, target, model$method)
    auroc = calculate_auroc(sens_spec$fpr, sens_spec$sensitivity)
    auprc = calculate_auprc(sens_spec$recall, sens_spec$precision)

    return(list(Metrics = sens_spec, AUC = list("AUROC" = auroc, "AUPRC" = auprc), Predictions = predict))
  }else{
    message("Testing set does not count with the same features as model")
  }
}
