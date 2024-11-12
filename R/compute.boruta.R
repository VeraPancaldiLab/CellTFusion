#' Compute boruta algorithm
#'
#' @param data A dataframe with the target column name as "target" and the features to test
#' @param seed A numeric value to set the seed for ensuring reproducibility
#' @param fix Parameter from Boruta(). It tests whether the features classified as 'Tentative' need to be judge with an additional test to either confirm then or not (See TentativeRoughFix() from Boruta package for more information)
#'
#' @return A list containing:
#'
#' - Matrix of feature importance
#' - A character vector with the decision corresponding to each feature
#'
#' @export
#'
#' @examples
#'
#' res = compute.boruta(data, seed = 123, fix = F)
#'
compute.boruta <- function(data, seed, fix = TRUE) {

  set.seed(seed)
  boruta_output <- Boruta(target ~ ., data = data, doTrace = 0)

  if (fix) {
    roughFixMod <- TentativeRoughFix(boruta_output)
    boruta_output <- roughFixMod
  }

  imps <- attStats(boruta_output)
  decision <- as.character(imps$decision)

  res <- imps %>%
    data.frame() %>%
    rownames_to_column("Variable") %>%
    dplyr::select(-decision)


  return(list(res, decision))
}
