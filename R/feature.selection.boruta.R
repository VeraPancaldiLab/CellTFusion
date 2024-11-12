#' Compute feature selection using repeated Boruta algorithm
#'
#' @param data A dataframe with the target column name as "target" and the features to test.
#' @param iterations Number of iterations to run Boruta.
#' @param fix Parameter from Boruta(). It tests whether the features classified as 'Tentative' need to be judge with an additional test to either confirm then or not (See TentativeRoughFix() from Boruta package for more information)
#' @param doParallel Whether to do or not parallelization.
#' @param workers Number of processes available to run on parallel. If no number is set, this will correspond to detectCores() - 1
#' @param file_name File name for the csv files and plots saved in the Results/ directory
#' @param threshold Threshold to consider the features as confirmed after several iterations of Boruta. If boruta_threshold = 0.8, features labeled as confirmed in more than 80% of the times will be finally considered as confirmed.
#' @param return Whether to save or not the plots in the Results/ directory.
#'
#' @return A list containing
#'
#' - Confirmed features
#' - Tentative features
#' - Matrix of feature importance
#'
#' @export
#'
#' @examples
#'
#' res_boruta = feature.selection.boruta(training_set, iterations = 100, fix = F, doParallel = T, workers = 4, threshold = 0.8, file_name = "Test", return = T)
#'
feature.selection.boruta <- function(data, iterations = NULL, fix, doParallel = F, workers=NULL, file_name = NULL, threshold = NULL, return) {
  if(doParallel){
    if(is.null(iterations) == T){
      stop("No iterations specified for running in parallel, please set a number. If you want to run feature selection once consider setting doParallel = F")
    }else{
      if(is.null(workers)==T){
        num_cores <- detectCores() - 1
      }else{
        num_cores <- workers
      }

      cl = parallel::makeCluster(num_cores) #Forking just copy the R session in its current state. - makeCluster() all must be exported (copied) to the clusters, which can add some overhead
      doParallel::registerDoParallel(cl)

      message("Running ", iterations, " iterations of the Boruta algorithm using ", num_cores, " cores")
      # arg_list <- replicate(iterations, list(data, sample.int(100000, 1), fix), simplify = FALSE)
      # system.time({
      #   res <- mclapply(arg_list, function(x) {
      #     do.call(compute.boruta, x)
      #   }, mc.cores = num_cores)
      # })

      res <- foreach(seed = sample.int(100000, iterations)) %dopar% {

        source("src/environment_set.R")

        tryCatch({
          # If successful, return the result and the seed
          list(result = compute.boruta(data, seed, fix),
               error = NULL,
               seed = seed)
        }, error = function(e) {
          # If an error occurs, return the error message and the seed for debugging
          list(result = NULL, error = e$message, seed = seed)
        })
      }

      parallel::stopCluster(cl)
      unregister_dopar() #Stop Dopar from running in the background

    }

    # Extract the first sublist of each element
    matrix_of_importance <- lapply(res, function(x) x[[1]])

    # Extract the second sublist of each element
    features_labels <- lapply(res, function(x) x[[2]])

    res = merge_boruta_results(matrix_of_importance, features_labels, file_name = file_name, iterations = iterations, threshold = threshold, return = return)
  }else{
    if(is.null(iterations) == T){
      stop("No iterations specified for running Boruta algorithm for feature selection")
    }else{
      message("Running ", iterations, " iterations of the Boruta algorithm")
      res = list()
      for (i in 1:iterations) {
        res[[i]] = compute.boruta(data, seed = sample.int(100000, 1), fix)
      }

      # Extract the first sublist of each element
      matrix_of_importance <- lapply(res, function(x) x[[1]])

      # Extract the second sublist of each element
      features_labels <- lapply(res, function(x) x[[2]])

      res = merge_boruta_results(matrix_of_importance, features_labels, file_name = file_name, iterations = iterations, threshold = threshold, return = return)
    }

  }

  return(res)

}
