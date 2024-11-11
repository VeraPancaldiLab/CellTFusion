compute.bootstrap.ML = function(raw.counts, normalized = F, clinical, trait, trait.positive, partition = 0.8, metric = "Accuracy", iterations, feature.selection = F, stack, deconv_methods = c("Quantiseq", "MCP", "xCell", "CBSX", "Epidish", "DeconRNASeq", "DWLS"), workers = NULL, file.name = NULL, return = F){

  dir.create(file.path(getwd(), "Results/ML_models"))

  if(is.null(iterations) == T){
    stop("No iterations specified, please set a number")
  }else{
    if(is.null(workers)==T){
      num_cores <- detectCores() - 1
    }else{
      num_cores <- workers
    }

    #options(cluster.timeout = 300)
    cl = parallel::makeCluster(num_cores) #Forking just copy the R session in its current state. - makeCluster() all must be exported (copied) to the clusters, which can add some overhead
    doParallel::registerDoParallel(cl)

    #future::plan("multicore", workers = num_cores)
    #future::plan()

    message("Running ", iterations, " splits for training and test using ", num_cores, " cores")

    #List of arguments inputs
    #arg_list <- replicate(iterations, list(norm.counts, clinical, trait, trait.positive, deconvolution, tfs.matrix, partition, sample.int(100000, 1)), simplify = FALSE)


    # Run in parallel with error handling
    # system.time({
    #   res <- mclapply(seq_along(arg_list), function(i) {
    #     args <- arg_list[[i]]
    #     result = tryCatch({
    #         list(result = do.call(compute.ML, args), error = NULL, data = args)
    #       }, error = function(e) {
    #         list(result = NULL, error = e$message, data = args)
    #       })
    #     rm(args) #Remove input arguments to free memory after each iteration
    #     invisible(gc()) #Do garbage collection
    #     return(result)
    #   }, mc.cores = num_cores)
    # })

    #options(future.globals.maxSize = 10 * 1024^3)  # Set limit to 10 GiB
    # res <- future.apply::future_lapply(sample.int(100000, iterations), function(seed) {
    #   compute.ML(norm.counts, clinical, trait, trait.positive, deconvolution, tfs.matrix, partition, seed, universe)
    # }, future.seed = TRUE)

    # res <- future.apply::future_lapply(sample.int(100000, iterations), function(seed) {
    #   # Use tryCatch to handle any errors during the function call
    #   result <- tryCatch({
    #     # If successful, return the result and the seed
    #     list(result = compute.ML(norm.counts, clinical, trait, trait.positive, deconvolution, tfs.matrix, partition, seed, universe),
    #          error = NULL,
    #          seed = seed)
    #   }, error = function(e) {
    #     # If an error occurs, return the error message and the seed for debugging
    #     list(result = NULL, error = e$message, seed = seed)
    #   })
    #
    #   return(result)
    # }, future.seed = TRUE)


    # Run foreach loop using each random seed directly
    foreach(iteration = seq_len(iterations), random.seed = sample.int(100000, iterations)) %dopar% {

      # Use absolute path for the source file to avoid path issues
      source("src/environment_set.R")

      tryCatch({
        # Compute the result with the current random seed
        result <- compute.ML(raw.counts, normalized, clinical, trait, trait.positive, partition, metric, stack,
                             feature.selection, doParallel = F, workers = NULL, seed = random.seed, deconv_methods = deconv_methods,
                             file_name = file.name, return = return)

        # Save result as RDS file with unique identifier based on iteration (random seed)
        saveRDS(list(result = result, seed = random.seed),
                file = file.path(paste0("Results/ML_models/ML_result_", iteration, ".rds")))

      }, error = function(e) {
        # Save error information as RDS file with random seed identifier
        saveRDS(list(result = NULL, error = e$message, seed = random.seed),
                file = file.path(paste0("Results/ML_models/ML_result_", iteration, ".rds")))
      })


    }

    parallel::stopCluster(cl)
    unregister_dopar() #Stop Dopar from running in the background

  }


  # # Extract the first sublist of each element
  # matrix_of_importance <- lapply(res, function(x) x[[1]])
  #
  # # Extract the second sublist of each element
  # features_labels <- lapply(res, function(x) x[[2]])
  #
  # res = merge_boruta_results(matrix_of_importance, features_labels, file_name, iterations, threshold = 0.7)

  #future::plan(future::sequential)

  ###################################################


  # #########ROC curve with a specific list of thresholds
  # x <- seq(0, 1, by = 0.01) #List of thresholds
  #
  # result <- lapply(res, function(mat) {
  #   # For each matrix, apply the thresholding function over x
  #   lapply(x, function(x_val) {
  #     mat[["result"]][["Prediction_metrics"]] %>%
  #       select(sensitivity, specificity, fpr, model) %>%
  #       filter(specificity - x_val >= 0) %>% #Define specific set of thresholds of specificities
  #       top_n(sensitivity, n = 1) %>% #Extract maximum value of sensitivity
  #       mutate(specificity = x_val, fpr = 1 - x_val) %>% #Add fpr value
  #       distinct()
  #   }) %>%
  #     bind_rows() # Combine results for each value of specificity into a single data frame for this matrix
  # }) %>%
  #   bind_rows() # Combine matrices from all iterations
  #
  # roc = result %>%
  #   group_by(model, specificity) %>%
  #   summarise(lquartile = quantile(sensitivity, prob = 0.25),
  #             uquartile = quantile(sensitivity, prob = 0.75),
  #             sensitivity = median(sensitivity),
  #             .groups = "drop") %>%
  #   mutate(fpr = 1 - specificity) %>%
  #   ggplot(aes(x=1-specificity, y = sensitivity,ymin = lquartile, ymax = uquartile)) +
  #   geom_ribbon(alpha = 0.25, aes(fill=model)) +
  #   geom_step(aes(color=model)) +
  #   theme_classic()
  #
  # pdf(paste0("Results/ROC_curve_iterations_", file.name, ".pdf"))
  # print(roc)
  # dev.off()
  #
  message("Analysis is done!")

  message("ML models are saved in Results/ML_models folder")

}
