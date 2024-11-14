#' Compute CellTFusion machine learning pipeline
#'
#' @param raw.counts A matrix with the raw counts with samples as columns and genes symbols as rows.
#' @param normalized Boolean value to specify if raw_counts need to be normalized (If no raw_counts are available and argument corresponds to already normalized counts this arguments needs to be set to False)
#' @param clinical A data frame with the clinical data containing the target variable to analyze
#' @param trait Column name of target variable
#' @param trait.positive Value of target consider as positive for prediction
#' @param partition Proportion of data going to train
#' @param metric Metric to based the methods choice during cross validation (CV). Either accuracy or AUC
#' @param stack If stacking needs to be applied
#' @param feature.selection If feature selection using Boruta algorithm needs to be computed
#' @param deconv_methods A character vector with the deconvolution methods to run. Default are "Quantiseq", "MCP", "xCell", "CBSX", "Epidish", "DeconRNASeq", "DWLS"
#' @param doParallel Whether to do or not parallelization. Only CBSX and DWLS methods will run in parallel.
#' @param workers Number of processes available to run on parallel. If no number is set, this will correspond to detectCores() - 1
#' @param seed Seed to ensure reproducibility regarding the train/test split
#' @param file_name File name for the csv files and plots saved in the Results/ directory
#' @param return Return the plots printed during each iteration
#'
#' @return A list containing
#'
#' - Trained model
#' - Features used to train the model
#' - Feature importance matrix
#' - Cell groups
#' - Prediction metrics
#' - AUROC and AUPRC
#' - Predictions
#'
#' @export
#'
#' @examples
#'
#' res_ml = compute.ML(raw.counts, normalized = T, clinical, trait = "Response",trait.positive = "CR", partition = 0.8, metric = "AUC", stack = T, feature.selection = F,seed = 1234, doParallel = T,  workers = 2, file_name = "Test", return = T)
#'
compute.ML = function(raw.counts, normalized = F, clinical, trait, trait.positive, partition, metric = "Accuracy", stack, feature.selection = F, deconv_methods = c("Quantiseq", "MCP", "xCell", "CBSX", "Epidish", "DeconRNASeq", "DWLS"), doParallel = F, workers = NULL, seed, file_name = NULL, return = F){
  set.seed(seed)

  # Do stratified partition
  index = createDataPartition(clinical[,trait], times = 1, p = partition, list = FALSE)

  # Normalize counts
  if(normalized == T){
    norm.counts = data.frame(ADImpute::NormalizeTPM(raw.counts, log = T))
  }else{
    norm.counts = raw.counts
  }

  # Train cohort
  traitData_train = clinical[index, ]
  raw.counts_train = raw.counts[,index]
  counts.normalized_train = norm.counts[,index]
  tfs_train = compute.TFs.activity(counts.normalized_train)
  deconv_train = compute.deconvolution(raw.counts_train, normalized = normalized, methods = deconv_methods, doParallel = doParallel, workers = workers, credentials.mail = "marcelo.hurtado@inserm.fr", credentials.token = "734212f6ad77fc4eea2bdb502792f294", return = return)

  # Test cohort
  traitData_test = clinical[-index, ]
  raw.counts_test = raw.counts[,-index]
  counts.normalized_test = norm.counts[,-index]
  tfs_test = compute.TFs.activity(counts.normalized_test)
  deconv_test = compute.deconvolution(raw.counts_test, normalized = normalized, methods = deconv_methods, doParallel = doParallel, workers = workers, credentials.mail = "marcelo.hurtado@inserm.fr", credentials.token = "734212f6ad77fc4eea2bdb502792f294", return = return)

  ###############################################################################################################################################################################

  # CellTFusion
  set.seed(seed)
  network = compute.WTCNA(tfs_train, corr_mod = 0.9, clustering.method = "ward.D2", return = return)
  pathways = compute.pathway.activity(counts.normalized_train, paths = NULL)
  tfs.modules.clusters = compute.TF.network.classification(network, pathways, return = return)
  dt = compute.deconvolution.analysis(deconv_train, corr = 0.7, seed = seed, return = return)
  corr_modules = compute.modules.relationship(network[[1]], dt[[1]], return = T, plot = return)
  cell_dendrograms = identify.cell.groups(corr_modules, tfs.modules.clusters, height = 20, return = return)
  cell.groups = cell.groups.computation(dt[[1]], tfs.module.network = network, cell.dendrograms = cell_dendrograms, return = return) #Identify cell groups with specific cut

  ###############################################################################################################################################################################

  ####################################################Training

  #Set training set
  train_data = cell.groups[[1]] %>%
    data.frame() %>%
    mutate(Trait = traitData_train[,trait],
           target = as.factor(ifelse(Trait == trait.positive, 'yes', 'no'))) %>%
    dplyr::select(-Trait)

  train_data$target <- factor(train_data$target, levels = c("no", "yes"))  # Order, just in case to ensure positive class is not well defined

  #Cross-validation training (5 k-folds and 100 repetitions)
  training = compute.k_fold_CV(train_data, k_folds = 5, n_rep = 100, metric = metric, stacking = stack, boruta = feature.selection, boruta_iterations = 100, fix_boruta = F, boruta_threshold = 0.8, file_name = file_name, return= return)

  ####################################################Predicting
  if(length(training)!=0){
    cell_groups = cell.groups #Save cell groups scores per partition
    features = training[["Features"]] #Save selected features per partition
    ####################### Testing set
    testing_set = compute_cell_groups_signatures(dt, network, cell.groups, features, deconv_test, tfs_test) #Cell groups projection
    #Extract target variable
    target = traitData_test %>%
      mutate(target = ifelse(traitData_test[,trait] == trait.positive, "yes", "no")) %>%
      pull(target)

    target = factor(target, levels = c("no", "yes"))

    if(stack){
      model = training[["Meta_learner"]]
      var_importance = calculate_feature_importance_stacking(training[["Variable_importance"]], training[["Base_models"]], model)
      prediction = compute.prediction.stacked(model, testing_set, target, training[["ML_models"]], training[["Base_models"]])

    }else{
      model = training[["Model"]] #Save best ML model based on the Accuracy/AUC from CV per partition
      var_importance = varImp(model, scale = F) #Retrieve variable importance
      prediction = compute.prediction(model, testing_set, target)
    }

    auc_roc_score = prediction[["AUC"]][["AUROC"]]
    auc_prc_score = prediction[["AUC"]][["AUPRC"]]

    metrics = prediction[["Metrics"]]
    predictions = prediction[["Predictions"]]

    if(return == T){
      get_curves(metrics, "specificity", "sensitivity", "recall", "precision", "model", auc_roc_score, auc_prc_score, file_name)
    }

    rm(network, pathways, tfs.modules.clusters, dt, corr_modules, cell_dendrograms, cell.groups,
       traitData_train, counts.normalized_train, tfs_train, deconv_train,
       traitData_test, counts.normalized_test, deconv_test, tfs_test,
       clinical, norm.counts) #Remove variables

    gc() #Clean garbage

    return(list(Model = model, Features = features, Variable_importance = var_importance, Cell_groups = cell_groups, Prediction_metrics = metrics, AUC = list(AUROC = auc_roc_score, AUPRC = auc_prc_score), Prediction = predictions))
  }else{  #No features are selected as predictive

    rm(network, pathways, tfs.modules.clusters, dt, corr_modules, cell_dendrograms, cell.groups,
       traitData_train, counts.normalized_train, tfs_train, deconv_train,
       traitData_test, counts.normalized_test, deconv_test, tfs_test,
       clinical, norm.counts) #Remove variables

    gc() #Clean garbage

    message("No features selected as predictive after Boruta runs. No model returned.")

    return(NULL)
  }

}

#' Compute machine learning models across different splits of train and test
#'
#' @param raw.counts A matrix with raw counts (genes symbols as rows and samples as columns)
#' @param normalized A boolean value to specify if raw counts have to be normalized (this should always be true unless the provided raw.counts are already the normalized ones)
#' @param clinical A data frame with the clinical data containing the target variable to analyze
#' @param trait Column name of target variable
#' @param trait.positive Value of target consider as positive for prediction
#' @param partition Proportion of data going to train
#' @param metric Metric to based the methods choice during cross validation (CV). Either accuracy or AUC
#' @param iterations Number of iterations to perform
#' @param feature.selection If feature selection using Boruta algorithm needs to be computed
#' @param stack If stacking needs to be applied
#' @param deconv_methods Deconvolution methods to run during CellTFusion
#' @param workers Number of workers to use for parallelization
#' @param file.name File name for plots
#' @param return Return the plots printed during each iteration
#'
#' @return Machine learning models are saved in Results/ folder as .rds files
#' @export
#'
#' @examples
#'
#' compute.bootstrap.ML(raw.counts, normalized = T, clinical, trait = "Response", trait.positive = "YES", partition = 0.8, metric = "Accuracy", iterations = 20, feature.selection = F, stack = T, workers = 4, file.name = "Test", return = F)
#'
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

#' Compute k fold cross validation
#'
#' Perfoms a repeated and stratified k fold cross validation on a dataset in order to train and tune hyperparameters of 13 machine learning methods based either on the Accuracy or AUC scores from their predictions.
#'
#' @param model A dataframe with the features and a target column named 'target' corresponding to the trait to predict
#' @param k_folds Number of k folds to perform during cross validation (Default is 5)
#' @param n_rep Number of repeated folds to perform during cross validation (Default is 100)
#' @param stacking Either to do or not stacking
#' @param metric Metric to be used for tuning the hyperparameters and selecting the base models (if stacking = T). Metrics supported are Accuracy and AUC.
#' @param boruta Whether to do or not Boruta for feature selection before training the model. Take into account that many ML models already give a weight to the features so no previous feature selection approach is necessarily (unless there are features with multicollinearity)
#' @param boruta_iterations Number of iterations Boruta needs to be run. Boruta applies random number generator in each run so to be consistent regarding the confirmed features we advised to do several iterations of the algorithm (Default is 100)
#' @param fix_boruta Parameter from Boruta(). See compute.boruta() for more information
#' @param tentative Whether to consider the tentative features as part of the dataframe to train or not.
#' @param boruta_threshold Threshold to consider the features as confirmed after several iterations of Boruta. If boruta_threshold = 0.8, features labeled as confirmed in more than 80% of the times will be finally considered as confirmed.
#' @param file_name File name for the plots to be saved in the Results/ directory.
#' @param return Whether to return and save the plots generated during the function
#'
#' @return A list containing
#' - Features used during training
#' - The selected machine learning model
#' - All the machine learning models trained.
#'
#' If stacking = T, list will also include
#' - Base models
#' - Meta-learner
#' - Matrix of weighted feature importance. See calculate_feature_importance_stacking() for more information.
#'
#' @export
#'
#' @examples
#'
#'  training = compute.k_fold_CV(train_data, k_folds = 5, n_rep = 100, metric = "Accuracy", stacking = T, boruta = F, boruta_iterations = 100, fix_boruta = F, boruta_threshold = 0.8, file_name = "Test", return= T)
#'
compute.k_fold_CV = function(model, k_folds, n_rep, stacking = F, metric = "Accuracy", boruta, boruta_iterations = NULL, fix_boruta = NULL, tentative = F, boruta_threshold = NULL, file_name = NULL, return){

  if(!(metric %in% c("AUC","Accuracy"))){
    stop("The metric assigned is not supported. Choose either accuracy or AUC.")
  }

  ######### Feature selection across folds
  if(boruta == T){
    cat("Feature selection using Boruta...............................................................\n\n")
    # Feature selection using Boruta
    res_boruta = feature.selection.boruta(model, iterations = boruta_iterations, fix = fix_boruta, file_name = file_name, doParallel = F, workers=NULL, threshold = boruta_threshold, return = return)

    if(tentative == F){
      if(length(res_boruta$Confirmed) <= 1){ #No enough features selected for training model
        message("No enough features selected for training a model")
        results = list()
        return(results)
      }else{
        cat("\nKeeping only features confirmed in more than 80% of the times for training...............................................................\n\n")
        cat("If you want to consider also tentative features, please specify tentative = T in the parameters.\n\n")
        train_data = model[,colnames(model)%in%res_boruta$Confirmed, drop = F] %>%
          mutate(target = model$target)
      }
    }else{
      sum_features = length(res_boruta$Confirmed) + length(res_boruta$Tentative)
      if(sum_features <= 1){
        message("No enough features selected for training a model")
        results = list()
        return(results)
      }else{
        cat("Keeping features confirmed and tentative in more than 80% of the times for training...............................................................\n\n")
        train_data = model[,colnames(model)%in%c(res_boruta$Confirmed, res_boruta$Tentative), drop = F] %>%
          mutate(target = model$target)
      }
    }

    rm(res_boruta) #Clean memory
    gc()

  }else{
    train_data = model
  }

  rm(model) #Clean memory
  gc()

  cat("Training machine learning model...............................................................\n\n")

  ######### Machine Learning models

  ######### Stratify K fold cross-validation
  #folds <- createFolds(train_data[,'target'], k = k_folds, returnTrain = T, list = T) #this for single folds
  multifolds <- createMultiFolds(train_data[,'target'], k = k_folds, times = n_rep) #repeated folds
  trainControl <- trainControl(index = multifolds, method="repeatedcv", number=k_folds, repeats=n_rep, verboseIter = F, allowParallel = F, classProbs = TRUE, savePredictions=T)

  ######### Feature selection across folds (DEPRECATED)
  # if(boruta == T){
  #   features = computed.boruta.kfolds(multifolds, model, boruta_iterations = boruta_iterations, fix_boruta = fix_boruta, tentative = tentative, threshold = boruta_threshold, file_name = file_name)
  #   if(is.null(features)==T){
  #     cat("No features selected across folds after Boruta. Try different parameter values")
  #     return(NULL)
  #   }else{
  #     train_data = model[,colnames(model)%in%features, drop = F] %>%
  #       mutate(target = model$target)
  #   }
  # }else{
  #   train_data = model
  # }


  ##################################################### ML models
  #To do: Re-calculate accuracy values based on tuning parameters optimized by the cv AUC - now the values are based on accuracy! be careful

  ################## Bagged CART
  fit.treebag <- train(target~., data = train_data, method = "treebag", metric = "Accuracy",trControl = trainControl)

  ################## RF
  require(randomForest)
  fit.rf <- train(target~., data = train_data, method = "rf", metric = "Accuracy",trControl = trainControl)

  ################## C5.0
  require(C50)
  fit.c50 <- train(target~., data = train_data, method = "C5.0", metric = "Accuracy",trControl = trainControl)

  ################## LG - Logistic Regression
  fit.glm <- train(target~., data = train_data, method="glm", metric="Accuracy",trControl=trainControl)

  ################## LDA - Linear Discriminate Analysis
  fit.lda <- train(target~., data = train_data, method="lda", metric="Accuracy",trControl=trainControl)

  ################## GLMNET - Regularized Logistic Regression (Elastic net)
  fit.glmnet <- train(target~., data = train_data, method="glmnet", metric="Accuracy",trControl=trainControl)

  ################## KNN - k-Nearest Neighbors
  fit.knn <- train(target~., data = train_data, method="knn", metric="Accuracy",trControl=trainControl)

  ################## CART - Classification and Regression Trees (CART),
  fit.cart <- train(target~., data = train_data, method="rpart", metric="Accuracy",trControl=trainControl)

  # NB - Naive Bayes (NB)
  #Grid = expand.grid(usekernel=TRUE,adjust=1,fL=c(0,0.2,0.5,0.8,1))
  #fit.nb <- train(target~., data = train_data, method="nb", metric="Accuracy",trControl=trainControl, tuneGrid=Grid)

  ################## Regularized Lasso
  fit.lasso <- train(target~., data = train_data, method="glmnet", metric="Accuracy",trControl=trainControl, tuneGrid = expand.grid(alpha = 1, lambda = seq(0.001, 1, length = 20)))

  ################## Ridge regression
  fit.ridge <- train(target~., data = train_data, method="glmnet", metric="Accuracy",trControl=trainControl, tuneGrid = expand.grid(alpha = 0, lambda = seq(0.001, 1, length = 20)))

  ################## Support Vector Machine with Radial Kernel
  fit.svm_radial <- train(target ~ ., data = train_data, method = "svmRadial", metric = "Accuracy", trControl = trainControl)

  ################## Support Vector Machine with Linear Kernel
  fit.svm_linear <- train(target ~ ., data = train_data, method = "svmLinear", metric = "Accuracy", trControl = trainControl)

  ####### Optimized based on metric (only AUC or Accuracy available)
  if(metric == "AUC"){

    ################################################Bagged CART

    ## Integrate AUCs into prediction matrix
    fit.treebag$pred = fit.treebag$pred %>%
      group_by(Resample) %>%
      mutate(AUC = calculate_auc_resample(obs, yes)) %>% #Calculate resamples AUC scores
      ungroup()

    ## Integrate AUCs into resamples matrix
    auc = c()
    for (i in 1:nrow(fit.treebag$resample)) {
      auc_val = fit.treebag$pred %>%
        filter(Resample == fit.treebag$resample$Resample[i]) %>%
        pull(AUC) %>% #AUC per resample is the same
        unique()

      auc = c(auc, auc_val)
    }

    fit.treebag$resample = fit.treebag$resample %>%
      mutate(AUC = auc) %>%
      select(AUC, everything())

    ## Integrate average CV AUCs into results
    fit.treebag$results = fit.treebag$results %>%
      mutate(AUC = mean(fit.treebag$resample$AUC))

    ################################################Random Forest

    ## Integrate AUCs into prediction matrix
    fit.rf$pred = fit.rf$pred %>%
      group_by(Resample, mtry) %>% #Parameters for tunning
      mutate(AUC = calculate_auc_resample(obs, yes)) %>%
      ungroup()

    ## Integrate AUCs into results per parameter
    auc_values = fit.rf$pred %>%
      group_by(mtry) %>%
      summarise(AUC = mean(AUC, na.rm = TRUE), .groups = 'drop')

    fit.rf[["results"]] <- fit.rf[["results"]] %>%
      left_join(auc_values, by = "mtry")

    #Tuning parameter (select combination with top AUC)
    tune = which.max(fit.rf$results$AUC)
    fit.rf$bestTune = fit.rf$bestTune %>%
      mutate(mtry = fit.rf$results$mtry[tune])

    #Configure resamples to have the AUCs only using tuned parameter
    fit.rf$resample = fit.rf$resample[order(fit.rf$resample$Resample),] #Order resamples (just in case) to match with correct AUCs from prediction object
    auc = c()
    for (i in 1:nrow(fit.rf$resample)){
      auc_val = fit.rf$pred %>%
        filter(mtry == as.numeric(fit.rf$bestTune),
               Resample == fit.rf$resample$Resample[i]) %>%
        pull(AUC) %>% #AUC per resample is the same
        unique()
      auc = c(auc, auc_val)
    }

    fit.rf$resample = fit.rf$resample %>%
      mutate(AUC = auc)

    ################################################C5.0

    ## Integrate AUCs into prediction matrix
    fit.c50$pred = fit.c50$pred %>%
      group_by(trials, model, winnow) %>% #Parameters for tunning
      mutate(AUC = calculate_auc_resample(obs, yes)) %>%
      ungroup()

    ## Integrate AUCs into results per parameter
    auc_values = fit.c50$pred %>%
      group_by(trials, model, winnow) %>%
      summarise(AUC = mean(AUC, na.rm = TRUE), .groups = 'drop')

    fit.c50[["results"]] <- fit.c50[["results"]] %>%
      left_join(auc_values, by = c("trials", "model", "winnow"))

    #Tuning parameter (select combination with top AUC)
    tune = which.max(fit.c50$results$AUC)
    fit.c50$bestTune = fit.c50$bestTune %>%
      mutate(trials = fit.c50$results$trials[tune],
             model = fit.c50$results$model[tune],
             winnow = fit.c50$results$winnow[tune])

    #Configure resamples to have the AUCs only using tuned parameter
    fit.c50$resample = fit.c50$resample[order(fit.c50$resample$Resample),] #Order resamples (just in case) to match with correct AUCs from prediction object
    auc = c()
    for (i in 1:nrow(fit.c50$resample)){
      auc_val = fit.c50$pred %>%
        filter(trials == as.numeric(fit.c50$bestTune$trials),
               model == as.character(fit.c50$bestTune$model),
               winnow == as.character(fit.c50$bestTune$winnow),
               Resample == fit.c50$resample$Resample[i]) %>%
        pull(AUC) %>% #AUC per resample is the same
        unique()
      auc = c(auc, auc_val)
    }

    fit.c50$resample = fit.c50$resample %>%
      mutate(AUC = auc)

    ################################################LG

    ## Integrate AUCs into prediction matrix
    fit.glm$pred = fit.glm$pred %>% #Calculate resamples AUC scores
      group_by(Resample) %>%
      mutate(AUC = calculate_auc_resample(obs, yes)) %>%
      ungroup()

    ## Integrate AUCs into resamples matrix
    auc = c()
    for (i in 1:nrow(fit.glm$resample)) {
      auc_val = fit.glm$pred %>%
        filter(Resample == fit.glm$resample$Resample[i]) %>%
        pull(AUC) %>% #AUC per resample is the same
        unique()

      auc = c(auc, auc_val)
    }

    fit.glm$resample = fit.glm$resample %>%
      mutate(AUC = auc)

    ## Integrate average CV AUCs into results
    fit.glm$results = fit.glm$results %>%
      mutate(AUC = mean(fit.glm$resample$AUC))

    ################################################LDA

    ## Integrate AUCs into prediction matrix
    fit.lda$pred = fit.lda$pred %>% #Calculate resamples AUC scores
      group_by(Resample) %>%
      mutate(AUC = calculate_auc_resample(obs, yes)) %>%
      ungroup()

    ## Integrate AUCs into resamples matrix
    auc = c()
    for (i in 1:nrow(fit.lda$resample)) {
      auc_val = fit.lda$pred %>%
        filter(Resample == fit.lda$resample$Resample[i]) %>%
        pull(AUC) %>% #AUC per resample is the same
        unique()

      auc = c(auc, auc_val)
    }

    fit.lda$resample = fit.lda$resample %>%
      mutate(AUC = auc)

    ## Integrate average CV AUCs into results
    fit.lda$results = fit.lda$results %>%
      mutate(AUC = mean(fit.lda$resample$AUC))

    ################################################GLMNET

    ## Integrate AUCs into prediction matrix
    fit.glmnet$pred = fit.glmnet$pred %>% #Calculate resamples AUC scores
      group_by(Resample, alpha, lambda) %>%
      mutate(AUC = calculate_auc_resample(obs, yes)) %>%
      ungroup()

    ## Integrate AUCs into results per parameter
    auc_values = fit.glmnet$pred %>%
      group_by(alpha, lambda) %>%
      summarise(AUC = mean(AUC, na.rm = TRUE), .groups = 'drop')

    fit.glmnet[["results"]] <- fit.glmnet[["results"]] %>%
      left_join(auc_values, by = c("alpha", "lambda"))

    #Tuning parameter (select combination with top AUC)
    tune = which.max(fit.glmnet$results$AUC)
    fit.glmnet$bestTune = fit.glmnet$bestTune %>%
      mutate(alpha = fit.glmnet$results$alpha[tune],
             lambda = fit.glmnet$results$lambda[tune])

    #Configure resamples to have the AUCs only using tuned parameter
    fit.glmnet$resample = fit.glmnet$resample[order(fit.glmnet$resample$Resample),] #Order resamples (just in case) to match with correct AUCs from prediction object
    auc = c()
    for (i in 1:nrow(fit.glmnet$resample)){
      auc_val = fit.glmnet$pred %>%
        filter(alpha == as.numeric(fit.glmnet$bestTune$alpha),
               lambda == as.numeric(fit.glmnet$bestTune$lambda),
               Resample == fit.glmnet$resample$Resample[i]) %>%
        pull(AUC) %>% #AUC per resample is the same
        unique()
      auc = c(auc, auc_val)
    }

    fit.glmnet$resample = fit.glmnet$resample %>%
      mutate(AUC = auc)

    ################################################KNN

    ## Integrate AUCs into prediction matrix
    fit.knn$pred = fit.knn$pred %>% #Calculate resamples AUC scores
      group_by(Resample, k) %>%
      mutate(AUC = calculate_auc_resample(obs, yes)) %>%
      ungroup()

    ## Integrate AUCs into results per parameter
    auc_values = fit.knn$pred %>%
      group_by(k) %>%
      summarise(AUC = mean(AUC, na.rm = TRUE), .groups = 'drop')

    fit.knn[["results"]] <- fit.knn[["results"]] %>%
      left_join(auc_values, by = "k")

    #Tuning parameter (select combination with top AUC)
    tune = which.max(fit.knn$results$AUC)
    fit.knn$bestTune = fit.knn$bestTune %>%
      mutate(k = fit.knn$results$k[tune])

    #Configure resamples to have the AUCs only using tuned parameter
    fit.knn$resample = fit.knn$resample[order(fit.knn$resample$Resample),] #Order resamples (just in case) to match with correct AUCs from prediction object
    auc = c()
    for (i in 1:nrow(fit.knn$resample)){
      auc_val = fit.knn$pred %>%
        filter(k == as.numeric(fit.knn$bestTune$k),
               Resample == fit.knn$resample$Resample[i]) %>%
        pull(AUC) %>% #AUC per resample is the same
        unique()
      auc = c(auc, auc_val)
    }

    fit.knn$resample = fit.knn$resample %>%
      mutate(AUC = auc)

    ################################################CART

    ## Integrate AUCs into prediction matrix
    fit.cart$pred = fit.cart$pred %>% #Calculate resamples AUC scores
      group_by(Resample, cp) %>%
      mutate(AUC = calculate_auc_resample(obs, yes)) %>%
      ungroup()

    ## Integrate AUCs into results per parameter
    auc_values = fit.cart$pred %>%
      group_by(cp) %>%
      summarise(AUC = mean(AUC, na.rm = TRUE), .groups = 'drop')

    fit.cart[["results"]] <- fit.cart[["results"]] %>%
      left_join(auc_values, by = "cp")

    #Tuning parameter (select combination with top AUC)
    tune = which.max(fit.cart$results$AUC)
    fit.cart$bestTune = fit.cart$bestTune %>%
      mutate(cp = fit.cart$results$cp[tune])

    #Configure resamples to have the AUCs only using tuned parameter
    fit.cart$resample = fit.cart$resample[order(fit.cart$resample$Resample),] #Order resamples (just in case) to match with correct AUCs from prediction object
    auc = c()
    for (i in 1:nrow(fit.cart$resample)){
      auc_val = fit.cart$pred %>%
        filter(cp == as.numeric(fit.cart$bestTune$cp),
               Resample == fit.cart$resample$Resample[i]) %>%
        pull(AUC) %>% #AUC per resample is the same
        unique()
      auc = c(auc, auc_val)
    }

    fit.cart$resample = fit.cart$resample %>%
      mutate(AUC = auc)

    ################################################Regularized Lasso

    ## Integrate AUCs into prediction matrix
    fit.lasso$pred = fit.lasso$pred %>% #Calculate resamples AUC scores
      group_by(Resample, alpha, lambda) %>%
      mutate(AUC = calculate_auc_resample(obs, yes)) %>%
      ungroup()

    ## Integrate AUCs into results per parameter
    auc_values = fit.lasso$pred %>%
      group_by(alpha, lambda) %>%
      summarise(AUC = mean(AUC, na.rm = TRUE), .groups = 'drop')

    fit.lasso[["results"]] <- fit.lasso[["results"]] %>%
      left_join(auc_values, by = c("alpha", "lambda"))

    #Tuning parameter (select combination with top AUC)
    tune = which.max(fit.lasso$results$AUC)
    fit.lasso$bestTune = fit.lasso$bestTune %>%
      mutate(alpha = fit.lasso$results$alpha[tune],
             lambda = fit.lasso$results$lambda[tune])

    #Configure resamples to have the AUCs only using tuned parameter
    fit.lasso$resample = fit.lasso$resample[order(fit.lasso$resample$Resample),] #Order resamples (just in case) to match with correct AUCs from prediction object
    auc = c()
    for (i in 1:nrow(fit.lasso$resample)){
      auc_val = fit.lasso$pred %>%
        filter(alpha == as.numeric(fit.lasso$bestTune$alpha),
               lambda == as.numeric(fit.lasso$bestTune$lambda),
               Resample == fit.lasso$resample$Resample[i]) %>%
        pull(AUC) %>% #AUC per resample is the same
        unique()
      auc = c(auc, auc_val)
    }

    fit.lasso$resample = fit.lasso$resample %>%
      mutate(AUC = auc)

    ################################################Ridge regression

    ## Integrate AUCs into prediction matrix
    fit.ridge$pred = fit.ridge$pred %>% #Calculate resamples AUC scores
      group_by(Resample, alpha, lambda) %>%
      mutate(AUC = calculate_auc_resample(obs, yes)) %>%
      ungroup()

    ## Integrate AUCs into results per parameter
    auc_values = fit.ridge$pred %>%
      group_by(alpha, lambda) %>%
      summarise(AUC = mean(AUC, na.rm = TRUE), .groups = 'drop')

    fit.ridge[["results"]] <- fit.ridge[["results"]] %>%
      left_join(auc_values, by = c("alpha", "lambda"))

    #Tuning parameter (select combination with top AUC)
    tune = which.max(fit.ridge$results$AUC)
    fit.ridge$bestTune = fit.ridge$bestTune %>%
      mutate(alpha = fit.ridge$results$alpha[tune],
             lambda = fit.ridge$results$lambda[tune])

    #Configure resamples to have the AUCs only using tuned parameter
    fit.ridge$resample = fit.ridge$resample[order(fit.ridge$resample$Resample),] #Order resamples (just in case) to match with correct AUCs from prediction object
    auc = c()
    for (i in 1:nrow(fit.ridge$resample)){
      auc_val = fit.ridge$pred %>%
        filter(alpha == as.numeric(fit.ridge$bestTune$alpha),
               lambda == as.numeric(fit.ridge$bestTune$lambda),
               Resample == fit.ridge$resample$Resample[i]) %>%
        pull(AUC) %>% #AUC per resample is the same
        unique()
      auc = c(auc, auc_val)
    }

    fit.ridge$resample = fit.ridge$resample %>%
      mutate(AUC = auc)

    ################################################SVM radial

    ## Integrate AUCs into prediction matrix
    fit.svm_radial$pred = fit.svm_radial$pred %>% #Calculate resamples AUC scores
      group_by(Resample, sigma, C) %>%
      mutate(AUC = calculate_auc_resample(obs, yes)) %>%
      ungroup()

    ## Integrate AUCs into results per parameter
    auc_values = fit.svm_radial$pred %>%
      group_by(sigma, C) %>%
      summarise(AUC = mean(AUC, na.rm = TRUE), .groups = 'drop')

    fit.svm_radial[["results"]] <- fit.svm_radial[["results"]] %>%
      left_join(auc_values, by = c("sigma", "C"))

    #Tuning parameter (select combination with top AUC)
    tune = which.max(fit.svm_radial$results$AUC)
    fit.svm_radial$bestTune = fit.svm_radial$bestTune %>%
      mutate(sigma = fit.svm_radial$results$sigma[tune],
             C = fit.svm_radial$results$C[tune])

    #Configure resamples to have the AUCs only using tuned parameter
    fit.svm_radial$resample = fit.svm_radial$resample[order(fit.svm_radial$resample$Resample),] #Order resamples (just in case) to match with correct AUCs from prediction object
    auc = c()
    for (i in 1:nrow(fit.svm_radial$resample)){
      auc_val = fit.svm_radial$pred %>%
        filter(sigma == as.numeric(fit.svm_radial$bestTune$sigma),
               C == as.numeric(fit.svm_radial$bestTune$C),
               Resample == fit.svm_radial$resample$Resample[i]) %>%
        pull(AUC) %>% #AUC per resample is the same
        unique()
      auc = c(auc, auc_val)
    }

    fit.svm_radial$resample = fit.svm_radial$resample %>%
      mutate(AUC = auc)

    ################################################SVM linear

    ## Integrate AUCs into prediction matrix
    fit.svm_linear$pred = fit.svm_linear$pred %>% #Calculate resamples AUC scores
      group_by(Resample, C) %>%
      mutate(AUC = calculate_auc_resample(obs, yes)) %>%
      ungroup()

    ## Integrate AUCs into results per parameter
    auc_values = fit.svm_linear$pred %>%
      group_by(C) %>%
      summarise(AUC = mean(AUC, na.rm = TRUE), .groups = 'drop')

    fit.svm_linear[["results"]] <- fit.svm_linear[["results"]] %>%
      left_join(auc_values, by = "C")

    #Tuning parameter (select combination with top AUC)
    tune = which.max(fit.svm_linear$results$AUC)
    fit.svm_linear$bestTune = fit.svm_linear$bestTune %>%
      mutate(C = fit.svm_linear$results$C[tune])

    #Configure resamples to have the AUCs only using tuned parameter
    fit.svm_linear$resample = fit.svm_linear$resample[order(fit.svm_linear$resample$Resample),] #Order resamples (just in case) to match with correct AUCs from prediction object
    auc = c()
    for (i in 1:nrow(fit.svm_linear$resample)){
      auc_val = fit.svm_linear$pred %>%
        filter(C == as.numeric(fit.svm_linear$bestTune$C),
               Resample == fit.svm_linear$resample$Resample[i]) %>%
        pull(AUC) %>% #AUC per resample is the same
        unique()
      auc = c(auc, auc_val)
    }

    fit.svm_linear$resample = fit.svm_linear$resample %>%
      mutate(AUC = auc)

  }

  ###Prediction with best tuned hyper-parameters

  ###Bagged CART

  predictions.bag <- data.frame(predict(fit.treebag, newdata = train_data, type = "prob")) %>% #Predictions using tuned model
    dplyr::select(yes) %>%
    dplyr::rename(BAG = yes)

  ###Random Forest

  predictions.rf = data.frame(predict(fit.rf, newdata = train_data, type = "prob"))[,"yes", drop=F]  %>%
    dplyr::select(yes) %>%
    dplyr::rename(RF = yes) #Predictions of model (already ordered)

  ###C5.0

  predictions.c50 = data.frame(predict(fit.c50$finalModel, newdata = train_data, type = "prob"))[,"yes", drop=F]  %>%
    dplyr::select(yes) %>%
    dplyr::rename(C50 = yes)  #Predictions of model (already ordered)

  ### LG

  predictions.glm = predict(fit.glm, newdata = train_data, type = "prob")[,"yes", drop=F]  %>%
    dplyr::select(yes) %>%
    dplyr::rename(GLM = yes)  #Predictions of model (already ordered)

  ### LDA

  predictions.lda = predict(fit.lda, newdata = train_data, type = "prob")[,"yes", drop=F]  %>%
    dplyr::select(yes) %>%
    dplyr::rename(LDA = yes)  #Predictions of model (already ordered)

  ### GLMNET

  predictions.glmnet = predict(fit.glmnet, newdata = train_data, type = "prob")[,"yes", drop=F]  %>%
    dplyr::select(yes) %>%
    dplyr::rename(GLMNET = yes)  #Predictions of model (already ordered)

  ### KNN

  predictions.knn = predict(fit.knn, newdata = train_data, type = "prob")[,"yes", drop=F]  %>%
    dplyr::select(yes) %>%
    dplyr::rename(KNN = yes) #Predictions of model (already ordered)

  ## CART

  predictions.cart = predict(fit.cart, newdata = train_data, type = "prob")[,"yes", drop=F]  %>%
    dplyr::select(yes) %>%
    dplyr::rename(CART = yes)  #Predictions of model (already ordered)

  ## Regularized Lasso

  predictions.lasso = predict(fit.lasso, newdata = train_data, type = "prob")[,"yes", drop=F]  %>%
    dplyr::select(yes) %>%
    dplyr::rename(LASSO = yes)  #Predictions of model (already ordered)

  ## Ridge regression

  predictions.ridge = predict(fit.ridge, newdata = train_data, type = "prob")[,"yes", drop=F]  %>%
    dplyr::select(yes) %>%
    dplyr::rename(RIDGE = yes)  #Predictions of model (already ordered)

  ## SVM radial

  predictions.svm_radial = predict(fit.svm_radial, newdata = train_data, type = "prob")[,"yes", drop=F]  %>%
    dplyr::select(yes) %>%
    dplyr::rename(SVM_radial = yes)  #Predictions of model (already ordered)

  ## SVM linear

  predictions.svm_linear = predict(fit.svm_linear, newdata = train_data, type = "prob")[,"yes", drop=F]  %>%
    dplyr::select(yes) %>%
    dplyr::rename(SVM_linear = yes)  #Predictions of model (already ordered)

  ############################################################## Save models

  ensembleResults <- list(BAG = fit.treebag,
                          RF = fit.rf,
                          C50 = fit.c50,
                          GLM = fit.glm,
                          LDA = fit.lda,
                          KNN = fit.knn,
                          CART = fit.cart,
                          GLMNET = fit.glmnet,
                          LASSO = fit.lasso,
                          RIDGE = fit.ridge,
                          SVM_radial = fit.svm_radial,
                          SVM_linear = fit.svm_linear)


  model_predictions = list(BAG = predictions.bag,
                           RF = predictions.rf,
                           C50 = predictions.c50,
                           GLM = predictions.glm,
                           LDA = predictions.lda,
                           KNN = predictions.knn,
                           CART = predictions.cart,
                           GLMNET = predictions.glmnet,
                           LASSO = predictions.lasso,
                           RIDGE = predictions.ridge,
                           SVM_radial = predictions.svm_radial,
                           SVM_linear = predictions.svm_linear)

  #Remove models with same predictions across samples (not able to make distinction)
  model_predictions <- lapply(model_predictions, function(df) {
    df = df %>%
      select(where(~ n_distinct(.) > 1))

    if(ncol(df) == 0){
      df = NULL
    }

    return(df)
  })

  model_predictions = Filter(Negate(is.null), model_predictions) #Discard not useful predictions
  ensembleResults = ensembleResults[names(model_predictions)] #Discard not useful models based on predictions

  model_predictions = do.call(cbind, model_predictions) #Join as data frame

  #Clean memory
  rm(fit.treebag, fit.rf, fit.c50, fit.glm, fit.lda, fit.knn, fit.cart, fit.glmnet, fit.lasso, fit.ridge, fit.svm_radial, fit.svm_linear, multifolds)
  gc()

  if(stacking){
    features = colnames(train_data)[colnames(train_data) != "target"]

    #Base models using ML models with best accuracy or AUC from each family
    if(metric == "Accuracy"){
      base_models = compute_cv_accuracy(ensembleResults, base_models = T, file_name = file_name, return = return)
    }else if(metric == "AUC"){
      base_models = compute_cv_AUC(ensembleResults, base_models = T, file_name = file_name, return = return)
    }

    #Save variable importance of each base model
    importance = list()
    for (i in 1:length(base_models$Base_models)) {
      importance[[i]] = varImp(ensembleResults[[base_models$Base_models[i]]], scale = F)
      names(importance)[i] = base_models$Base_models[i]
    }

    features_predictions = model_predictions %>%
      t() %>%
      data.frame() %>%
      rownames_to_column("Models") %>%
      filter(grepl(paste0("\\b(", paste(base_models$Base_models, collapse = "|"), ")\\b"), Models)) %>%
      column_to_rownames("Models") %>%
      t() %>%
      data.frame()

    meta_features = cbind(features_predictions, "true_label" = train_data$target)

    meta_learner <- train(true_label ~ ., data = meta_features, method = "glmnet", trControl = trainControl) #Staking based on simple logistic regression

    #Base models using ALL ML models
    meta_features_all = cbind(model_predictions, "true_label" = train_data$target)

    meta_learner_all <- train(true_label ~ ., data = meta_features_all, method = "glmnet", trControl = trainControl) #Staking based on simple logistic regression

    cat("Meta-learners ML model based on GLM\n")
    output = list("Features" = features, "Meta_learner" = meta_learner, "Base_models" = base_models$Base_models, "ML_models" = ensembleResults, "Variable_importance" = importance)

  }else{
    features = colnames(train_data)[colnames(train_data) != "target"] #Extract features used for model training

    #Top model with best accuracy or AUC
    if(metric == "Accuracy"){
      metrics = compute_cv_accuracy(ensembleResults, file_name, file_name = file_name, return = return)
    }else if(metric == "AUC"){
      metrics = compute_cv_AUC(ensembleResults, file_name, file_name = file_name, return = return)
    }

    top_model = metrics[["Top_model"]]

    model = ensembleResults[[top_model]]

    cat("Best ML model found: ", top_model, "\n")

    cat("Returning model trained\n")

    output = list("Features" = features, "Model" = model, "ML_Models" = ensembleResults)
  }


  return(output)

}

#' Compute Leaving-one-dataset-out approach
#'
#'If we are training different cohorts and trying to predict in independent ones, CellTFusion needs to be compute individually to avoid data leakage. For this, the implemented function eases this workflow to the user.
#'
#' @param raw.counts A matrix with the raw counts with samples as columns and genes symbols as rows.
#' @param normalized Boolean value to specify if raw_counts need to be normalized (If no raw_counts are available and argument corresponds to already normalized counts this arguments needs to be set to False)
#' @param clinical A data frame with the clinical data containing the target variable to analyze
#' @param trait Column name of target variable
#' @param trait.positive Value of target consider as positive for prediction
#' @param trait.out Column name indicating the column from which dataset need to be subset (e.g. cohort = c(A, B, C, D, A, A) this parameter should be "cohort")
#' @param out Value to use for subsetting the matrix. If cohort = "A" is going to be leave out, this parameter should be "A".
#' @param metric Metric to based the methods choice during cross validation (CV). Either accuracy or AUC
#' @param stack If stacking needs to be applied
#' @param feature.selection If feature selection using Boruta algorithm needs to be computed
#' @param doParallel Whether to do or not parallelization. Only CBSX and DWLS methods will run in parallel.
#' @param workers Number of processes available to run on parallel. If no number is set, this will correspond to detectCores() - 1
#' @param deconv_methods A character vector with the deconvolution methods to run. Default are "Quantiseq", "MCP", "xCell", "CBSX", "Epidish", "DeconRNASeq", "DWLS"
#' @param file_name File name for the csv files and plots saved in the Results/ directory
#' @param return Return the plots printed during each iteration
#'
#' @return A list containing
#'
#' - Trained model
#' - Features used to train the model
#' - Cell groups
#' - Prediction metrics
#' - AUROC and AUPRC
#' - Predictions
#'
#' @export
#'
#' @examples
#'
#' res_ml = compute.LODO.ML(raw.counts, normalized = T, clinical, trait = "Response",trait.positive = "R", trait.out = "Cohort", out = "Dupont", metric = "Accuracy", stack = T, feature.selection = F, doParallel = T, workers = 4, file_name = "Test", return = F)
#'
compute.LODO.ML = function(raw.counts, normalized = F, clinical, trait, trait.positive, trait.out, out, metric = "Accuracy", stack, feature.selection = T, doParallel = F, workers = NULL, deconv_methods = c("Quantiseq", "MCP", "xCell", "CBSX", "Epidish", "DeconRNASeq", "DWLS"), file_name = NULL, return = T){

  clinical = clinical %>%
    mutate(trait.out = clinical[,trait.out]) ##Just a way to make work "filter" cause it does not allow "" variables (might change after)

  # Normalize counts
  if(normalized == T){
    norm.counts = data.frame(ADImpute::NormalizeTPM(raw.counts))
  }else{ #If they are already normalized
    norm.counts = raw.counts
  }

  # Test cohort
  traitData_test = clinical %>%
    filter(trait.out == out)
  raw.counts_test = raw.counts[,colnames(raw.counts)%in%rownames(traitData_test)]
  counts.normalized_test = norm.counts[,colnames(norm.counts)%in%rownames(traitData_test)]

  tfs_test = compute.TFs.activity(counts.normalized_test)
  deconv_test = compute.deconvolution(raw.counts_test, normalized = normalized, methods = deconv_methods, doParallel = doParallel, workers = workers, credentials.mail = "marcelo.hurtado@inserm.fr", credentials.token = "734212f6ad77fc4eea2bdb502792f294")

  # Train cohort
  traitData_train = clinical %>%
    filter(trait.out != out)
  raw.counts_train = raw.counts[,colnames(raw.counts)%in%rownames(traitData_train)]
  counts.normalized_train = norm.counts[,colnames(norm.counts)%in%rownames(traitData_train)]

  tfs_train = compute.TFs.activity(counts.normalized_train)
  deconv_train = compute.deconvolution(raw.counts_train, normalized = normalized, methods = deconv_methods, doParallel = doParallel, workers = workers, credentials.mail = "marcelo.hurtado@inserm.fr", credentials.token = "734212f6ad77fc4eea2bdb502792f294")

  ###############################################################################################################################################################################

  # CellTFusion
  network = compute.WTCNA(tfs_train, corr_mod = 0.9, clustering.method = "ward.D2", return = return)
  pathways = compute.pathway.activity(counts.normalized_train, paths = NULL)
  tfs.modules.clusters = compute.TF.network.classification(network, pathways, return = return)
  dt = compute.deconvolution.analysis(deconv_train, corr = 0.8, seed = 123)
  corr_modules = compute.modules.relationship(network[[1]], dt[[1]], return = T, plot = return)
  cell_dendrograms = identify.cell.groups(corr_modules, tfs.modules.clusters, height = 20, return = return)
  cell.groups = cell.groups.computation(dt[[1]], tfs.module.network = network, cell.dendrograms = cell_dendrograms) #Identify cell groups with specific cut

  ###############################################################################################################################################################################

  ####################################################Training

  #Set training set
  train_data = cell.groups[[1]] %>%
    data.frame() %>%
    mutate(Trait = traitData_train[,trait],
           target = as.factor(ifelse(Trait == trait.positive, 'yes', 'no'))) %>%
    dplyr::select(-Trait)

  train_data$target <- factor(train_data$target, levels = c("no", "yes"))  # Just in case positive class is not well defined

  #Cross-validation training (5 k-folds and 100 repetitions)
  training = compute.k_fold_CV(train_data, k_folds = 5, n_rep = 100,  metric = metric, stacking = stack, boruta = feature.selection, boruta_iterations = 100, fix_boruta = F, boruta_threshold = 0.8, file_name = file_name, return= return)

  ####################################################Predicting
  if(length(training)!=0){
    cell_groups = cell.groups #Save cell groups scores per partition
    features = training[["Features"]] #Save selected features per partition
    ####################### Testing set
    testing_set = compute_cell_groups_signatures(dt, network, cell.groups, features, deconv_test, tfs_test) #Compute features in testing set
    #Extract target variable
    target = traitData_test %>%
      mutate(target = ifelse(traitData_test[,trait] == trait.positive, "yes", "no")) %>%
      pull(target)

    target = factor(target, levels = c("no", "yes"))

    if(stack){
      model = training[["Meta_learner"]]
      prediction = compute.prediction.stacked(model, testing_set, target, training[["ML_models"]], training[["Base_models"]])
    }else{
      model = training[["Model"]] #Save best ML model based on the Accuracy from CV per partition
      prediction = compute.prediction(model, testing_set, target)
    }

    auc_roc_score = prediction[["AUC"]][["AUROC"]]
    auc_prc_score = prediction[["AUC"]][["AUPRC"]]

    metrics = prediction[["Metrics"]]
    predictions = prediction[["Predictions"]]

    if(return == T){
      get_curves(metrics, "specificity", "sensitivity", "recall", "precision", "model", auc_roc_score, auc_prc_score, file_name)
    }

    rm(network, pathways, tfs.modules.clusters, dt, corr_modules, cell_dendrograms, cell.groups,
       traitData_train, counts.normalized_train, tfs_train, deconv_train,
       traitData_test, counts.normalized_test, deconv_test, tfs_test,
       clinical, norm.counts) #Remove variables

    gc() #Clean garbage

    return(list(Model = model, Features = features, Cell_groups = cell_groups, Prediction_metrics = metrics, AUC = list(AUROC = auc_roc_score, AUPRC = auc_prc_score), Prediction = predictions))
  }else{  #No features are selected as predictive

    rm(network, pathways, tfs.modules.clusters, dt, corr_modules, cell_dendrograms, cell.groups,
       traitData_train, counts.normalized_train, tfs_train, deconv_train,
       traitData_test, counts.normalized_test, deconv_test, tfs_test,
       clinical, norm.counts) #Remove variables

    gc() #Clean garbage

    return(NULL)
  }

}


#' Calculates accuracy values from prediction
#'
#' @param metrics A Dataframe with metrics obtained using get_sensitivity_specificity()
#' @param target A character vector containing the true values from the target variable
#'
#' @return A numeric vector with the accuracy values
#' @export
#'
#' @examples
#'
#'observed_values = c("yes", "yes", "no", "yes")
#'accuracy = calculate_accuracy(metrics, observed_values)
#'
calculate_accuracy <- function(metrics, target) {
  sensitivity = metrics[,"sensitivity"]
  specificity = metrics[,"specificity"]
  total_positives = sum(target == "yes")
  total_negatives = sum(target == "no")
  TP <- sensitivity * total_positives
  FN <- total_positives - TP
  TN <- specificity * total_negatives
  FP <- total_negatives - TN

  accuracy <- (TP + TN) / (TP + TN + FP + FN)

  return(accuracy)
}

#' Calculates AUC values for each resampling during cross-validation (CV)
#'
#' @param obs A character vector with observed values (grountruth)
#' @param pred A character vector with predicted values
#'
#' @return A numeric value corresponding to the AUC score
#' @export
#'
#' @examples
#'
#' obs = c("yes", "no", "yes", "no")
#' yes = c("yes", "yes", "yes", "no")
#' AUC = calculate_auc_resample(obs, yes)
#'
calculate_auc_resample = function(obs, pred){

  prob_obs = data.frame("yes" = pred, "obs" = obs)

  prob_obs = prob_obs %>%
    arrange(desc(pred)) %>% #need to be arrange for apply cumulative sum
    mutate(is_yes = (obs == "yes"),
           tp = cumsum(is_yes), #true positive above the threshold - cumulative sum to refer to the threshold
           fp = cumsum(!is_yes), #false positive above the threshold - cumulative sum to refer to the threshold
           fpr = fp/sum(obs == 'no'),
           tpr = tp/sum(obs == 'yes'))

  auc_value = calculate_auc(prob_obs$fpr, prob_obs$tpr)

  return(auc_value)
}

#' Calculates AUC from precision-recall curve
#'
#' @param recall A numeric column with the recall values obtained from get_sensitivity_specificity()
#' @param precision A numeric column with the precision values obtained from get_sensitivity_specificity()
#'
#' @return A numeric value corresponding to the AUC from the precision-recall curve
#' @export
#'
#' @examples
#'
#' sens_spec = get_sensitivity_specificity(predict, target, model)
#' auprc = calculate_auprc(sens_spec$recall, sens_spec$precision)
#'
calculate_auprc <- function(recall, precision) {
  # Sort by Recall to ensure trapezoidal rule is correctly applied
  ordered <- order(recall)
  recall <- recall[ordered]
  precision <- precision[ordered]

  auprc <- 0
  for (i in 1:(length(recall) - 1)) { # -1 to avoid NA from the last terms
    # Trapezoidal rule: (precision[i] + precision[i+1]) / 2 * (recall[i+1] - recall[i])
    auprc <- auprc + ((precision[i+1] + precision[i]) / 2) * (recall[i+1] - recall[i])
  }
  return(auprc)
}

#' Calculates AUC from ROC curve
#'
#' @param fpr A numeric column with the false positive rate values obtained from get_sensitivity_specificity()
#' @param sensitivity A numeric column with the sensitivity values obtained from get_sensitivity_specificity()
#'
#' @return A numeric value corresponding to the AUC from the precision-recall curve
#' @export
#'
#' @examples
#'
#' sens_spec = get_sensitivity_specificity(predict, target, model)
#' auroc = calculate_auroc(sens_spec$fpr, sens_spec$sensitivity)
#'
calculate_auroc <- function(fpr, sensitivity) {
  #tpr is sensitivity

  # Sort by FPR to ensure trapezoidal rule is correctly applied (Already ordered)
  ordered <- order(fpr)
  fpr <- fpr[ordered]
  sensitivity <- sensitivity[ordered]

  auc <- 0
  for (i in 1:(length(fpr) - 1)) { #-1 to avoid NA cause last terms are TPR = 1 and FPR = 1
    # Trapezoidal rule: (TPR_i + TPR_{i+1}) / 2 * (FPR_{i+1} - FPR_i)
    auc <- auc + ((sensitivity[i+1] + sensitivity[i]) / 2) * (fpr[i+1] - fpr[i])
  }
  return(auc)

}

#' Calculate confusion values
#'
#' @param metrics A Dataframe with metrics obtained using get_sensitivity_specificity()
#' @param target A character vector containing the true values from the target variable
#'
#' @return A list containing the True positives (TP), False negatives (FN), True negatives (TN) and False positives (FP)
#' @export
#'
#' @examples
#'
#' target = c("yes", "no", "yes", "no", "no")
#' metrics = get_sensitivity_specificity(predict, target, model)
#' confusion_values = calculate_confusion_values(metrics, target)
#'
calculate_confusion_values <- function(metrics, target) {
  sensitivity <- metrics[,"sensitivity"]
  specificity <- metrics[,"specificity"]

  # Count total positives and negatives
  total_positives <- sum(target == "yes")
  total_negatives <- sum(target == "no")

  # Calculate confusion matrix values
  TP <- sensitivity * total_positives
  FN <- total_positives - TP
  TN <- specificity * total_negatives
  FP <- total_negatives - TN

  return(list(TP = TP, FN = FN, TN = TN, FP = FP))
}

#' Compute weighted feature importance from base models and meta-learner for stacking models
#'
#' @param base_importance A matrix of feature importance from base models
#' @param base_models A character vector with chosen base models
#' @param meta_learner A caret object with the meta-learner model trained using base models
#'
#' @return A matrix of feature importance weighted from the base models and the meta-learner
#' @export
#'
#' @examples
#'
#' var_importance = calculate_feature_importance_stacking(variable_importance, base_models, meta-learner)
#'
calculate_feature_importance_stacking = function(base_importance, base_models, meta_learner){

  #Extract features importance values within each base model for the meta-learner
  base_importance_list = list()
  for (i in 1:length(base_models)) {
    check = ncol(base_importance[[base_models[i]]][["importance"]])
    if(check > 1){ #Means importance is given for each class
      base_importance_list[[i]] = base_importance[[base_models[i]]][["importance"]] %>%
        rownames_to_column("features") %>%
        dplyr::select(features, yes) %>% #Take only importance for positive class
        dplyr::rename(importance = yes)
    }else{
      base_importance_list[[i]] = base_importance[[base_models[i]]][["importance"]] %>%
        rownames_to_column("features") %>%
        dplyr::rename(importance = Overall)
    }
    names(base_importance_list)[i] = base_models[i]
  }

  #Combine all base model importances in one data frame and add the model name
  combined_importance <- bind_rows(
    lapply(names(base_importance_list), function(model) {
      base_importance_list[[model]] %>%
        data.frame() %>%
        mutate(model = model)
    })
  )

  #Calculate base-models importance for the meta-learner
  meta_importance = varImp(meta_learner, scale = F)$importance %>%
    rownames_to_column("model")

  #Normalize the meta-learner's importance scores so they sum to 1
  meta_importance$Overall <- meta_importance$Overall / sum(meta_importance$Overall)

  #Combine features importance within base models with the overall importance for meta-learner
  combined_importance <- combined_importance %>%
    left_join(meta_importance, by = "model") %>%
    mutate(weighted_importance = importance * Overall) # importance is from base, Overall is from meta

  #Sum the weighted importance by feature across all models
  final_importance <- combined_importance %>%
    group_by(features) %>%
    summarise(final_importance = sum(weighted_importance, na.rm = TRUE)) %>%
    arrange(desc(final_importance))

  return(final_importance)

}

#' Calculate precision values
#'
#' @param metrics A Dataframe with metrics obtained using get_sensitivity_specificity()
#' @param target A character vector containing the true values from the target variable
#'
#' @return A numeric vector with precision values
#' @export
#'
#' @examples
#'
#' observed = c("yes", "no", "yes", "no", "no")
#' metrics = get_sensitivity_specificity(predict, target, model)
#' precision = calculate_precision(metrics, observed)
#'
calculate_precision <- function(metrics, target) {
  confusion_values <- calculate_confusion_values(metrics, target)
  TP <- confusion_values$TP
  FP <- confusion_values$FP

  precision <- TP / (TP + FP)

  return(precision)
}

#' Cakculate recall values
#'
#' @param metrics A Dataframe with metrics obtained using get_sensitivity_specificity()
#' @param target A character vector containing the true values from the target variable
#'
#' @return A numeric vector with recall values
#' @export
#'
#' @examples
#'
#' observed = c("yes", "no", "yes", "no", "no")
#' metrics = get_sensitivity_specificity(predict, target, model)
#' recall = calculate_recall(metrics, observed)
#'
calculate_recall <- function(metrics, target) {
  confusion_values <- calculate_confusion_values(metrics, target)
  TP <- confusion_values$TP
  FN <- confusion_values$FN

  # Calculate recall (sensitivity)
  recall <- TP / (TP + FN)

  return(recall)
}

#' Choose three base models for stacking based either in Accuracy or AUC scores
#'
#' @param models List with trained machine learning models
#' @param metric Metric to choose the top base models (either Accuracy or AUC). Default is Accuracy
#'
#' @return A character vector with the base models
#' @export
#'
#' @examples
#'
#' models <- list(BAG = fit.treebag,RF = fit.rf, C50 = fit.c50, GLM = fit.glm, LDA = fit.lda, KNN = fit.knn, CART = fit.cart)
#' base_models = choose_base_models(models, metric = "Accuracy")
#'
choose_base_models = function(models, metric = "Accuracy"){

  #Bind accuracy values from each model
  resample_df = list()
  for (i in 1:length(models)){
    resample_df[[i]] = models[[i]]$resample %>%
      mutate(model = names(models)[i])
    names(resample_df)[i] = names(models)[i]
  }
  resample_df = do.call(rbind, resample_df)

  if(metric == "Accuracy"){
    #Prepare data frame for ploting
    resample_df <- resample_df %>%
      group_by(model) %>%
      summarise(Accuracy = mean(Accuracy))
  }else if(metric == "AUC"){
    #Prepare data frame for ploting
    resample_df <- resample_df %>%
      group_by(model) %>%
      summarise(AUC = mean(AUC))
  }

  resample_df <- resample_df %>%
    mutate(Category = case_when(
      model %in% c("BAG", "C50", "CART", "RF") ~ "Tree-based Methods",
      model %in% c("GLM", "LDA", "GLMNET", "LASSO", "RIDGE") ~ "Linear Models",
      model %in% c("KNN", "SVM_linear", "SVM_radial") ~ "Instance-based Methods",
      TRUE ~ "Other"  # In case there are models not in the above lists
    ))

  if(metric == "Accuracy"){
    groupped_df <- resample_df %>%
      group_by(Category) %>%
      filter(Accuracy == max(Accuracy)) %>%
      ungroup()
  }else if(metric == "AUC"){
    groupped_df <- resample_df %>%
      group_by(Category) %>%
      filter(AUC == max(AUC)) %>%
      ungroup()
  }else{
    stop("No metric defined")
  }

  #Retrieve top model based on accuracy/auc
  base_models <- groupped_df %>%
    pull(model)

  return(base_models)
}

#' Compute Accuracy values for each machine learning model
#'
#' @param models List of trained models
#' @param file_name (Optional) File name for plot
#' @param base_models Boolean value to specify if base models for stacking need to be chosen or not
#' @param return Boolean value to specify if plot the accuracy values across models should be saved or not
#'
#' @return A list containing
#'
#' - Accuracy values for each model
#' - Top model with best accuracy
#' - If base_models = T, it returns a character vector with the chosen base models (see choose_base_models())
#'
#' @export
#'
#' @examples
#'
#' res = compute_cv_accuracy(ml_models, base_models = T, file_name = "Test", return = T)
#'
compute_cv_accuracy = function(models, file_name = NULL, base_models = F, return = T){

  #Bind accuracy values from each model
  accuracy = list()
  for (i in 1:length(models)){
    accuracy[[i]] = models[[i]]$resample %>%
      mutate(model = names(models)[i])
    names(accuracy)[i] = names(models)[i]
  }
  accuracy_data = do.call(rbind, accuracy)

  #Retrieve top model based on accuracy
  res_accuracy <- accuracy_data %>%
    group_by(model) %>%
    summarise(Accuracy = mean(Accuracy)) %>%
    arrange(desc(Accuracy))

  top_model = res_accuracy %>%
    slice(1) %>%
    pull(model)

  if(return){
    pdf(paste0("Results/Accuracy_CV_methods_", file_name, ".pdf"), width = 10)
    plot(ggplot(accuracy_data, aes(x = model, y = Accuracy, fill = model)) +
           geom_boxplot() +
           labs(title = "Distribution of Accuracy Values by Model",
                x = "Model",
                y = "Accuracy") +
           theme_minimal() +
           theme(legend.position = "none"))
    dev.off()
  }

  if(base_models == T){
    cat("Choosing base models for stacking.......................................\n\n")
    base_models = choose_base_models(models, metric = "Accuracy")
    cat("Models chosen are:", paste0(base_models, collapse = ", "), "\n\n")
    return(list("Accuracy" = res_accuracy, "Top_model" = top_model, "Base_models" = base_models))
  }else{
    return(list("Accuracy" = res_accuracy, "Top_model" = top_model))
  }

}

#' Compute AUC values for each machine learning model
#'
#' @param models List of trained models
#' @param file_name (Optional) File name for plot
#' @param base_models Boolean value to specify if base models for stacking need to be chosen or not
#' @param return Boolean value to specify if plot the AUC values across models should be saved or not
#'
#' @return A list containing
#'
#' - AUC values for each model
#' - Top model with best AUC
#' - If base_models = T, it returns a character vector with the chosen base models (see choose_base_models())
#'
#' @export
#'
#' @examples
#'
#' res = compute_cv_AUC(ml_models, base_models = T, file_name = "Test", return = T)
#'
compute_cv_AUC = function(models, file_name = NULL, base_models = F, return = T){

  #Bind AUC values from each model
  auc = list()
  for (i in 1:length(models)){
    auc[[i]] = models[[i]]$resample %>% #we use the resample matrix and not directly the results matrix as some have hyperparameters so we will need to define best on the tuned parameter (=more code) - resample matrix is made based on the best tuning
      mutate(model = names(models)[i])
    names(auc)[i] = names(models)[i]
  }
  auc_data = do.call(rbind, auc)

  #Retrieve top model based on accuracy
  res_auc <- auc_data %>%
    group_by(model) %>%
    summarise(AUC = mean(AUC))  %>%
    arrange(desc(AUC))

  top_model = res_auc %>%
    slice(1) %>%
    pull(model)

  if(return){
    pdf(paste0("Results/AUC_CV_methods_", file_name, ".pdf"), width = 10)
    plot(ggplot(auc_data, aes(x = model, y = AUC, fill = model)) +
           geom_boxplot() +
           labs(title = "Distribution of AUC scores by Model",
                x = "Model",
                y = "AUC") +
           theme_minimal() +
           theme(legend.position = "none"))
    dev.off()
  }

  if(base_models == T){
    cat("Choosing base models for stacking.......................................\n\n")
    base_models = choose_base_models(models, metric = "AUC")
    cat("Models chosen are:", paste0(base_models, collapse = ", "), "\n\n")
    return(list("AUC" = res_auc, "Top_model" = top_model, "Base_models" = base_models))
  }else{
    return(list("AUC" = res_auc, "Top_model" = top_model))
  }

}

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

#' Pool performance curves
#'
#' Plot boxplots from pooled metrics from different iterations, mean AUC-ROC and AUC-PRC are calculated.
#'
#' @param file.name String with the file name to saved the plots
#'
#' @return Boxplots in the Results/ directory.
#' @export
#'
#' @examples
#'
#' get_pooled_roc_curves("Test")
#'
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

#' Calculate prediction metrics
#'
#' Calculate sensitivity and specificity values from prediction values
#'
#' @param predictions A character vector with the prediction values
#' @param observed A character vector with the true labels
#' @param ml.model Trained machine learning model
#'
#' @return A matrix with the prediction metrics
#' @export
#'
#' @examples
#'
#' prediction_metrics = get_sensitivity_specificity(predict, target, model)
#'
get_sensitivity_specificity = function(predictions, observed, ml.model){
  prob_obs = bind_cols(predictions, observed = observed)

  prob_obs = prob_obs %>%
    arrange(desc(yes)) %>% #need to be arrange for apply cumulative sum
    mutate(is_yes = (observed == "yes"),
           tp = cumsum(is_yes), #true positive above the threshold - cumulative sum to refer to the threshold
           fp = cumsum(!is_yes), #false positive above the threshold - cumulative sum to refer to the threshold
           sensitivity = tp/sum(observed == 'yes'),
           fpr = fp/sum(observed == 'no'),
           specificity = 1 - fpr) %>%
    select(sensitivity, specificity, fpr) %>%
    mutate(model = ml.model)

  # starts_at_zero <- any(prob_obs$sensitivity == 0 & prob_obs$fpr == 0)

  # ##Add dummy row if it doesnt start at 0
  # if(!starts_at_zero){
  #   dummy_row <- data.frame(
  #     sensitivity = 0,
  #     specificity = 1,
  #     fpr = 0,
  #     model = ml.model
  #   )
  #
  #   prob_obs = rbind(dummy_row, prob_obs)
  # }

  prob_obs = prob_obs %>%
    mutate(Accuracy = calculate_accuracy(., observed),
           precision = calculate_precision(., observed),
           recall = calculate_recall(., observed))


  return(prob_obs)

}

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
