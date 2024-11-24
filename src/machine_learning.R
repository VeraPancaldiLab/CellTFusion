
##Machine learning pipeline functions adapted for CellTFusion (https://github.com/VeraPancaldiLab/CellTFusion)

##Basic function to compute boruta algorithm
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

##Merge results from boruta iterations 
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

##Function for iteratively running boruta (parallelization available)
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

# Create folds for each repetition with different seeds (Not used anymore, replace by Multifolds)
create_folds_for_repetitions <- function(data, k_folds, n_rep) {
  all_folds <- list()
  for (rep in 1:n_rep) {
    set.seed(sample.int(100000, 1)) # Change seed for each repetition
    folds <- createFolds(data$target, k = k_folds, returnTrain = TRUE, list = TRUE)
    all_folds[[rep]] <- folds
  }
  return(all_folds)
}

##Compute boruta during each fold (DEPRECATED - might give overfitting)
computed.boruta.kfolds = function(folds, data_model, boruta_iterations, fix_boruta, tentative, threshold, file_name){
  
  folds_threshold = 0.8*length(folds) 
  features_folds = list()
  
  for (i in 1:length(folds)) {
    message("Feature selection using Boruta...............................................................\n\n")
    training_set = data_model[folds[[i]],]
    res_boruta = feature.selection.boruta(training_set, iterations = boruta_iterations, fix = fix_boruta, thres = threshold, file_name = file_name, return = T)
    
    if(tentative == F){
      if(length(res_boruta$Confirmed) <= 1){
        features_folds[[i]] = list()
        message("\nNo features were confirmed in more than ", round(threshold*100) ,"% of the times for training in this specific fold.......................\n\n")
      }
      message("\nKeeping only features confirmed in more than", round(threshold*100) ,"% of the times for training in this specific fold......................\n\n")
      message("If you want to consider also tentative features, please specify tentative = T in the parameters.\n\n")
      features_folds[[i]] = res_boruta$Confirmed
    }else{
      sum_features = length(res_boruta$Confirmed) + length(res_boruta$Tentative)
      if(sum_features <= 1){
        features_folds[[i]] = list()
        message("\nNo features were confirmed in more than ", round(thresh*100) ,"% of the times for training in this specific fold.......................\n\n")
      }
      message("\nKeeping only features confirmed and tentative in more than", round(thresh*100) ,"% of the times for training in this specific fold............................\n\n")
      features_folds[[i]] = c(res_boruta$Confirmed, res_boruta$Tentative)
    }
  }
  
  all_features <- unlist(features_folds)
  feature_freq <- table(all_features)
  selected_features <- names(feature_freq[feature_freq >= folds_threshold])
  
  if(length(selected_features)<=1){
    message("No features selected meet the requirements. Try with different parameter values.")
  }else{
    return(selected_features)
  }
  
}

#Main function for CV training using 13 ML models
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

#Main ML Pipeline for one partition 
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

#Main ML Pipeline for doing Leaving-one-dataset-out (LODO)

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

unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
  gc()
}

#Main ML pipeline for several partitions
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
    
    cl = parallel::makeCluster(num_cores) #Forking just copy the R session in its current state. - makeCluster() all must be exported (copied) to the clusters, which can add some overhead
    doParallel::registerDoParallel(cl)

    message("\nRunning ", iterations, " splits for training and test using ", num_cores, " cores")

    # Run foreach loop using each random seed directly
    foreach(iteration = seq_len(iterations), random.seed = sample.int(100000, iterations)) %dopar% {
      
      # Use absolute path for the source file to avoid path issues
      source("src/environment_set.R")
      
      # Run foreach loop using each random seed directly
      tryCatch({
        # Compute the result with the current random seed
        result <- compute.ML(
          raw.counts, normalized, clinical, trait, trait.positive, partition, 
          metric, stack, feature.selection, doParallel = F, workers = NULL, 
          seed = random.seed, deconv_methods = deconv_methods, 
          file_name = file.name, return = return
        )
        
        # Save result as RDS file with unique identifier based on iteration (random seed)
        saveRDS(list(result = result, seed = random.seed), 
                file = file.path("Results/ML_models", paste0("ML_result_", iteration, ".rds")))
        
      }, error = function(e) {

        # Save error information as RDS file with random seed identifier
        saveRDS(list(result = NULL, error = e$message, seed = random.seed), 
                file = file.path(paste0("Results/ML_models/ML_result_", iteration, ".rds")))
      })
      
    }
        
    #Stop cluster after all runs
    parallel::stopCluster(cl)
    unregister_dopar() #Stop Dopar from running in the background
    
  }
    
  message("Analysis is done!")
  
  message("ML models are saved in Results/ML_models folder")
  
}

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
    dplyr::summarize(medianAUROC = median(AUC_roc))
  
  mean_auc_prc = data %>%
    group_by(Cohort) %>%
    dplyr::summarize(medianAUPRC = median(AUC_prc))
  
  # Plot boxplot with mean AUC annotations
  plot_roc = ggplot(data, aes(x = Cohort, y = AUC_roc, fill = Cohort)) +
    geom_boxplot() +
    labs(title = paste0("Distribution of AUROC values across ", iterations, " splits"),
         x = "Model",
         y = "AUROC") +
    theme_minimal() +
    theme(legend.position = "right") +
    geom_text(data = mean_auc_roc, aes(x = Cohort, y = max(data$AUC_roc), 
                                   label = paste("Median AUROC:", round(meanAUROC, 3))),
              size = 4, color = "black", vjust = -0.5)
  
  pdf(paste0("Results/Boxplot_AUROC_performance_", file.name, ".pdf"))
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
                                       label = paste("Median AUPRC:", round(meanAUPRC, 3))),
              size = 4, color = "black", vjust = -0.5)
  
  pdf(paste0("Results/Boxplot_AUPRC_performance_", file.name, ".pdf"))
  print(plot_prc)
  dev.off()
  
  
  
}

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

#Take sensitivities values based on values of specificities
get_sensitivity = function(x, data){
  data %>%
    filter(specificity - x >= 0)%>% #Take specificity values above threshold x
    top_n(sensitivity, n=1) %>% #Take highest sensitivity from that threshold
    mutate(specificity = x, fpr = 1-x) %>% #Define sensitivity based on the specified threshold
    distinct() #If multiple thresholds have same sensitivity values take only one
}

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

calculate_precision <- function(metrics, target) {
  confusion_values <- calculate_confusion_values(metrics, target)
  TP <- confusion_values$TP
  FP <- confusion_values$FP
  
  precision <- TP / (TP + FP)

  return(precision)
}

calculate_recall <- function(metrics, target) {
  confusion_values <- calculate_confusion_values(metrics, target)
  TP <- confusion_values$TP
  FN <- confusion_values$FN
  
  # Calculate recall (sensitivity)
  recall <- TP / (TP + FN)
  
  return(recall)
}

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


######### THIS FUNCTION BELONGS TO MLeval (version 0.3) PACKAGE https://www.rdocumentation.org/packages/MLeval/versions/0.3 (it only has a slight modification to retrieve auc scores)

if(getRversion() >= "2.15.1")  utils::globalVariables(c('Group','PG','PREC','RG','SENS','predmean',
                                                        'realmean'))

#' evalm: Evaluate Machine Learning Models in R
#'
#' evalm is for machine learning model evaluation in R. The function can accept the Caret 'train' 
#' function results to evaluate machine learning predictions or a data frame of probabilities and 
#' ground truth labels can be passed in to evaluate. Probability data must be column1: probability 
#' group1 (column named as your group name 1), column2: probability group2 (column named as your group name 2), 
#' column3: observation labels (column named 'obs'), column4: Group, e.g. different models (column
#' named 'Group'), optional to include if different models are combined horizontally.
#'
#' @param list1 List or data frame: List of Caret results objects from train, or a single train results object, or a data frame of probabilities and observed labels
#' @param gnames Character vector: A vector of group names for the fit objects 
#' @param title Character string: A title for the ROC plot 
#' @param cols Character vector: A vector of colours for the group or groups 
#' @param rlinethick Numerical value: Thickness of the ROC curve line
#' @param fsize Numerical value: Font size for the ROC curve plots
#' @param dlinecol Character string: Colour of the diagonal line
#' @param optimise Character string: Metric by which to select the operating point (INF, MCC, or F1)
#' @param positive Character string: Name of the positive group (will effect PR metrics)
#' @param dlinethick Numerical value: Thickness of the diagonal line
#' @param bins Numerical value: Number of bins for calibration curve
#' @param showplots Logical flag: whether to show plots or not
#' @param plots Character vector: which plots to show: r = roc, pr = proc, prg = precision recall gain, cc = calibration curve 
#' @param percent Numerical value: percentage for the confidence intervals (default = 95)
#' @param silent Logical flag: whether to hide messages (default=FALSE)
#' 
#' @return 
#' List containing: 1) A ggplot2 ROC curve object for printing
#' 2) A ggplot2 PROC object for printing
#' 3) A ggplot2 PRG curve for printing
#' 4) Optimised results according to defined metric
#' 5) P cut-off of 0.5 standard results
#' @export
#'
#' @examples
#' r <- evalm(fit)
evalm <- function(list1,gnames=NULL,title='',cols=NULL,silent=FALSE,
                  rlinethick=1.25,fsize=12.5,dlinecol='grey',
                  dlinethick=0.75,bins=6,optimise='INF',percent=95,
                  positive=NULL){
  
  if(silent==FALSE){
    message('***MLeval: Machine Learning Model Evaluation***')
  }
  
  ## if not a list convert to list
  if (class(list1)[1] != 'list'){
    list1 <- list(list1)
    list1b <- list1
  }
  
  ## read whether input is data frame or caret object
  if (class(list1[[1]])[1] == 'train'){
    input <- 'caret'
    if(silent==FALSE){
      message('Input: caret train function object')
    }
    ## decide whether to average over reps or not
    if (list1[[1]]$control$method == 'cv' | list1[[1]]$control$method == 'LOOCV'){
      mode <- 'cv'
      if(silent==FALSE){
        message('Not averaging probs.')
      }
    }else{
      mode <- 'rep'
      if(silent==FALSE){
        message('Averaging probs.')
      }
    }
  }else if (class(list1[[1]])[1] == 'data.frame'){
    input <- 'normal'
    if(silent==FALSE){
      message('Input: data frame of probabilities of observed labels')
    }
  }else{
    stop('Data frame or Caret train object required please.')
  }
  
  ## if w/o group names, make these
  if (is.null(gnames) & input == 'caret'){
    gnames <- c()
    #
    for (dt in seq(1,length(list1))){
      gnames <- c(gnames,paste('Group',dt))
    }
  }else if (is.null(gnames) & input == 'normal'){
    gnames <- levels(as.factor(list1[[1]]$Group))
  }
  
  ## if not caret get group names here
  if (input == 'normal'){
    df <- list1[[1]]
    names <- colnames(df)[1:2]
  }
  
  ## if not custom colours then use this palette
  if (is.null(cols)){
    cols <- c("red","slateblue","grey","gold","orange","forestgreen", 
              "violetred","violet","skyblue","darkorchid")
  }
  
  ## error handling
  if ( is.null(list1[[1]]$pred) & input == 'caret') {
    stop("No probabilities found in Caret output")
  }
  
  ## set positive names manually if necessary
  if (is.null(positive) == FALSE){
    if (input == 'caret'){
      names <- as.character(list1[[1]]$levels)
      if(which(names==positive)==1){
        namesr <- rev(names)
      }else{
        namesr <- names
      }
    }else if (input == 'normal'){
      names <- names
      if(which(names==positive)==1){
        namesr <- rev(names)
      }else{
        namesr <- names
      }
    }
  }
  
  ## average probabilities if this is a caret object
  if (input == 'caret'){
    
    ## diagnostics for input object
    for (dt in seq(1,length(list1))){
      fit1z <- list1[[dt]]
      if(silent==FALSE){
        message(paste('Group',dt,'type:',fit1z$control$method))
      }
    }
    
    ## for caret fit object input
    names <- as.character(list1[[1]]$levels)
    gnames <- as.factor(gnames)
    
    if (is.null(positive) == FALSE){
      names <- namesr
    }
    G1 <- names[1]
    G2 <- names[2]
    
    ## start caret code (get optimal indices, get mean of repeated cv, combine data)
    list2 <- list()
    ## get caret optimal parameter indices
    for (dt in seq(1,length(list1))){
      fit1z <- list1[[dt]]
      mind <- list()
      for (ii in seq(1,length(fit1z$bestTune))){
        pp <- names(fit1z$bestTune)[ii]
        toadd <- which(fit1z$pred[[pp]] == fit1z$bestTune[[ii]])
        mind[[ii]] <- toadd
      }
      # keep only those indices matching for all parameters
      mind <- Reduce(intersect, mind)
      list2[[dt]] <- mind
    }
    
    ## get data out
    if (mode == 'rep'){
      # start KF
      myl <- list()
      for (dt in seq(1,length(list1))){ # for each fit in list
        ## get data 1
        output <- list1[[dt]]
        indices <- list2[[dt]]
        # code to get a mean series of probabilities out
        preds <- output$pred[indices, ]
        preds <- preds[,c('obs',G1,G2,'rowIndex','Resample')]
        preds$replicates <- sapply(strsplit(preds$Resample,'\\.'), `[`, 2)
        ## resamples is a list of data frames each with an individual sample's probabilities
        ## in for all resamples
        resamples <- split(preds,preds$rowIndex)
        finalres <- matrix(nrow=length(resamples),ncol=4)
        for (i in seq(1,length(resamples))){
          dfx <- resamples[[i]]
          finalres[i,1] <- dfx[1,4] # row index
          finalres[i,2] <- mean(dfx[[G1]]) # mean group 1
          finalres[i,3] <- mean(dfx[[G2]]) # mean group 2
          finalres[i,4] <- as.character(dfx[1,1])
        }
        finalres <- data.frame(finalres)
        finalres$X2 <- as.numeric(as.character(finalres$X2))
        finalres$X3 <- as.numeric(as.character(finalres$X3))
        finalres <- finalres[,c(2,3,4)] # select G1, G2, obs
        colnames(finalres) <- c(G1,G2,'obs')
        finalres$obs <- as.character(finalres$obs)
        finalres$obs <- as.factor(finalres$obs)
        finalres$Group <- gnames[dt]
        #
        myl[[dt]] <- finalres
      }
      # end KF
    }else if (mode == 'cv') {
      # start loocv/ cv
      myl <- list()
      for (dt in seq(1,length(list1))){ # for each fit in list
        ## get data 1
        output <- list1[[dt]]
        indices <- list2[[dt]]
        preds <- output$pred[indices, ]
        finalres <- preds[c(G1,G2,'obs')]
        finalres$obs <- as.character(finalres$obs)
        finalres$obs <- as.factor(finalres$obs)
        finalres$Group <- gnames[dt]
        myl[[dt]] <- finalres
      }
    }
    ### bind them
    fres <- matrix(ncol=4,nrow=0)
    for (dt in seq(1,length(list1))){
      dd <- myl[[dt]]
      fres <- rbind(fres,dd)
    }
    finalres <- fres
    ## end caret code
  }else if (input == 'normal'){
    ## for a non caret input
    # column1: prob G1
    # column2: prob G2
    # column3: obs labels
    # column4: Group (optional)
    # get obs labels
    if (is.null(positive) == FALSE){
      names <- namesr
    }
    G1 <- names[1]
    G2 <- names[2]
    # process
    finalres <- list1[[1]]
    colnames(finalres)[3] <- 'obs'
    finalres$obs <- as.character(finalres$obs)
    finalres$obs <- as.factor(finalres$obs)
    if("Group" %in% colnames(finalres)){
      if(silent==FALSE){
        message('Group column exists.')
      }
    }else{
      finalres$Group <- 'Group1'
      gnames <- factor(c('Group1'))
      if(silent==FALSE){
        message('Group does not exist, making column.')
      }
    }
  }
  
  #dim(finalres)
  #print(finalres$Group)
  #Sys.sleep(100000)
  
  if(silent==FALSE){
    ## sample probabilities diagnostics
    message(paste('Observations:',nrow(finalres)))
    message(paste('Number of groups:',length(gnames)))
    message(paste('Observations per group:',nrow(finalres)/length(gnames)))
    
    ### inform which is positive and negative
    message(paste('Positive:',G2))
    message(paste('Negative:',G1))
  }
  
  #print(finalres$Group)
  gszp<-c()
  gszn<-c()
  ## get N of each group
  for (group in seq(1,length(gnames))){
    #print(gnames[group])
    #print(finalres$Group)
    finalres2 <- subset(finalres, finalres$Group == as.character(gnames[group]))
    #print(head(finalres2))
    if(silent==FALSE){
      message(paste('Group:',as.character(gnames[group])))
      message(paste('Positive:',sum(as.character(finalres2$obs)==G2)))
    }
    gszp[group] <- sum(as.character(finalres2$obs)==G2)
    if(silent==FALSE){
      message(paste('Negative:',sum(as.character(finalres2$obs)==G1)))
    }
    gszn[group] <- sum(as.character(finalres2$obs)==G1)
  }
  
  ## calculate metrics for each threshold do PR and ROC
  finalres$Group <- as.character(finalres$Group)
  aucs <- c()
  aucprs <- c()
  rocm <- NULL # hold the roc curves for each group
  # loop over to get multiple rocs from different groups
  for (group in seq(1,length(gnames))){
    finalres2 <- subset(finalres, finalres$Group == gnames[group])
    ## calculate ml metrics
    mlout <- mlmetrics(finalres2,G1=G1,G2=G2)
    ## append the 0,0 and 1,1 coordinates
    SPECi <- c(0,1-mlout$SPEC)
    SENSi <- c(0,mlout$SENS)
    ## calculate AUC
    auc <- AUCc(SPECi,SENSi,method = c("trapezoid"))
    aucpr <- AUCc(mlout$SENS,mlout$PREC,method = c("trapezoid"))
    ##
    aucs <- c(aucs,round(auc,2))
    aucprs <- c(aucprs,round(aucpr,2))
    rocm <- rbind(rocm, mlout)
  }
  rocm <- data.frame(rocm)
  
  ## reformatting
  rocm[[G1]] <- as.numeric(as.character(rocm[[G1]]))
  rocm[[G2]] <- as.numeric(as.character(rocm[[G2]]))
  # rocm$SENS <- as.numeric(as.character(rocm$SENS))
  # rocm$SPEC <- as.numeric(as.character(rocm$SPEC))
  # rocm$Informedness <- as.numeric(as.character(rocm$Informedness))
  # rocm$PREC <- as.numeric(as.character(rocm$PREC))
  # rocm$NPV <- as.numeric(as.character(rocm$NPV))
  # rocm$F1 <- as.numeric(as.character(rocm$F1))
  rocm$FPR <- 1-rocm$SPEC
  
  ## plot ROC curves on one plot
  gnames <- as.character(gnames)
  llabels <- paste(gnames,c('\n AUC-ROC =','\n AUC-ROC ='),aucs)
  llabels2 <- paste(gnames,c('\n AUC-PR =','\n AUC-PR ='),aucprs)
  
  #
  rocm$Group <- factor(rocm$Group, levels = gnames)
  
  ## prg calculations
  cc<-NULL
  aucprgs <- c()
  for (dt in seq(1,length(gnames))){
    rocmt <- subset(rocm, rocm$Group==gnames[dt])
    # PRG computations 
    # https://papers.nips.cc/paper/5867-precision-recall-gain-curves-pr-analysis-done-right.pdf
    P <- sum(rocmt$obs==G2)
    PN <- length(rocmt$obs)
    bl <- P/PN
    rocmt$PG <- (rocmt$PREC-bl)/(1-bl)*rocmt$PREC
    rocmt$RG <- (rocmt$SENS-bl)/(1-bl)*rocmt$SENS
    rocmt$RG[rocmt$RG<0] <- 0
    rocmt$PG[rocmt$PG<0] <- 0
    # AUC
    aucprg <- AUCc(rocmt$PG,rocmt$RG,method = c("trapezoid"))
    aucprg <- round(aucprg,2)
    aucprgs <- c(aucprgs,aucprg)
    #
    cc<-rbind(cc,rocmt)
  }
  llabels3 <- paste(gnames,c('\n AUC-PRG =','\n AUC-PRG ='),aucprgs)
  rocm<-cc
  
  ## set NA to zero
  rocm$NPV[is.na(rocm$NPV)]<-0
  rocm$PREC[is.na(rocm$PREC)]<-0
  
  #print(rocm)
  #print()
  #Sys.sleep(1000)
  
  ## end metric calcs
  if(silent==FALSE){
    message('***Performance Metrics***')
  }
  
  ## ggplot2 functions
  g <- mroc(rocm,title=title,cols=cols,
            rlinethick=rlinethick,fsize=fsize,dlinecol=dlinecol,
            dlinethick=dlinethick,llabels = llabels) 
  g2 <- proc(rocm,title=title,cols=cols,
             rlinethick=rlinethick,fsize=fsize,dlinecol=dlinecol,
             dlinethick=dlinethick,llabels = llabels2,G1=G2,
             G2=G2) 
  g3 <- prg(rocm,title=title,cols=cols,
            rlinethick=rlinethick,fsize=fsize,
            llabels=llabels3,G1=G2,
            G2=G2)
  g4 <- cc(rocm,title=title,cols=cols,
           rlinethick=rlinethick,fsize=fsize,
           llabels=gnames,G1=G2,
           G2=G2,bins=bins)
  
  if (input == 'caret'){
    ## get group metrics and print
    probl <- list()
    ## split results by group
    for (dt in seq(1,length(list1))){
      probl[[dt]] <- subset(rocm, rocm$Group==gnames[dt])
    }
    names(probl) <- gnames
  }else if (input == 'normal'){ ### probl = problem atm
    ## get group metrics and print
    probl <- list()
    ## split results by group
    for (dt in seq(1,length(gnames))){
      probl[[dt]] <- subset(rocm, rocm$Group==gnames[dt])
    }
    names(probl) <- gnames
  }
  
  ## find and print optimised metrics
  ## extract optimised metrics for results
  ## extract p=0.5 metrics for results also
  optrl <- list()
  stdrl <- list()
  for (dt in seq(1,length(gnames))){
    temp <- probl[[dt]]
    ### get optimised performance (by metric X)
    if (optimise == 'INF'){
      if(silent==FALSE){
        message(paste(gnames[dt],'Optimal Informedness =',temp$Informedness[which.max(temp$Informedness)]))
      }
      optres <- temp[which.max(temp$Informedness),,drop=FALSE]
    }else if (optimise == 'MCC'){
      if(silent==FALSE){
        message(paste(gnames[dt],'Optimal MCC =',temp$MCC[which.max(temp$MCC)]))
      }
      optres <- temp[which.max(temp$MCC),,drop=FALSE]
    }else if (optimise == 'F1'){
      if(silent==FALSE){
        message(paste(gnames[dt],'Optimal F1 score =',temp$F1[which.max(temp$F1)]))
      }
      optres <- temp[which.max(temp$F1),,drop=FALSE]
    }
    optres <- optres[,c('SENS','SPEC','MCC','Informedness','PREC','NPV','FPR','F1',
                        'TP','FP','TN','FN')] # add here for new metrics
    optresb <- optres[,c('SENS','SPEC','MCC','Informedness','PREC','NPV','FPR','F1',
                         'TP','FP','TN','FN')] # add here for new metrics
    optres <- t(optres)
    optres <- data.frame(optres)
    colnames(optres)[1] <- 'Score'
    optres[,1] <- as.numeric(as.character(optres[,1]))
    new <- data.frame(Score=c(aucs[dt],aucprs[dt],aucprgs[dt]))
    ## AUC always go at end
    optres <- rbind(optres,new)
    s<-nrow(optres)-2
    e<-nrow(optres)
    row.names(optres)[s:e] <- c('AUC-ROC','AUC-PR','AUC-PRG') # bump for new metrics +1
    ## get CIs
    optres$CI <- NA
    #print(ciroc)
    optres[s,2] <- ciauc(optres[s,1],gszp[dt],gszn[dt],percent)
    optres['SENS',2] <- wci(optres['SENS',1],gszp[dt],percent)
    optres['SPEC',2] <- wci(optres['SPEC',1],gszn[dt],percent)
    optres['PREC',2] <- wci(optres['PREC',1],(optresb$TP+optresb$FP),percent)
    optres['NPV',2] <- wci(optres['NPV',1],(optresb$TN+optresb$FN),percent)
    #
    optres[,1] <- round(optres[,1],3)
    #
    optrl[[dt]] <- optres
    ### get p=0.5 performance (default)
    stdres <- mlmetricsb(temp,G1=G1,G2=G2)
    stdresb <- stdres
    stdres <- t(stdres)
    stdres <- data.frame(stdres)
    colnames(new) <- 'stdres'
    stdres <- rbind(stdres,new)
    # add CI to standard
    stdres$CI <- NA
    stdres['SENS',2] <- wci(stdres['SENS',1],gszp[dt],percent)
    stdres['SPEC',2] <- wci(stdres['SPEC',1],gszn[dt],percent)
    stdres['PREC',2] <- wci(stdres['PREC',1],(stdresb$TP+stdresb$FP),percent)
    stdres['NPV',2] <- wci(stdres['NPV',1],(stdresb$TN+stdresb$FN),percent)
    stdres[s,2] <- optres[s,2]
    colnames(stdres)[1] <- 'Score'
    row.names(stdres)[s:e] <- c('AUC-ROC','AUC-PR','AUC-PRG')
    stdres[,1] <- as.numeric(as.character(stdres[,1]))
    stdres[,1] <- round(stdres[,1],3)
    stdrl[[dt]] <- stdres
  }
  names(optrl) <- gnames
  names(stdrl) <- gnames
  
  ## reformatting
  for (dt in seq(1,length(gnames))){
    probl[[dt]] <- probl[[dt]][ , -which(names(probl[[dt]]) %in% c("predt"))]
  }
  
  ## print AUCs
  for (n in seq(1,length(gnames))){
    if(silent==FALSE){
      message(paste(gnames[n],'AUC-ROC =',aucs[n]))
    }
  }
  
  ## output
  return(list('roc'=g, 'auc' = aucs))
}

mlmetrics <- function(mlr,G1=G1,G2=G2){
  ## for each group compute...
  p <- mlr
  p <- p[order(-p[,G2]), ]
  i = 1
  for (val in p[[G2]]){
    p$predt <- NA
    p$predt[p[[G2]]<val] <- G1 # for cut-off if less than label as
    p$predt[p[[G2]]>=val] <- G2
    # do calcs
    pred.pos <- p$predt == G2
    pred.neg <- p$predt != G2
    truth.pos <- p$obs == G2
    truth.neg <- p$obs != G2
    p$TP[i] <- sum(pred.pos & truth.pos)
    p$TN[i] <- sum(pred.neg & truth.neg)
    p$FP[i] <- sum(pred.pos & truth.neg)
    p$FN[i] <- sum(pred.neg & truth.pos)
    i = i + 1
  }
  p$SENS <- p$TP/(p$TP + p$FN)
  p$SPEC <- p$TN/(p$TN + p$FP)
  p$Informedness <- p$SENS + p$SPEC - 1
  p$PREC <- p$TP/(p$TP + p$FP)
  p$NPV <- p$TN/(p$TN + p$FN)
  p$MARK <- p$PREC + p$NPV - 1
  p$F1 <- 2*p$PREC*p$SENS/(p$PREC + p$SENS)
  p$F1[is.na(p$F1)] <- 0
  p$MCC <- sign(p$Informedness)*sqrt(p$Informedness*p$MARK) # ((p$TP*p$TN)-(p$FP*p$FN)) / sqrt((p$TP+p$FP)*(p$TP+p$FN)*(p$TN+p$FP)*(p$TN+p$FN))
  # 10 column output
  return(p)
}

mlmetricsb <- function(mlr,G1=G1,G2=G2){
  #print(head(mlr))
  ## for each group compute...
  p <- mlr
  p <- p[order(-p[,G2]), ]
  p$predt[p[[G2]]<0.5] <- G1 # for p=0.5 cut-off
  p$predt[p[[G2]]>=0.5] <- G2
  pred.pos <- p$predt == G2
  pred.neg <- p$predt != G2
  truth.pos <- p$obs == G2
  truth.neg <- p$obs != G2
  TP <- sum(pred.pos & truth.pos)
  TN <- sum(pred.neg & truth.neg)
  FP <- sum(pred.pos & truth.neg)
  FN <- sum(pred.neg & truth.pos)
  SENS <- TP/(TP + FN)
  SPEC <- TN/(TN + FP)
  FPR <- 1-SPEC
  INF <- SENS + SPEC - 1
  PPV <- TP/(TP + FP)
  NPV <- TN/(TN + FN)
  F1 <- 2*PPV*SENS/(PPV + SENS)
  if (is.na(F1)){
    F1 <- 0
  }
  MARK <- PPV + NPV - 1
  MCC <- sign(INF)*sqrt(INF*MARK) # ((TP*TN)-(FP*FN)) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  stdres <- data.frame(SENS=SENS,SPEC=SPEC,MCC=MCC,Informedness=INF,
                       PREC=PPV,NPV=NPV,FPR=FPR,F1=F1,TP=TP,FP=FP,
                       TN=TN,FN=FN) # add here for new metrics must be in order
  # 
  return(stdres)
}

AUCc <- function(x, y, method = "trapezoid"){
  idx <- order(x)
  x <- x[idx]
  y <- y[idx]
  if (method == 'trapezoid'){
    auc <- sum((rowMeans(cbind(y[-length(y)], y[-1]))) * (x[-1] - x[-length(x)]))
  }else if (method == 'step'){
    auc <- sum(y[-length(y)] * (x[-1] - x[-length(x)]))
  }else if (method == 'spline'){
    auc <- integrate(splinefun(x, y, method = "natural"), lower = min(x), upper = max(x))
    auc <- auc$value
  }
  return(auc)
}

mroc <- function(rocm,title='',cols=NULL,
                 rlinethick=1.5,fsize=15,dlinecol='grey',
                 dlinethick=0.75,llabels = llabels){
  
  rocm2 <- NULL
  for (group in unique(rocm$Group)){
    rocmt <- subset(rocm,rocm$Group==group)
    ## add 0,0 co ordinates to ROC and 1,1
    FPR <- rocmt$FPR
    SENS <- rocmt$SENS
    FPR <- c(0,FPR,1)
    SENS <- c(0,SENS,1)
    rocmt <- data.frame(SENS=SENS,FPR=FPR,Group=group)
    rocm2 <- rbind(rocm2,rocmt)
  }
  
  ##
  g <- ggplot(rocm2, aes(x = FPR, y = SENS, color = Group)) + 
    geom_line(size = rlinethick) +
    geom_abline(intercept = 0, slope = 1, colour = dlinecol, linetype = 1,
                size = dlinethick) +
    coord_equal() +
    theme_bw() +
    theme(plot.title = element_text(size = 15, colour = 'black', hjust = 0.5),
          axis.text.y = element_text(size = fsize, colour = 'black'),
          axis.text.x = element_text(size = fsize, colour = 'black', 
                                     angle = 90, vjust = 0.5, hjust=1),
          legend.title=element_blank(),
          legend.text=element_text(size=fsize*0.75),
          legend.key.size = unit(0.5, "lines"),  # Reduce the size of the legend keys
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title.x = element_text(size = fsize, colour = 'black'),
          axis.title.y = element_text(size = fsize, colour = 'black')) +
    xlab('False positive rate') +
    ylab('True positive rate') +
    ggtitle(title) +
    scale_color_manual(values=cols,name = 'Group', labels = llabels)
  return(g)
}

proc <- function(rocm,title='',cols=NULL,
                 rlinethick=1.5,fsize=15,dlinecol='grey',
                 dlinethick=0.75,llabels = llabels,G1=G2,
                 G2=G2){
  ## compute baseline
  P <- sum(rocm$obs==G2)
  PN <- length(rocm$obs)
  bl <- P/PN
  ##
  g <- ggplot(rocm, aes(x = SENS, y = PREC, color = Group)) + 
    geom_line(size = rlinethick) +
    geom_abline(intercept = bl, slope = 0, colour = dlinecol, linetype = 1,
                size = dlinethick) +
    coord_equal() +
    theme_bw() +
    theme(plot.title = element_text(size = fsize, colour = 'black', hjust = 0.5),
          axis.text.y = element_text(size = fsize, colour = 'black'),
          axis.text.x = element_text(size = fsize, colour = 'black', 
                                     angle = 90, vjust = 0.5, hjust=1),
          legend.title=element_blank(),
          legend.text=element_text(size=fsize),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title.x = element_text(size = fsize, colour = 'black'),
          axis.title.y = element_text(size = fsize, colour = 'black')) +
    xlab('Recall(sensitivity)') +    
    ylab('Precision') +
    ggtitle(title) +
    scale_color_manual(values=cols,name = 'Group', labels = llabels) +
    scale_y_continuous(limits = c(0,1))
  return(g)
}

prg <- function(rocm,title='',cols=NULL,
                rlinethick=1.5,fsize=15,
                llabels = llabels,G1=G2,
                G2=G2){
  ## negatives are ignored
  rocm$RG[rocm$RG<0] <- 0
  rocm$PG[rocm$PG<0] <- 0
  ##
  g <- ggplot(rocm, aes(x = RG, y = PG, color = Group)) + 
    geom_line(size = rlinethick) +
    coord_equal() +
    theme_bw() +
    theme(plot.title = element_text(size = fsize, colour = 'black', hjust = 0.5),
          axis.text.y = element_text(size = fsize, colour = 'black'),
          axis.text.x = element_text(size = fsize, colour = 'black', 
                                     angle = 90, vjust = 0.5, hjust=1),
          legend.title=element_blank(),
          legend.text=element_text(size=fsize),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title.x = element_text(size = fsize, colour = 'black'),
          axis.title.y = element_text(size = fsize, colour = 'black')) +
    xlab('Recall gain') +    
    ylab('Precision gain') +
    ggtitle(title) +
    scale_color_manual(values=cols,name = 'Group', labels = llabels) +
    scale_y_continuous(limits = c(0,1)) +
    scale_x_continuous(limits = c(0,1))
  return(g)
}

cc <- function(rocm,title='',cols=NULL,G1=G1,G2=G2,
               rlinethick=1.5,fsize=15,dlinecol='grey',
               dlinethick=0.75,llabels = llabels,
               bins=bins){
  
  ## code to bin pred probabilities and compare with mean
  p2 <- NULL
  for (group in unique(rocm$Group)){
    p <- subset(rocm,rocm$Group==group)
    # process data
    p$obs <- as.character(p$obs)
    p$obs[p$obs!=G2] <- '0'
    p$obs[p$obs==G2] <- '1'
    p <- data.frame(pred=p[[G2]],obs=as.numeric(p$obs))
    # make bins
    p$bin <- cut(p$pred, bins)
    # get bin means
    predmeans <- data.frame(tapply(p$pred, cut(p$pred, bins), mean))
    colnames(predmeans)[1] <- 'predmean'
    realmeans <- data.frame(tapply(p$obs, cut(p$pred, bins), mean))
    colnames(realmeans)[1] <- 'realmean'
    predmeans$bin <- row.names(predmeans)
    realmeans$bin <- row.names(realmeans)
    d1 <- merge(p,predmeans,by='bin')
    d2 <- merge(d1,realmeans,by='bin')
    d2 <- d2[,c('bin','predmean','realmean')]
    deduped.data <- unique( d2[ , 1:3 ] )
    deduped.data$Group <- group
    ##
    p2 <- rbind(p2,deduped.data)
  }
  
  ## plot code
  cc <- ggplot(p2, aes(x = predmean, y = realmean, color = Group)) + 
    geom_abline(intercept = 0, slope = 1,colour = dlinecol, linetype = 1,
                size = dlinethick) +
    geom_line(size = rlinethick) +
    geom_point() +
    scale_y_continuous(limits = c(0,1)) +
    scale_x_continuous(limits = c(0,1)) +
    xlab('Predicted probability') +    
    ylab('True probability in each bin') +
    theme_bw() +
    theme(plot.title = element_text(size = fsize, colour = 'black', hjust = 0.5),
          axis.text.y = element_text(size = fsize, colour = 'black'),
          axis.text.x = element_text(size = fsize, colour = 'black', 
                                     angle = 90, vjust = 0.5, hjust=1),
          legend.title=element_blank(),
          legend.text=element_text(size=fsize),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.title.x = element_text(size = fsize, colour = 'black'),
          axis.title.y = element_text(size = fsize, colour = 'black')) +
    ggtitle(title) +
    scale_color_manual(values=cols,name = 'Group', labels = llabels)
  return(cc)
}

ciauc <- function(auc,N1,N2,ci=95){
  # Hanley, James A., and Barbara J. McNeil. "The meaning and use of the area under a 
  # receiver operating characteristic (ROC) curve." Radiology 143.1 (1982): 29-36.
  # https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/PASS/
  # Confidence_Intervals_for_the_Area_Under_an_ROC_Curve.pdf
  z <- qnorm(1-((1-(ci/100))/2))
  Q1 <- auc/(2-auc)
  Q2 <- (2*auc^2)/(1+auc)
  s1 <- (auc*(1-auc)+(N1-1)*(Q1-auc^2)+(N2-1)*(Q2-auc^2))/(N1*N2)
  se_auc <- sqrt(s1)
  ui <- auc+(z*se_auc)
  li <- auc-(z*se_auc)
  ciauc <- c(round(li,2),round(ui,2))
  ciauc <- paste(ciauc,collapse='-')
  return(ciauc)
}

wci <- function(p,n,ci=95){
  # n = total in fraction e.g. TP+FN or TP+FP
  # p = fraction value in sample, i.e. p^
  # z = z
  # wilson score CI
  # https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/
  # pdf/Procedures/PASS/Confidence_Intervals_for_One-Sample_Sensitivity_and_Specificity.pdf
  z <- qnorm(1-((1-(ci/100))/2))
  num <- (2*n*p+z^2)+z*sqrt((z^2)+4*n*p*(1-p))
  denom <- 2*(n+z^2)
  ui <- num/denom
  num <- (2*n*p+z^2)-z*sqrt((z^2)+4*n*p*(1-p))
  denom <- 2*(n+z^2)
  li <- num/denom
  ciw <- c(round(li,2),round(ui,2))
  ciw <- paste(ciw,collapse='-')
  return(ciw)
}

#' brier_score: A Brier score function
#'
#' Calculates the Brier score to evaluate probabilities. A data frame of probabilities and ground truth labels must
#' be passed in to evaluate. Raw probability data must be column1: prob G1, column2: prob G2,
#' column3: obs labels, column4: Group (optional). Zero is optimal and more positive is less.
#'
#' @param preds Data frame: Data frame of probabilities and ground truth labels.
#' @param positive Character vector: The name of the positive group, must equal a column name consisting of probabilities.

#' @return
#' Brier score
#' @export
#'
#' @examples
#' r2 <- brier_score(preds)
brier_score <- function(preds,positive=colnames(preds)[2]){
  # Zero is optimal and more positive is less
  bs <- mean((preds[[positive]]-as.integer(preds$obs == positive))^2)
  return(bs)
}

#' LL: Log-likelihood function
#'
#' Calculates the Log-likelihood to evaluate probabilities. A data frame of probabilities and ground truth labels must
#' be passed in to evaluate. Raw probability data must be column1: prob G1, column2: prob G2,
#' column3: obs labels, column4: Group (optional). Zero is optimal and more negative is less.
#'
#' @param preds Data frame: Data frame of probabilities and ground truth labels.
#' @param positive Character vector: The name of the positive group, must equal a column name consisting of probabilities.

#' @return
#' Log-likelihood
#' @export
#'
#' @examples
#' r1 <- LL(preds)
LL <- function(preds,positive=colnames(preds)[2]){
  # Zero is optimal and more negative is less
  y <- as.integer(preds$obs == positive) # ground truth
  p <- preds[[positive]]
  L <- sum(log(p*y+(1-y)*(1-p)))
  return(L)
}

