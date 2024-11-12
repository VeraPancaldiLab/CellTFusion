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
