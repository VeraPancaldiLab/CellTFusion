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
