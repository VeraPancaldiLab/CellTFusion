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
