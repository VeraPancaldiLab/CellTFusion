#' Anova test using cell groups scores
#'
#' @param cell.groups A list of cell groups obtained from cell.groups.computation()
#' @param coldata A data frame containing the trait column to test with Anova
#' @param trait Name of the column to test from coldata
#' @param pval pvalue to consider significant (default 0.05)
#'
#' @return A list containing the significant cell groups after Anova test. Additionally it saves in the Results/ folder its respective plots
#' @export
#'
#' @examples
#'
#' res_anova = cell.groups.anova.test(cell.groups, clinical.data, "Response")
#'
cell.groups.anova.test = function(cell.groups, coldata, trait, pval = 0.05){

  sig = c()
  coldata[,trait] = as.factor(coldata[,trait])

  for (j in 1:ncol(cell.groups[[1]])) {
    data = data.frame("Value" = cell.groups[[1]][,j], "Trait" = coldata[,trait])
    model  <- stats::lm(.data$Value ~ .data$Trait, data = data)
    res.aov <- data %>% rstatix::anova_test(.data$Value ~ .data$Trait)

    ##Extract only significant features
    if(round(res.aov$p, 5) <= pval){
      cat("Significant pval after doing Anova test for", colnames(cell.groups[[1]])[j], "\n")
      pdf(paste0("Results/Anova_", trait, "_", colnames(cell.groups[[1]])[j]), width = 12, height = 9)
      print(ggplot(data, aes(x=Trait, y=Value, fill=Trait)) +
              geom_violin(width=0.6) +
              geom_boxplot(width=0.07, color="black", alpha=0.2) +
              scale_fill_brewer() +
              geom_smooth(aes(x=Trait, y=Value), method = "loess") +
              xlab(paste0("Clinical trait: ", trait)) +
              labs(title= paste0("Dendrogram_", colnames(cell.groups[[1]])[j]),
                   subtitle = rstatix::get_test_label(res.aov, detailed = TRUE)) +
              theme(axis.text.x = element_text(angle = 0),
                    axis.title.y = element_text(size = 8, angle = 90)))
      dev.off()

      sig = c(sig, j)
    }
  }

  cell.groups.sig = list()
  cell.groups.sig[[1]] = cell.groups[[1]][,sig]
  cell.groups.sig[[2]] = cell.groups[[2]][sig]

  if(length(sig)==0){
    message("No significant cell groups (pvalue < ", pval, ") after Anova test")
  }else{
    return(cell.groups.sig)
  }

}

#' Fisher test using cell groups scores
#'
#' @param cell.groups A list of cell groups obtained from cell.groups.computation()
#' @param coldata A data frame containing the trait column to test with Fisher
#' @param trait Name of the column to test from coldata
#' @param pval pvalue to consider significant (default 0.05)
#'
#' @return A list containing the significant cell groups after Fisher test. Additionally it saves in the Results/ folder its respective plots
#' @export
#'
#' @examples
#'
#' res_fisher = cell.groups.fisher.test(cell.groups, clinical.data, "Response")
#'
cell.groups.fisher.test = function(cell.groups, coldata, trait, pval = 0.05){

  sig = c()
  coldata[,trait] = as.factor(coldata[,trait])

  for (j in 1:ncol(cell.groups[[1]])) {
    coldata = coldata %>%
      dplyr::mutate(level = .data$cell.groups[[1]][,j],
                    Cells_level = ifelse(level > summary(level)[3], 'High', 'Low'))

    contingency = table(coldata[,"Cells_level"], coldata[,trait])
    test = stats::fisher.test(contingency)

    ##Extract only significant features
    if(round(test$p.value, 5) <= pval){
      cat("Significant pval after doing Fisher test for", colnames(cell.groups[[1]])[j], "\n")
      df = data.frame("Cells_level" = coldata[,"Cells_level"], "Trait" = coldata[,trait])
      pdf(paste0("Results/Fisher_", trait, "_", colnames(cell.groups[[1]])[j]), width = 12, height = 9)
      print(ggbarstats(df, Cells_level, Trait, results.subtitle = F,
                       title= paste0("Dendrogram_", colnames(cell.groups[[1]])[j]),
                       subtitle = paste0("Fisher's exact test, p-value = ", ifelse(test$p.value < 0.001, "< 0.001", round(test$p.value, 5))))+
              ggplot2::theme(plot.title = ggplot2::element_text(size=15), axis.text = ggplot2::element_text(size=14), legend.title = ggplot2::element_text(size=14)))
      dev.off()

      sig = c(j, sig)
    }
  }

  cell.groups.sig = list()
  cell.groups.sig[[1]] = cell.groups[[1]][,sig]
  cell.groups.sig[[2]] = cell.groups[[2]][sig]

  if(length(sig)==0){
    message("No significant cell groups (pvalue < ", pval, ") after Fisher test")
  }else{
    return(cell.groups.sig)
  }

}

#' Extract cells from cell type groups
#'
#' @param groups Cell type groups obtained from cell.groups.computation()
#'
#' @return Cell types composition ignoring method and signature
#' @export
#'
#' @examples
#'
#' cells_names = extract_cells(cell_group)
#'
extract_cells = function(groups){
  names_cells = c("B.cells", "B.naive", "B.memory", "Macrophages.cells", "Macrophages.M0", "Macrophages.M1", "Macrophages.M2", "Monocytes", "Neutrophils", "NK.cells", "NK.activated",
                  "NK.resting", "NKT.cells", "CD4.cells", "CD4.memory.activated", "CD4.memory.resting", "CD4.naive", "CD8.cells", "T.cells.regulatory", "T.cells.non.regulatory","T.cells.helper",
                  "T.cells.gamma.delta", "Dendritic.cells", "Dendritic.activated", "Dendritic.resting", "Cancer", "Endothelial", "Eosinophils", "Plasma.cells", "Myocytes", "Fibroblasts",
                  "Mast.cells", "Mast.activated", "Mast.resting", "CAF")

  regex_pattern <-  paste0("(", paste(names_cells, collapse = "|"), ")")

  extracted_names <- sapply(groups, function(x) {
    match <- regexpr(regex_pattern, x)
    if (match != -1) {
      return(regmatches(x, match))
    } else {
      return(NA)
    }
  })

  extracted_names <- unname(extracted_names)
  extracted_names <- unique(stats::na.omit(extracted_names))
  return(extracted_names)
}

#' Extract colors
#'
#' Extract TF module colors from cell type groups
#'
#' @param module_colors A character vector with the TF module colors. This can be found as an element from the output of compute.WTCNA()
#' @param cell_group_name Cell type group name
#'
#' @return A character vector with the module colors
#' @export
#'
#' @examples
#'
#' color = extract_colors(module_colors, cell.dendrogram.group.1)
#'
extract_colors <- function(module_colors, cell_group_name) {
  module_colors = c(module_colors, "all")
  matches <- c() # For storing the matches
  for (color in module_colors) {
    match <- regexpr(color, cell_group_name) # Find the position of the match
    if (match != -1) {
      matches <- c(matches, regmatches(cell_group_name, match))
    }
  }


  if (length(matches) > 0) {
    order <- sapply(matches, function(m) regexpr(m, cell_group_name)) # Sort matches based on their position in the original string to ensure names are the same
    matches <- matches[order(order)]  # Order the matches based on their position
    return(matches)  # Return the ordered matches
  } else {
    return(NA)  # If no matches are found, return NA
  }

}

#' Create TFs modules
#'
#' This function re-create existing TF modules on a different TF activity matrix.
#'
#' @param TF.matrix TFs activity matrix with samples as rows and TFs as columns.
#' @param network_tfs A TF network object obtained from compute.WTCNA() from which TF modules need to be re-create.
#'
#' @return A matrix of TFs modules scores across samples
#' @export
#'
#' @examples
#'
#' TF.matrix_simulated = create_tfs_modules(TFs_test, network_res)
#'
create_tfs_modules = function(TF.matrix, network_tfs){

  tfs.modules = TF.matrix %>%
    t() %>%
    data.frame() %>%
    mutate(Module = "na")

  for (i in 1:length(network_tfs[[3]])) {
    tfs.modules$Module[which(rownames(tfs.modules) %in% network_tfs[[3]][[i]])] = names(network_tfs[[3]])[i]
  }

  tfs_colors = tfs.modules %>%
    dplyr::pull(Module)

  MEList = WGCNA::moduleEigengenes(TF.matrix, colors = tfs_colors, scale = F) #Data already scale
  MEs = MEList$eigengenes
  MEs = WGCNA::orderMEs(MEs)

  return(MEs)
}

#' Find maximum iteration from subgroups
#'
#' @param cells.groups Cell groups corresponding to a specific cell type.
#'
#' @return Maximum subgroupping iteration
#'
#' @examples
#'
#' iterations = find.maximum.iteration(deconv_subgroups)
#'
find.maximum.iteration = function(cells.groups){
  max_iteration = c()
  for (i in 1:length(cells.groups)){
    if(is.null(names(cells.groups[[i]]))==F){
      iterations <- sapply(names(cells.groups[[i]]), function(x) {
        as.numeric(sub(".*\\.Iteration\\.(\\d+)", "\\1", x))
      })
      local_max = max(iterations)
      max_iteration = c(max_iteration, local_max)
    }
  }

  return(max(max_iteration))
}

#' Merge TFs modules
#'
#' Identify high correlated TFs modules and merge.
#'
#' @param data TFs modules matrix
#' @param colors TFs modules colors
#' @param corr Correlation value above which two modules are merge.
#'
#' @return A list containing
#'
#' - Merge modules
#' - TFs module colors
#'
#'
#' @examples
#'
#' merge = mergeModules(TFs_modules, colors, 0.9)
#'
mergeModules = function(data, colors, corr){
  df = correlation(data)
  idx = which(round(df$r,2) > corr)
  if(length(idx)>0){
    for(i in seq(1, length(idx), by=2)){
      module1 = df$measure1[idx[i]]
      module2 = df$measure2[idx[i]]
      if((module1 %in% colnames(data)) && (module2 %in% colnames(data))){
        colors[which(colors%in%c(substring(module1, 3), substring(module2, 3)))] = substring(module1, 3)
        data <- data %>%
          dplyr::mutate(new_column = rowMeans(dplyr::select(., module1, module2))) %>%
          dplyr::rename(module1 = new_column) %>%
          dplyr::select(., -module1, -module2)
      }
    }
  }

  return(list(data, colors))
}

#' Remove cell groups with equal composition
#'
#' @param cell.values Cell groups scores
#' @param cell.composition Cell groups composition
#'
#' @return A list containing
#'
#' - Cell groups scores after removal of equal cell groups
#' - Cell groups composition after removal of equal cell groups
#'
#' @examples
#'
#' cell.groups = remove_equal(cell.groups.values, cell.groups.composition)
#'
remove_equal = function(cell.values, cell.composition){

  #Sorted list to avoid no recognizing vectors with equal composition but different order of cells
  for(i in 1:length(cell.composition)){
    for (j in 1:length(cell.composition[[i]])) {
      cell.composition[[i]][[j]] = sort(cell.composition[[i]][[j]])
    }
  }

  #Remove cell groups
  for(i in 1:length(cell.composition)){
    exist = c() #Initialize vector of equalities
    rang = seq(1, length(cell.composition))[-i] #Create sequence to iterate all list elements except the one being analyzed
    for (j in rang){
      idx = which(cell.composition[[i]] %in% cell.composition[[j]] == TRUE) #Map all cell groups which already existed
      exist = c(exist, idx) #Save index cluster
    }
    if(length(exist)!=0){
      cell.composition[[i]] = cell.composition[[i]][-unique(exist)] #Remove cell groups that already exist
      cell.values[[i]] = cell.values[[i]][-unique(exist)] #Remove cell groups that already exist
    }
  }

  #Remove dendrograms without elements (length equal 0)
  vec = c()
  for(i in 1:length(cell.composition)){
    if(length(cell.composition[[i]]) == 0){
      vec = c(vec, i)
    }
  }

  if(length(vec)>0){
    cell.composition = cell.composition[-vec]
    cell.values = cell.values[-vec]
  }

  cell.composition = unlist(cell.composition, recursive = FALSE)
  cell.values = unlist(cell.values, recursive = FALSE)

  return(list(cell.values, cell.composition))
}

#' Remove cell groups with only one feature
#'
#' @param cell.values Cell groups scores
#' @param cell.composition Cell groups composition
#'
#' @return A list containing
#'
#' - Cell groups scores after removal of single cell groups
#' - Cell groups composition after removal of single cell groups
#'
#'
#' @examples
#'
#' cell.groups = remove_single_groups(cell.values, cell.composition)
#'
remove_single_groups = function(cell.values, cell.composition){

  message("Removing cell groups composed of one single feature..............................................................................")
  vec = c()
  for (i in 1:length(cell.composition)) {
    if(length(cell.composition[[i]])==1){
      vec = c(vec, i)
    }
  }

  if(length(vec)>0){
    cell.composition = cell.composition[-vec]
    cell.values = cell.values[-vec]
  }

  if(length(cell.composition)==0){
    return(NULL)
  }else{
    return(list(cell.values, cell.composition))
  }

}

#' Calculate dendrogram cuts
#'
#' @param cell.group.dendrogram List with the cell dendrograms corresponding to each TF module obtained from identify.cell.groups()
#' @param n_cuts Optional parameter to limit the number of cuts the dendrogram needs to be cut (Default is NULL). If no parameter is set, number of cuts will be proportional to the height of the dendrogram.
#'
#' @return A list with the sequence of numbers where each dendrogram will be cut
#'
#' @examples
#'
#' cuts = calculate_dendrogram_cuts(cell.dendrograms, n_cuts = NULL)
#' cuts = calculate_dendrogram_cuts(cell.dendrograms, n_cuts = 5)
#'
calculate_dendrogram_cuts = function(cell.group.dendrogram, n_cuts = NULL){

  cuts = list()
  for (i in 1:length(cell.group.dendrogram)) {
    dend_heights <- dendextend::heights_per_k.dendrogram(stats::as.dendrogram(cell.group.dendrogram[[i]])) #Calculate dendrogram heights

    sorted_heights <- sort(dend_heights) # Sort heights
    height_diffs <- diff(sorted_heights) # Calculate differences between consecutive heights

    buffer <- max(stats::median(height_diffs[height_diffs > 0]), 1) # Buffer: take the median of the non-zero differences

    # Add and rest the buffer to the minimum and maximum height respectively to avoid trivial cuts (clusters with 1 feature and cluster with all features)
    min_height <- min(sorted_heights) + buffer
    max_height <- max(sorted_heights) - buffer

    # Ensure the buffer-adjusted min height is valid
    if (min_height >= max_height) {
      stop("Buffered minimum height exceeds maximum height for dendrogram")
    }

    if(is.null(n_cuts)==T){
      number_cuts = floor(max_height) #Truncate based on highest height of dendrogram
    }else{
      number_cuts = n_cuts
    }

    cut_sequence <- seq(min_height, max_height, length.out = number_cuts)  # Generate a sequence of cut heights between the buffered min and max
    #cut_sequence = cut_sequence[-c(1, length(cut_sequence))] #Remove first and last cut to avoid trivial groups (either clusters of 1 feature or a single cluster with all of them)
    cut_sequence <- round(cut_sequence, 2)  # Round the cut heights for cleaner values (2 decimal place)

    cuts[[i]] = cut_sequence
  }

  #Give format to the list to have a sequence of cuts per dendrogram
  if(is.null(n_cuts)==F){
    #Adjust to return the combinations of cuts across all dendrograms
    combined_cuts <- matrix(NA, nrow = n_cuts, ncol = length(cell.group.dendrogram))

    for (i in 1:length(cuts)) {
      combined_cuts[1:length(cuts[[i]]), i] <- cuts[[i]]
    }

    combined_cuts_list <- split(combined_cuts, row(combined_cuts))
    return(combined_cuts_list)
  }else{
    return(cuts)
  }


}

#' Remove high correlated cell groups
#'
#' @param data List with cell groups containing cell groups values and cell groups composition respectively.
#' @param colors TF module colors
#' @param threshold Threshold for defined high correlated features
#'
#' @return A list containing
#'
#' - Cell groups scores after removal of high correlated cell groups
#' - Cell groups composition after removal of high correlated cell groups
#'
#'
#' @examples
#'
#' cell.groups = remove.cell.groups.corr(cell.groups, module_colors, threshold = 0.9)
#'
remove.cell.groups.corr <- function(data, colors, threshold = 0.9) {

  features_high_corr = c()
  # Compute correlation matrix
  corr_matrix <- stats::cor(data[[1]])
  # Find highly correlated features
  contador = 1
  while(nrow(corr_matrix)>0){
    color_features = c()
    feature = data.frame(corr_matrix[1, , drop = FALSE]) #Extract first row feature
    feature = feature %>%                                #Take only high corr above threshold
      dplyr::mutate_all(~ifelse(. > threshold, ., NA)) %>%
      dplyr::select_if(~all(!is.na(.)))

    corr_matrix = corr_matrix[-which(rownames(corr_matrix)%in%colnames(feature)),-which(colnames(corr_matrix)%in%colnames(feature)), drop = F] #Remove already joined features
    color_features = list()
    if(ncol(feature)>1){
      for (m in 1:ncol(feature)) {
        color_group = extract_colors(colors, colnames(feature)[m])
        color_features[[m]] = color_group
      }

      all_equal <- all(sapply(color_features, function(x) identical(x, color_features[[1]]))) #Check whether highly correlated features belong to the same group of TFs modules

      if(all_equal == TRUE){
        new_group_composition = unique(unlist(unname(data[[2]][colnames(feature)])))
        new_group_value = rowMeans(data[[1]][,colnames(data[[1]])%in%colnames(feature)])

        message("Highly correlated features (r>", threshold,"): ", paste(colnames(feature), collapse = ', '), ". Combining.")

        if(contador==1){
          #Remove features from original data
          new_data <- data[[1]][, -which(colnames(data[[1]])%in%colnames(feature)), drop = F]
          new_groups = data[[2]][-which(names(data[[2]]) %in% colnames(feature))]
        }else{
          new_data <- new_data[, -which(colnames(new_data)%in%colnames(feature)), drop = F]
          new_groups = new_groups[-which(names(new_groups) %in% colnames(feature))]
        }

        #Add new combined features
        new_name = paste0("Dendrogram_",  paste0(color_group, collapse = "_"), ".group_combined_", contador)
        new_data = cbind(new_data, new_group_value)
        colnames(new_data)[length(new_data)] = new_name

        new_groups[[length(new_groups)+1]] = new_group_composition
        names(new_groups)[length(new_groups)] = new_name

        contador = contador + 1
      }else{
        message("Highly correlated features between different clusters of TFs modules (r>", threshold,"): ", paste(colnames(feature), collapse = ', '), ". Not joining.")
      }
    }else{
      if(contador == 1){
        new_data = data[[1]]
        new_groups = data[[2]]
      }
    }
  }

  res = list(new_data, new_groups)

  return(res)
}

#' Module enrichment
#'
#' @param tpm.counts A matrix with normalized counts (genes as rows and samples as columns)
#' @param module_color A character vector with TF module colors.
#' @param hub_genes List of hub TFs per module.
#' @param tfs_universe A matrix with TF-gene interactions
#'
#' @return Reactome results
#'
#' @examples
#'
#' reactome = module_enrich(as.matrix(RNA.tpm), color, hub_tfs, interactions)
#'
module_enrich = function(tpm.counts, module_color, hub_genes, tfs_universe){
  # genes = colnames(TFs.matrix)
  # inModule = is.finite(match(module_colors,module))
  # modGenes = genes[inModule]
  targets = tfs_universe[tfs_universe$source %in% hub_genes[[1]][[module_color]],] #Extract targets from TFs
  targets = unique(targets$target) #Keep only unique targets from TFs

  targets_genes = tpm.counts[rownames(tpm.counts)%in%targets,] #Extract gene expression from targets
  targets_genes = targets_genes[order(matrixStats::rowVars(targets_genes), decreasing = T),][1:round(0.2*nrow(targets_genes)),] #Keep only highly variable targets (20%)

  entrz <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = rownames(targets_genes), columns = "ENTREZID", keytype = "SYMBOL") #Change to EntrezID
  universe = AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = rownames(tpm.counts), columns = "ENTREZID", keytype = "SYMBOL") #Change to EntrezID

  reac <- ReactomePA::enrichPathway(gene    = entrz$ENTREZID,
                                    organism     = 'human',
                                    universe = universe$ENTREZID,
                                    pvalueCutoff = 0.05)

  reac@result = reac@result[reac@result$p.adjust<0.05,]

  if(nrow(reac@result)!=0){
    return(reac)
  }

}

#' Compute composite score for cell groups
#'
#' @param cell_group A matrix with the cell deconvolution features from cell group
#' @param color_group A character vector with the TF module group colors corresponding to the cell group (parameter can be obtained with extract_colors())
#' @param tfs.module.matrix A matrix with the TF module matrix. It corresponds to the first element of the output from compute.WTCNA()
#' @param prop_var A numeric value with the minimum variance that should be explained by the PCs (Default is 0.7)
#'
#' @return A numeric vector with the scores across samples
#' @export
#'
#' @examples
#'
#' color = extract_colors(module_colors, names(cell.dendrograms)[1])
#' compute_composite_score(pca_group, color, tfs.modules)
#'
compute_composite_score = function(cell_group, color_group, tfs.module.matrix, prop_var = 0.7){

  module_group = paste0("ME", color_group) #To match with columns of TFs modules
  pca_group <- cbind(scale(tfs.module.matrix[,module_group, drop=F]), scale(cell_group)) #Combined TF module corresponding to each cell group
  svd_result <- svd(pca_group, nu = min(nrow(pca_group), ncol(pca_group)), nv = min(nrow(pca_group), ncol(pca_group))) #Performs Singular Value Decomposition (SVD)
  singular_values <- svd_result$d #Extract singular values
  nPCs <- sum(cumsum(singular_values^2) / sum(singular_values^2) < prop_var) + 1 #Extract number of PCs necessary to cover the prop_var
  variance_explained <- (singular_values^2/sum(singular_values^2)) #Extract the % variance explained by each PC
  weights <- variance_explained[1:nPCs]
  PCs <- svd_result$u #Take components
  aveg_group = rowMeans(pca_group) #Consider average expression to aligned direction of components

  for (i in 1:ncol(PCs)) {
    cor_group = WGCNA::cor(aveg_group, PCs[,i], use = "p")

    # Check if the correlation is finite; if not, set it to 0
    if (!is.finite(cor_group)) {
      cor_group = 0
    }

    if(cor_group<0){
      PCs[,i] = PCs[,i]*-1
    }
  }

  #selected_components = minMax(PCs[,1:nPCs]) # avoiding negatives (refers only to the direction)
  selected_components = PCs[,1] # Taking eigenvalue

  # if(nPCs != 1){
  #   composite_score = rowSums(selected_components * weights)
  # }else{
  #   composite_score = selected_components * weights
  # }

  #composite_score = scale(composite_score)
  composite_score = scale(selected_components)

  return(composite_score)
}

#' Unregister workers
#'
#' @return Clean parallelization
#'
#' @examples
#'
#' unregister_dopar()
#'
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
  gc()
}

