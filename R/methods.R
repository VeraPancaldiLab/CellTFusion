#' Compute one-step CellTFusion
#'
#' @param raw.counts A matrix with the raw counts with samples as columns and genes symbols as rows.
#' @param coldata (Optional) A data frame with the clinical data to use if the user wants to compute the metadata association of the TF modules and the clinical data.
#' @param cbsx.mail (Optional) Credential email for running CibersortX. If not provided, cibersortX method will not be run.
#' @param cbsx.token  (Optional) Credential token for running CibersortX. If not provided, cibersortX method will not be run.
#' @param file_name  (Optional) File name for plots saved in Results/ folder.
#'
#' @return A list containing as elements
#' - The deconvolution matrix (samples as rows, cell deconvolution features as columns)
#' - A matrix with the TFs activity (samples as rows, TFs as columns)
#' - A list with the TFs module network
#' - A matrix with the pathway scores
#' - An object with the processed deconvolution
#' - A list with the cell dendrograms corresponding to each TF-module group
#' - A matrix with the cell groups scores
#' @export
#'
#' @examples
#'
#' res = CellTFusion(raw.counts, cbsx.mail = "XXXXXXX", cbsx.token = "XXXXX", file_name = "Test")
#'

CellTFusion = function(raw.counts, coldata = NULL, cbsx.mail = NULL, cbsx.token = NULL, file_name = NULL){

  #Normalize counts
  counts.norm = data.frame(NormalizeTPM(raw.counts, log = T))
  #Deconvolution
  cat("Calculating cell type deconvolution............................................................\n")
  if(is.null(cbsx.mail)==T || is.null(cbsx.token)==T){
    cat("No CibersortX credentials given, deconvolution will excluded this method.....................\n")
    deconv = compute.deconvolution(raw.counts, normalized = T, methods = c("Quantiseq", "MCP", "xCell", "Epidish", "DeconRNASeq", "DWLS"), file_name = file_name)
  }else{
    deconv = compute.deconvolution(raw.counts, normalized = T, credentials.mail = cbsx.mail, credentials.token = cbsx.token, file_name = file_name)
  }
  #TF activity
  cat("Calculating TF activity............................................................\n")
  tfs = compute.TFs.activity(counts.norm)

  # 1. TFs network construction
  cat("Constructing TF network............................................................\n")
  network = compute.WTCNA(tfs, corr_mod = 0.9, clustering.method = "ward.D2", return = T)
  if(is.null(coldata) == F){
    compute.metada.association(network[[1]], coldata, pval = 0.05, width = 10)
  }
  # 1.2. Modules characterization
  cat("Performing TF module characterization............................................................\n")
  hub_tfs = identify_hub_TFs(t(tfs), network, MM_thresh = 0.8, degree_thresh = 0.9)
  compute.modules.enrichment(counts.norm, hub_tfs)
  # 2. Pathways activity inference
  cat("Calculating pathway activities............................................................\n")
  pathways = compute.pathway.activity(counts.norm)
  compute.modules.relationship(network[[1]], pathways, "Pathways_Progeny-TFs_Modules", width = 15)
  # 3. Deconvolution analysis
  cat("Performing deconvolution analysis............................................................\n")
  dt = compute.deconvolution.analysis(deconv, corr = 0.7, seed = 123)
  compute.modules.relationship(network[[1]], dt[[1]], "Deconvolution-TFs_Modules", vertical = T, height = 30, width = 10, pval = 0.05)
  # 4. Cell groups identification
  cat("Cell groups identification............................................................\n")
  tfs.modules.clusters = compute.TF.network.classification(network, pathways, return = T)
  corr_modules = compute.modules.relationship(network[[1]], dt[[1]], return = T, plot = F)
  cell_dendrograms = identify.cell.groups(corr_modules, tfs.modules.clusters, height = 20, return = T)
  # 5. Cell groups scores
  cat("Computing cell groups scores............................................................\n")
  cell.groups = cell.groups.computation(dt[[1]], tfs.module.network = network, cell.dendrograms = cell_dendrograms)

  res = list("Deconvolution" = deconv, "TFs_matrix" = tfs, "TF_network" = network, "Pathways_scores" = pathways,
             "Processed_deconvolution" = dt, "Cell_dendrograms" = cell_dendrograms, "Cell_groups" = cell.groups)

  return(res)

}

#' Computation of cell groups scores
#'
#' @param deconvolution A dataframe containing the deconvolution features after processing. It corresponds to the first element of the list obtained from compute.deconvolution.analysis()
#' @param tfs.module.network A list with the network information of TF modules obtained from compute.WTCNA()
#' @param cell.dendrograms A list with the cell dendrograms corresponding to each TF module obtained from identify.cell.groups()
#' @param return Boolean value whether to return or not the .csv files with cell groups composition
#'
#' @return A matrix with the cell groups scores across samples
#' @export
#'
#' @examples
#'
#' cell.groups = cell.groups.computation(deconvolution, network, cell_dendrograms)
#'
cell.groups.computation = function(deconvolution, tfs.module.network, cell.dendrograms, return = T){

  cuts = calculate_dendrogram_cuts(cell.dendrograms) #Identify dendrograms cuts
  module_colors = unique(tfs.module.network[[2]]) #Identify module colors

  cell.groups = list()
  for(k in 1:length(cell.dendrograms)){
    groups_cut = list()
    for (i in 1:length(cuts[[k]])){
      clusters <- dendextend::cutree(cell.dendrograms[[k]], h = cuts[[k]][i], order_clusters_as_data=FALSE)
      y = list() #Store cell groups composition
      x = list() #Store cell groups scores
      for (j in 1:length(table(clusters))) {
        cells = names(clusters)[clusters==j]
        y[[j]] = cells
        ###################################################Compute score for each cell group
        if(length(cells)>1){
          pca_group = deconvolution[,colnames(deconvolution) %in% cells, drop = F]
          color = extract_colors(module_colors, names(cell.dendrograms)[k])
          x[[j]] <- compute_composite_score(pca_group, color, tfs.module.network[[1]])
        }else{
          x[[j]] = deconvolution[,colnames(deconvolution) %in% cells]
        }
      }

      ###################################################Remove groups with one feature
      groups_cut[[i]] = remove_single_groups(x, y)
    }

    cell.groups.values <- lapply(groups_cut, function(x) x[[1]])
    cell.groups.names <- lapply(groups_cut, function(x) x[[2]])

    ###################################################Remove groups with equal composition (IMPORTANT TO MAKE IT WITHIN DENDROGRAM - Cell groups can have same composition but belong to a different TF module
    cell.groups[[k]] = remove_equal(cell.groups.values, cell.groups.names)
  }

  ###################################################Naming of cell groups
  for (i in 1:length(cell.groups)) {
    for (j in 1:length(cell.groups[[i]][[1]])) {
      names(cell.groups[[i]][[1]])[j] = paste0("Dendrogram_", names(cell.dendrograms)[i], "_group_", j)
      names(cell.groups[[i]][[2]])[j] = paste0("Dendrogram_", names(cell.dendrograms)[i], "_group_", j)
    }
  }

  cell.groups.values <- lapply(cell.groups, function(x) x[[1]])
  cell.groups.names <- lapply(cell.groups, function(x) x[[2]])

  cell.groups = list(unlist(cell.groups.values, recursive = F), unlist(cell.groups.names, recursive = F))
  cell.groups[[1]] = cell.groups[[1]] %>%
    data.frame() %>%
    mutate("Samples" = rownames(deconvolution)) %>%
    tibble::column_to_rownames("Samples")

  cell.groups = remove.cell.groups.corr(cell.groups, module_colors, threshold = 0.9)

  cat("Removing low variance features (if present)........................\n")
  zero = caret::nearZeroVar(cell.groups[[1]], saveMetrics = TRUE)
  cell.groups[[1]] = cell.groups[[1]][, !zero$nzv]
  cell.groups[[2]] = cell.groups[[2]][!zero$nzv]

  ###################################################Export clusters in a table and save
  clusters = data.frame(matrix(nrow = length(cell.groups[[2]]), ncol = 2))
  for (i in 1:length(cell.groups[[2]])) {
    clusters[i,1] = names(cell.groups[[2]])[[i]]
    clusters[i,2] = paste(cell.groups[[2]][[i]], collapse ="\n")
  }
  colnames(clusters) = c("Cell groups", "Methods-signatures")

  if(return){
    write.csv(clusters, "Results/Cell.groups.composition.csv", row.names = F)
    write.csv(cell.groups[[1]], "Results/Cell.groups.scores.csv")
  }

  return(cell.groups)
}

#' Compute cell groups projection
#'
#' @param deconv_res Output from compute.deconvolution.analysis()
#' @param network_res Output from compute.WTCNA()
#' @param cell_groups Output from cell.groups.computation()
#' @param features Cell groups to compute
#' @param deconvolution_test A matrix of unprocessed deconvolution from independent cohort
#' @param TFs_test A matrix of TF activities scores from independent cohort
#'
#' @return A matrix with the projected cell group scores
#' @export
#'
#' @examples
#'
#' cell.groups.projected = compute_cell_groups_signatures(dt, network, cell.groups, features, deconv_test, tfs_test)
#'
compute_cell_groups_signatures = function(deconv_res, network_res, cell_groups, features, deconvolution_test, TFs_test){

  #Remove colors indicatives from deconvolution features to be able to project them in the raw deconvolution results
  #pattern_colors <- paste0("_(", paste(unique(network_res[[2]]), collapse = "|"), ")$")

  # Use the regular expression to remove the suffix colors from the deconvolution features
  # for (i in 1:length(cell_groups[[2]])) {
  #   cell_groups[[2]][[i]] = gsub(pattern_colors, "", cell_groups[[2]][[i]])
  # }

  ################################################################################Simulate TFs module scores
  TF.matrix_simulated = create_tfs_modules(TFs_test, network_res)
  module_colors = unique(network_res[[2]])

  ################################################################################Simulate cell subgroups
  #Scale deconvolution features by columns for making them comparable between cell types (0-1).
  for (i in 1:ncol(deconvolution_test)) {
    deconvolution_test[,i] = deconvolution_test[,i]/max(deconvolution_test[,i])
  }

  deconv_subgroups <- mapply(c, deconv_res[[3]], deconv_res[[4]], SIMPLIFY = FALSE) #Join cell groups
  iterations = find.maximum.iteration(deconv_subgroups)

  # Create same groups composition
  for (m in 1:iterations) {
    base_groups = list()
    for (i in 1:length(deconv_subgroups)){
      if(length(deconv_subgroups[[i]])!=0){
        idy = grep(paste0("Iteration.",m), names(deconv_subgroups[[i]]))
        if(length(idy)!=0){
          base_groups = append(base_groups, deconv_subgroups[[i]][idy])
        }
      }
    }

    deconv_subgroups_values = c()
    for (i in 1:length(base_groups)) {
      deconv_subgroups_values = cbind(deconv_subgroups_values, Biobase::rowMedians(as.matrix(deconvolution_test[,base_groups[[i]]]))) #Compute median using base groups
    }
    colnames(deconv_subgroups_values) = names(base_groups)
    deconvolution_test = cbind(deconv_subgroups_values, deconvolution_test) # Join cell subgroups and deconv features

  }

  deconvolution_test = deconvolution_test[,colnames(deconvolution_test)%in%colnames(deconv_res[[1]])]

  # Compute composite scores
  idx = which(names(cell_groups[[2]]) %in% features)
  cell_dendrogram = c()
  names = c()
  for (i in 1:length(idx)) {
    pca_cells = deconvolution_test[,cell_groups[[2]][[idx[i]]]]
    pca_cells <- pca_cells[, apply(pca_cells, 2, function(x) all(!is.na(x)) && var(x, na.rm = TRUE) != 0), drop = FALSE] #Drop zero-columns or NAs
    name_cell_group = names(cell_groups[[2]][idx[i]])
    color = extract_colors(module_colors, name_cell_group)

    if(ncol(pca_cells) > 1){ #Check if there are more than 2 columns
      cell_dendrogram = cbind(cell_dendrogram, compute_composite_score(pca_cells, color, TF.matrix_simulated))
      names = c(names, names(cell_groups[[2]])[idx[i]])
    }
  }

  if(is.null(cell_dendrogram)==T){
    print("No composite scores because all features have zero variance.")
  }else{
    colnames(cell_dendrogram) = names
    rownames(cell_dendrogram) = rownames(TF.matrix_simulated)
  }

  return(data.frame(scale(cell_dendrogram)))
}

#' Compute metadata association between TF-modules and clinical data
#'
#' It tests if there is an association between the TF modules scores and the available clinical traits. It does a pearson correlation for quantitative traits and an anova test for the qualitative ones.
#'
#' @param tfs.modules A matrix of TF modules scores across samples obtained from compute.WTCNA()
#' @param coldata A matrix with the clinical data to test
#' @param pval A numeric value to set if a test it is significant or not (default is 0.05)
#' @param width A numeric value to set the width of the labeled heatmap to plot.
#' @param height A numeric value to set the height of the labeled heatmap to plot.
#'
#' @return A labeled heatmap with the pearson correlations and violin plots for the Anova tests are saved in the Results/ directory.
#' @export
#'
#' @examples
#'
#' compute.metada.association(tf_modules, traitData, pval = 0.05)
#'
compute.metada.association = function(tfs.modules, coldata, pval = 0.05, width = 20, height = 8){
  ###Association with categorical variables
  coldata_categorical = coldata %>%
    dplyr::select(dplyr::where(is.character)|dplyr::where(is.factor))

  if(ncol(coldata_categorical)!=0){
    data = cbind(tfs.modules, coldata_categorical)
    pvals = data.frame()
    fvals = data.frame()
    for(i in 1:ncol(tfs.modules)){
      contador = 1
      for (j in (ncol(tfs.modules)+1):ncol(data)) {
        module <- names(data[i])
        trait <- names(data[j])
        avz <- broom::tidy(stats::aov(data[,i] ~ data[,j], data = data))
        pvals[i,contador] = avz$p.value[1]
        fvals[i,contador] = avz$statistic[1]
        contador = contador + 1
        if(avz$p.value[1] < pval) {
          pdf(paste0("Results/ANOVA_", module, "-", trait))
          print(ggplot(data, aes(x=data[,j], y=data[,i], fill=data[,j])) +
                  geom_violin(width=0.6) +
                  geom_boxplot(width=0.07, color="black", alpha=0.2) +
                  scale_fill_brewer() +
                  geom_smooth(aes(x=data[,j], y=data[,i]), method = "loess") +
                  ylab(paste0("Values for ", module)) +
                  xlab(paste0("Clinical trait: ", trait)) +
                  labs(title="One way ANOVA test",
                       subtitle=paste0("F statistic: ", round(avz$statistic[1],3), "\npvalue: ", round(avz$p.value[1], 3))) +
                  theme(axis.text.x = element_text(angle = 0),
                        axis.title.y = element_text(size = 8, angle = 90))+
                  scale_fill_discrete(name = trait))
          dev.off()
        }
      }
    }
    rownames(pvals) = colnames(tfs.modules)
    colnames(pvals) = colnames(coldata_categorical)

    rownames(fvals) = colnames(tfs.modules)
    colnames(fvals) = colnames(coldata_categorical)

    pvals = as.matrix(pvals)
    textMatrix2 = paste("ANOVA\n(", signif(pvals, 2), ")", sep = "")
    dim(textMatrix2) = dim(pvals)
  }

  ###Association with quantitative variables
  coldata_quantitative = coldata %>%
    dplyr::select(dplyr::where(is.numeric))

  if(ncol(coldata_quantitative)!=0){
    moduleTraitCor = WGCNA::cor(tfs.modules, coldata_quantitative, method = "p");
    moduleTraitPvalue = WGCNA::corPvalueStudent(moduleTraitCor, nrow(tfs.modules))

    textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 2), ")", sep = "")
    dim(textMatrix) = dim(moduleTraitCor)

    if(ncol(coldata_categorical)!=0){
      textMatrix = cbind(textMatrix, textMatrix2)
      moduleTraitPvalue = cbind(moduleTraitPvalue, pvals)
      simulated_corr = matrix(stats::runif(n=nrow(pvals)*ncol(pvals), min=-0.1, max=0.1), nrow = nrow(pvals), ncol = ncol(pvals))
      colnames(simulated_corr) = colnames(pvals)
      moduleTraitCor = cbind(moduleTraitCor, simulated_corr)
    }

    idx = which(round(moduleTraitPvalue,2)>pval)
    for (i in idx) {
      textMatrix[i] = NA
    }

    pdf("Results/TF.modules_metadata", width = width, height = height)
    par(mar = c(25, 15, 3, 3))
    WGCNA::labeledHeatmap(Matrix = moduleTraitCor,
                          xLabels = colnames(moduleTraitCor),
                          yLabels = rownames(moduleTraitCor),
                          ySymbols = rownames(moduleTraitCor),
                          colorLabels = FALSE,
                          colors = WGCNA::blueWhiteRed(50),
                          textMatrix = textMatrix,
                          setStdMargins = FALSE,
                          cex.text = 0.5,
                          zlim = c(-1,1),
                          main = paste0("Clinical associations with TFs modules\nOnly showing significant associations (pvalue < ", pval, ")"))
    dev.off()
  }

}

#' Computes TFs enrichment using directed targets for each TFs module.
#'
#'The function computes TF module enrichment using the target genes from each hub TF. Briefly, it extracts the TFs belonging to a specific module, it maps the hub TFs and their corresponding target genes and then performs an over representation analysis (ORA) using Reactome on these genes. To avoid common enrichment results because of the potential overlap of several target genes, a step to filter and only keep unique pathways by module is done.
#' @param RNA.tpm A matrix with the normalized counts (samples as columns and genes symbols as rows)
#' @param hub_tfs List of hub TFs per module, it can be get with the identify_hub_TFs() function.
#'
#' @return
#'
#' Dotplots with enrichment result per module saved in the Results/ directory. If no significant enrichment (pval < 0.05) is found, no plot is saved for this module.
#'
#' @export
#'
#' @examples
#'
#' hub_tfs = identify_hub_TFs(t(tfs), network, MM_thresh = 0.8, degree_thresh = 0.9)
#' compute.modules.enrichment(counts.norm, hub_tfs)
#'
compute.modules.enrichment <- function(RNA.tpm, hub_tfs){
  net = decoupleR::get_collectri(organism = 'human', split_complexes = F) #Get universe
  res = list()
  pathways = list()

  #Pathway enrichment using target genes
  for (i in 1:length(hub_tfs[[1]])) {
    color = names(hub_tfs[[1]])[i]
    res[[i]] = module_enrich(as.matrix(RNA.tpm), color, hub_tfs, net)
    names(res)[i] = color
    pathways[[i]] = res[[i]]@result[["Description"]]
    names(pathways)[i] = color
  }

  ItemsList <- gplots::venn(pathways, show.plot = FALSE)
  x = attributes(ItemsList)

  #Keep only unique pathways
  for (i in 1:length(hub_tfs[[1]])) {
    color = names(hub_tfs[[1]])[i]
    if(is.null(res[[i]]) == T){
      print(paste0("No enrichment for module ", color))
    }else{
      res[[i]]@result = res[[i]]@result[which(res[[i]]@result[["Description"]]%in%x$intersections[[color]]),] #take unique pathways per module
      if(nrow(res[[i]]@result)==0){
        print(paste0("No enrichment for module ", color))
      }else{
        pdf(paste0("Results/Enrichment by Reactome\nModule ", color))
        print(enrichplot::dotplot(res[[i]],  title=paste0("Enrichment by Reactome\nModule ", color)))
        dev.off()
      }
    }
  }

}

#' Compute modules relationship
#'
#' Performs linear correlation between two matrices of features across samples
#'
#' @param tfs_network TFs module matrix (samples x modules).
#' @param matB matrix B to correlate (samples x features).
#' @param file_name string indicating the name of the figure to save.
#' @param width integer indicating the width of the figure, default is 8.
#' @param height integer indicating the height of the figure, default is 8.
#' @param pval p value to use as threshold to differentiate between significant and no significant variables. Default is 0.05.
#' @param padj logical value indicating correction of p-value by Bonferroni method has to be applied. Default is false.
#' @param cor_type type of correlation to be used. Default is pearson “p”.
#' @param return logical value indicating if significant features have to be returned. Default is False.
#' @param vertical logical value indicating if function should return a horizontal or vertical plot (this can be useful when looking at several deconvolution features).
#' @param plot whether to saved or not the plots.
#'
#' @return
#'
#' A heatmap showing the level of correlation between TFs modules and corresponding features. Only significant features are being shown. If return = TRUE it will return a list with two elements: the correlation matrix between the TFs modules and the other features and a character vector containing the names of the significant associated features. Note that features not significantly associated with any module are not returned.
#' @export
#'
#' @examples
#'
#' pathways = compute.pathway.activity(counts.norm)
#' compute.modules.relationship(tfs_modules, pathways, "Pathways_Progeny-TFs_Modules", width = 15)
#' dt = compute.deconvolution.analysis(deconv, corr = 0.7, seed = 123)
#' compute.modules.relationship(tfs_modules, deconvolution_matrix, "Deconvolution-TFs_Modules", vertical = T, height = 30, width = 10, pval = 0.05)
#'
#'
compute.modules.relationship <- function(tfs_network, matB, file_name, width = 8, height = 8, pval=0.05, padj = F, cor_type = "p", return = F, vertical=F, plot = T){

  tfs_network = data.frame(tfs_network)
  matB = data.frame(matB)

  ##check if names from both features are the same
  if(all(rownames(tfs_network)==rownames(matB)) == F){
    stop("No equal names, verify the input objects")
  }


  moduleTraitCor = WGCNA::cor(tfs_network, matB, method = cor_type)
  moduleTraitPvalue = WGCNA::corPvalueStudent(moduleTraitCor, nrow(tfs_network))

  rev = which(colSums(moduleTraitPvalue > pval)==nrow(moduleTraitPvalue)) #check if there are features no significant with any module

  if(length(rev)>0){
    moduleTraitCor = moduleTraitCor[,-rev]
    moduleTraitPvalue = moduleTraitPvalue[,-rev]
  }

  if(padj == T){
    for (i in 1:ncol(moduleTraitPvalue)) {
      moduleTraitPvalue[,i] = stats::p.adjust(moduleTraitPvalue[,i], method = 'bonferroni')
    }
  }

  ##Plot in vertical
  if(vertical == T){
    if(ncol(data.frame(moduleTraitCor))>1){
      #Extract significant features per trait
      sig = list()
      for (i in 1:nrow(moduleTraitPvalue)) {
        sig[[i]] = names(which(signif(moduleTraitPvalue[i,],2)<=pval))
      }
      names(sig) = substring(rownames(moduleTraitPvalue), 3)

      if(return == T){
        retu = list(moduleTraitCor, sig)
        return(retu)
      }else{
        d <- stats::dist(t(moduleTraitCor), method = "manhattan")
        hc1 <- stats::hclust(d, method = "ward.D2")
        vec = hc1[["order"]]
        textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 2), ")", sep = "")
        dim(textMatrix) = dim(moduleTraitCor)
        idx = which(round(moduleTraitPvalue,2)>pval)
        for (i in idx) {
          textMatrix[i] = NA
        }
        textMatrix = t(textMatrix)
        moduleTraitCor = data.frame(t(moduleTraitCor))
        if(plot){
          pdf(paste0("Results/",file_name), width = width, height = height)
          par(mar = c(3, 25, 5, 3))
          WGCNA::labeledHeatmap(Matrix = moduleTraitCor[vec,],
                                xLabels = colnames(moduleTraitCor),
                                yLabels = rownames(moduleTraitCor[vec,]),
                                xLabelsPosition = "top",
                                colors = WGCNA::blueWhiteRed(50),
                                textMatrix = textMatrix[vec,],
                                setStdMargins = F,
                                cex.text = 0.5,
                                zlim = c(-1,1))
          dev.off()
        }
      }}else{
        #Extract significant features per trait
        sig = list()
        for (i in 1:nrow(moduleTraitPvalue)) {
          sig[[i]] = colnames(moduleTraitPvalue)[which(signif(moduleTraitPvalue[i,],2)<=pval)]
        }
        names(sig) = substring(rownames(moduleTraitPvalue), 3)
        if(return == T){
          retu = list(moduleTraitCor, sig)
          return(retu)
        }else{
          textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 2), ")", sep = "")
          idx = which(round(moduleTraitPvalue,2)>pval)
          for (i in idx) {
            textMatrix[i] = NA
          }
          textMatrix = t(textMatrix)
          moduleTraitCor = data.frame(t(moduleTraitCor))
          colnames(moduleTraitCor)[1] = colnames(matB)
          if(plot){
            pdf(paste0("Results/",file_name), width = width, height = height)
            par(mar = c(25, 15, 3, 3))
            WGCNA::labeledHeatmap(Matrix = moduleTraitCor,
                                  xLabels = colnames(moduleTraitCor),
                                  yLabels = rownames(moduleTraitCor),
                                  xLabelsPosition = "top",
                                  colors = WGCNA::blueWhiteRed(50),
                                  textMatrix = textMatrix,
                                  setStdMargins = F,
                                  cex.text = 0.5,
                                  zlim = c(-1,1))
            dev.off()
          }
        }}}


  ###Plot in horizontal
  if(vertical == F){
    if(ncol(data.frame(moduleTraitCor))>1){
      #Extract significant features per trait
      sig = list()
      for (i in 1:nrow(moduleTraitPvalue)) {
        sig[[i]] = names(which(signif(moduleTraitPvalue[i,],2)<=pval))
      }
      names(sig) = substring(rownames(moduleTraitPvalue), 3)
      if(return==T){
        retu = list(moduleTraitCor, sig)
        return(retu)
      }else{
        d <- stats::dist(t(moduleTraitCor), method = "manhattan")
        hc1 <- stats::hclust(d, method = "ward.D2")
        vec = hc1[["order"]]
        textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 2), ")", sep = "")
        dim(textMatrix) = dim(moduleTraitCor)
        idx = which(round(moduleTraitPvalue,2)>pval)
        for (i in idx) {
          textMatrix[i] = NA
        }
        moduleTraitCor = data.frame(moduleTraitCor)
        if(plot){
          pdf(paste0("Results/",file_name), width = width, height = height)
          par(mar = c(25, 15, 3, 3))
          WGCNA::labeledHeatmap(Matrix = moduleTraitCor[,vec],
                                xLabels = names(moduleTraitCor[,vec]),
                                yLabels = rownames(moduleTraitCor),
                                ySymbols = rownames(moduleTraitCor),
                                colorLabels = FALSE,
                                colors = WGCNA::blueWhiteRed(50),
                                textMatrix = textMatrix[,vec],
                                setStdMargins = FALSE,
                                cex.text = 0.5,
                                zlim = c(-1,1),
                                main = paste("Module-trait relationships"))
          dev.off()
        }
      }}else{
        #Extract significant features per trait
        sig = list()
        for (i in 1:nrow(moduleTraitPvalue)) {
          sig[[i]] = colnames(moduleTraitPvalue)[which(signif(moduleTraitPvalue[i,],2)<=pval)]
        }
        names(sig) = substring(rownames(moduleTraitPvalue), 3)

        if(return == T){
          retu = list(moduleTraitCor, sig)
          return(retu)
        }else{
          textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 2), ")", sep = "")
          idx = which(round(moduleTraitPvalue,2)>pval)
          for (i in idx) {
            textMatrix[i] = NA
          }
          moduleTraitCor = data.frame(moduleTraitCor)
          colnames(moduleTraitCor)[1] = colnames(matB)
          if(plot){
            pdf(paste0("Results/",file_name), width = width, height = height)
            par(mar = c(25, 15, 3, 3))
            WGCNA::labeledHeatmap(Matrix = moduleTraitCor,
                                  xLabels = names(moduleTraitCor),
                                  yLabels = rownames(moduleTraitCor),
                                  ySymbols = rownames(moduleTraitCor),
                                  colorLabels = FALSE,
                                  colors = WGCNA::blueWhiteRed(50),
                                  textMatrix = textMatrix,
                                  setStdMargins = FALSE,
                                  cex.text = 0.5,
                                  zlim = c(-1,1),
                                  main = paste("Module-trait relationships"))
            dev.off()
          }
        }}}
}

#' Computes TF-modules pathway activities
#'
#'Pathways activity computation is done by using a multivariate linear model (mlm) and a resource that leverages a large compendium of publicly available signaling perturbation experiments to yield a common core of pathway responsive genes for human PROGENy (Schubert M et al, 2018). Users can also select to do an additional Gene Set Variation Analysis (GSVA) using hallmark signatures or predefined gene sets.
#'
#' @param RNA.tpm gene expression matrix normalized with genes as rows and samples as columns.
#' @param gene_sets either hallmark signatures extracted from the Molecular Signature Database (MSigDB) or predefined gene sets signatures. Default is NULL.
#' @param paths whether external database wants to be provided for computing the pathway activities
#'
#' @return A matrix with PRGENy pathway activities (samples as rows and pathways as columns). If gene_sets is not NULL it will also return a matrix containing the signatures scores (samples as rows and signatures as columns).
#'
#' @export
#'
#' @references
#'
#' Schubert M, Klinger B, Klünemann M, Sieber A, Uhlitz F, Sauer S, Garnett MJ, Blüthgen N, Saez-Rodriguez J. 2018. Perturbation-response genes reveal signaling footprints in cancer gene expression. Nature Communications: 10.1038/s41467-017-02391-6.
#'
#'
#' @examples
#'
#' pathways = compute.pathway.activity(counts.norm)
#'
compute.pathway.activity <- function(RNA.tpm, gene_sets = NULL, paths = NULL){

  RNA.tpm = as.matrix(RNA.tpm)
  #Get universe
  if(is.null(paths)){
    paths <- decoupleR::get_progeny(organism = 'human', top = 500)
  }

  # Run mlm
  progeny <- decoupleR::run_mlm(mat=RNA.tpm, net=paths, .source='source', .target='target', .mor='weight', minsize = 5)

  #Remove variable
  rm(paths)
  gc()

  # Transform to wide matrix
  sample_acts_progeny <- progeny %>%
    tidyr::pivot_wider(id_cols = 'condition', names_from = 'source',
                       values_from = 'score') %>%
    tibble::column_to_rownames('condition') %>%
    as.matrix()

  if(is.null(gene_sets)==F){

    cat("Computing GSVA analysis using provided gene sets.....................................................\n")

    gsva_results <- GSVA::gsva(
      RNA.tpm,
      gene_sets,
      method = "gsva",
      kcdf = "Gaussian",
      min.sz = 1,
      mx.diff = TRUE,
      verbose = TRUE
    )
    sample_acts_hallmarks <- data.frame(scale(t(gsva_results)))
    return(list(sample_acts_progeny, sample_acts_hallmarks))
  }

  # Scale per feature
  sample_acts_progeny <- data.frame(scale(sample_acts_progeny))

  message("Pathways scores computed")
  return(sample_acts_progeny)

}

#' Compute survival analysis
#'
#' Computes univariate and multivariate cox proportional hazards (coxPH) models and Kaplan-Meier curves to evaluate across cell type groups. After ﬁtting the CoxPH models to different cell type groups combinations, patient are stratiﬁed based on the linear predictors of the model (risk scores) from which we deﬁne as ‘high’ the patients with risk scores above the median value of the cox model’s linear predictors and as ‘low’ the patients below it. It then performes a Kaplan Meier analysis and plotted the survival curves for each risk group. Finally both survival curves are assessed via a log rank test.
#'
#' @param features A matrix with the cell groups scores
#' @param survival.data Survival information containing progression free survival information (named as 'PFS') and the binary value indicating death/progression/recurrence (named as 'DRP_st')
#' @param time_unit Time unit for variables (either days/months/years)
#' @param p.value pvalue for significance. Default is 0.05
#' @param thres Threshold for classifying the risk scores. Default is 0.5 (median).
#' @param max_factors Maximum number of covariables in the lineal model to consider. Default is Inf, meaning it will test all possible combinations and keep only the significant ones.
#'
#' @return A list with the significant combinations formulas along with their survival curves saved in the Results/ directory.
#'
#' @export
#'
#' @examples
#'
#' survival_groups = compute.survival.analysis(features = cell_groups, survival_data, time_unit = "days", p.value = 0.01, max_factors = 3)
#'
compute.survival.analysis = function(features, survival.data, time_unit, p.value = 0.05, thres = 0.5, max_factors = Inf) {
  n_features <- ncol(features)
  significant_combinations <- list() # To store significant feature combinations

  # Generate all possible combinations of the features
  contador = 1
  for (n in 1:min(n_features, max_factors)) {
    combinations <- utils::combn(1:n_features, n, simplify = FALSE)

    for (comb in combinations) {
      # Create a formula dynamically based on the combination
      formula <- as.formula(paste("Surv(time, status) ~", paste(colnames(features)[comb], collapse = " + ")))

      # Prepare the data for survival analysis
      data_for_model <- data.frame("time" = survival.data$PFS,
                                   "status" = survival.data$DRP_st)

      data_for_model = cbind(data_for_model, features[,comb, drop=F])

      # Fit the Cox PH model with the combination of features (cox PH take into account covariates and measure the impact of each variable in the survival time)
      cox <- rms::cph(formula, data = data_for_model)
      data_for_model$CoxPredictors <- cox$linear.predictors #linear predictors is the risk score for each individual in the dataset

      # Check that the model is significant as a predictor (maybe not useful?, it gives the same linear.predictos - to be check)
      cphmodel <- survival::coxph(survival::Surv(time, status) ~ CoxPredictors, data = data_for_model)
      data_for_model$CoxPredictors <- cphmodel$linear.predictors

      quantiles <- stats::quantile(data_for_model$CoxPredictors, thres)

      # Binarize the Cox model output to draw two KM lines (linear predictors are used to stratify between high-risk and low-risk groups)
      data_for_model$coxHL <- ifelse(cphmodel$linear.predictors >= quantiles, 'High', "Low")

      # Perform Kaplan-Meier based on coxHL
      km_fit <- survival::survfit(survival::Surv(time, status) ~ coxHL, data = data_for_model)

      pval <- survminer::surv_pvalue(km_fit, data = data_for_model)$pval #Performs log-rank test to see whether both survival curves are significantly different

      if (!is.na(pval) && pval < p.value) {
        significant_combinations[[contador]] <- formula
        names(significant_combinations)[contador] = paste0("Formula_", contador)

        pdf(paste0("Results/SurvPlot_", names(significant_combinations)[contador]), width = 10, height = 5, onefile = FALSE)
        print(survminer::ggsurvplot(km_fit,
                                    data = data_for_model,
                                    size = 1,
                                    palette = c("#E7B800", "#2E9FDF"),
                                    conf.int.style = "step",
                                    pval = TRUE,
                                    risk.table = TRUE,
                                    risk.table.col = "strata",
                                    legend.labs = c("High", "Low"),
                                    risk.table.height = 0.3,
                                    ggtheme = theme_grey(),
                                    title = paste0("Cox PH for ", names(significant_combinations)[contador]),
                                    xlab = paste0("Time to death/recurrence/progression (", time_unit, ")")
        ))
        dev.off()
        contador = contador + 1
      }
    }
  }

  if (length(significant_combinations) == 0) {
    print("No significant combinations found.")
  } else {
    return(significant_combinations)
  }

}

#' Compute TF network classification
#'
#' Use the pathways information to classify TF modules into different clusters based on their association values across samples.
#'
#' @param tf.network List with the TF network information obtained from compute.WTCNA()
#' @param pathways.features A matrix with the pathway activities obtained from compute.pathway.activity()
#' @param return Whether to return the intermediate plots generated
#'
#' @return A list containing the clusters of TF modules
#' @export
#'
#' @examples
#'
#' tfs.modules.clusters = compute.TF.network.classification(network, pathways, return = T)
#'
compute.TF.network.classification = function(tf.network, pathways.features, return = T){

  tf.network = data.frame(tf.network[[1]])
  pathways.features = data.frame(pathways.features)

  moduleTraitCor = WGCNA::cor(tf.network, pathways.features, method = "p")
  # moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(tf.network))
  # significance_threshold <- 0.05
  #
  # # Replace non-significant correlations with zero
  # moduleTraitCor[moduleTraitPvalue > significance_threshold] <- 0
  #
  # # Identify columns that have variance (i.e., are not constant or all zeros)
  # non_constant_columns <- apply(moduleTraitCor, 2, function(x) var(x) != 0)
  # moduleTraitCor <- moduleTraitCor[, non_constant_columns] # Filter out the columns that are constant or all zeros

  ### Find clusters
  silhouette = factoextra::fviz_nbclust(moduleTraitCor, hcut, method = "silhouette", k.max = nrow(moduleTraitCor)-1)
  k_cluster = as.numeric(silhouette$data$clusters[which.max(silhouette$data$y)])

  if(return){
    pdf(paste0("Results/TFs_modules_Silhouette_scores"))
    print(silhouette)
    dev.off()
  }

  hc_modules = stats::hclust(stats::dist(moduleTraitCor), method = "ward.D2")
  dend_pathways = stats::as.dendrogram(hc_modules)

  if(return){
    pdf(paste0("Results/TFs_modules_clusters"))
    plot(dend_pathways, cex = 0.6)
    stats::rect.hclust(hc_modules, k = k_cluster, border = 2:5)
    dev.off()
  }

  ### Extract clusters
  sub_grp <- dendextend::cutree(hc_modules, k = k_cluster)

  ### Plot PCA and biplot
  p = factoextra::fviz_cluster(list(data = moduleTraitCor, cluster = sub_grp))

  if(return){
    pdf(paste0("Results/PCA_TFs_modules_clusters"))
    print(p)
    dev.off()
  }

  res.pca <- stats::prcomp(moduleTraitCor,  scale = F)
  p = factoextra::fviz_pca_biplot(res.pca, label="all", select.var = list(contrib = 6), addEllipses=TRUE, ellipse.level=0.75)

  if(return){
    pdf(paste0("Results/PCA_Biplot_TFs_modules_clusters"))
    print(p)
    dev.off()
  }

  # Extract the loadings
  loadings <- res.pca$rotation
  contribution <- (loadings^2)*100
  features = contribution %>%
    data.frame() %>%
    tibble::rownames_to_column("Features") %>%
    dplyr::arrange(desc(PC1))

  if(return){
    pdf("Results/Pathways_contribution_pca.pdf", width = 12, height = 8)
    p = ggplot(features, aes(x = reorder(Features, -PC1), y = PC1)) +
      geom_bar(stat = "identity", fill = "skyblue") +
      labs(title = "Contribution of pathways",
           x = "Feature",
           y = "Contribution (%)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size=15))
    print(p)
    dev.off()
  }

  groups = list()
  for (i in 1:length(unique(sub_grp))) {
    groups[[i]] = names(sub_grp)[sub_grp == i]
    names(groups)[i] = paste0(gsub("ME", "", groups[[i]]), collapse = "_")
  }

  return(groups)
}

#' Compute TF activity
#'
#' The function infers TFs activity based on a gene expression matrix using the VIPER algorithm from Alvarez et al. and a collection of TF-genes interactions (CollecTRI or Dorothea) from the OmnipathR package (Türei D et al. 2016). Users can also input TFs-gene interactions from ARACNE [4,5].
#'
#' @param RNA.counts Gene expression matrix with genes as rows and samples as columns, normalized either by TPM or other type of normalization.
#' @param TF.collection TF-gene interactions collection matrix used to infer the TFs activity. Default argument is “CollecTRI”, other arguments accepted are “Dorothea” (Garcia-Alonso L et al. 2019) or “ARACNE” (Lachmann A et al. 2016, Margolin AA et al. 2006) for input your own TF-gene collection.
#'
#' If TF.collection = “ARACNE”, the user will need to specify the path where the network file is located. This should be a 3 columns text file, with the “regulator” in the first column, “target” in the second and “mutual information” in the third one. If TF.collection = “Dorothea” the list of confidence levels to return from regulons when selecting will be A, B, C and D.
#'
#' @param min_targets_size Integer indicating the minimum number of targets allowed per regulon. Default value is 5.
#' @param tfs.pruned Logical value indicating if TFs need to be pruned, meaning, limiting the maximum size of the regulons. Default is FALSE.
#'
#' Pruned the regulons will help to deal with statistical bias produced by TFs that regulate hundreds of genes compared to others that regulate only a few ones. This argument must be carefully chosen as it will impact TFs scores and module formation.
#'
#' @return A scaled matrix of inferred protein activity with samples as rows and TFs as columns.
#' @export
#'
#' @references
#' Alvarez, M., Shen, Y., Giorgi, F. et al. Functional characterization of somatic mutations in cancer using network-based inference of protein activity. Nat Genet 48, 838–847 (2016). https://doi.org/10.1038/ng.3593
#'
#' Türei D, Korcsmáros T, Saez-Rodriguez J. OmniPath: guidelines and gateway for literature-curated signaling pathway resources. Nat Methods. 2016 Nov 29;13(12):966-967. doi: 10.1038/nmeth.4077. PMID: 27898060.
#'
#' Garcia-Alonso L, Holland CH, Ibrahim MM, Turei D, Saez-Rodriguez J. Benchmark and integration of resources for the estimation of human transcription factor activities. Genome Research. 2019. DOI: 10.1101/gr.240663.118.
#'
#' Lachmann A, Giorgi FM, Lopez G, Califano A. ARACNe-AP: gene network reverse engineering through adaptive partitioning inference of mutual information. Bioinformatics. 2016 Jul 15;32(14):2233-5. doi: 10.1093/bioinformatics/btw216. Epub 2016 Apr 23.
#'
#' Margolin AA, Nemenman I, Basso K, Wiggins C, Stolovitzky G, Dalla Favera R, Califano A. ARACNE: an algorithm for the reconstruction of gene regulatory networks in a mammalian cellular context. BMC Bioinformatics. 2006 Mar 20;7 Suppl 1:S7. doi: 10.1186/1471-2105-7-S1-S7.
#'
#' @examples
#'
#' tfs = compute.TFs.activity(counts.normalized)
#'
#'

compute.TFs.activity <- function(RNA.counts, TF.collection = "CollecTRI", min_targets_size = 5, tfs.pruned = FALSE){

  tfs2viper_regulons <- function(df){
    regulon_list <- split(df, df$source)
    regulons <- lapply(regulon_list, function(regulon) {
      tfmode <- stats::setNames(regulon$mor, regulon$target)
      list(tfmode = tfmode, likelihood = rep(1, length(tfmode)))
    })
    return(regulons)}

  if(tfs.pruned==T){
    cat("Pruned TFs is set to TRUE. Please specify the maximun size of targets allowed/n")
    max_size_targets = as.numeric(readline(prompt = "Maximun size of TFs-targets: "))
  }

  if(TF.collection == "CollecTRI"){
    net = decoupleR::get_collectri(organism = 'human', split_complexes = F)
    net_regulons = tfs2viper_regulons(net)
  } else if(TF.collection == "Dorothea"){
    net = decoupleR::get_dorothea(organism = 'human', levels = c("A", "B", "C", "D"))
    net_regulons = tfs2viper_regulons(net)
  }

  if(TF.collection == "ARACNE"){
    cat("For ARACNE analysis you need to specify the path of your network file. Remember this file should be a 3 columns text file, with regulator in the first column, target in the second and mutual information in the third column")
    network_file = readline(prompt = "Path for network file from aracne (no quotes): ")
    net_regulons <- viper::aracne2regulon(network_file, as.matrix(RNA.counts), format = "3col")
  }

  if(tfs.pruned == TRUE){
    net_regulons = viper::pruneRegulon(net_regulons, cutoff = max_size_targets)
  }

  sample_acts <- viper::viper(as.matrix(RNA.counts), net_regulons, minsize = min_targets_size, verbose=F, method = "scale")
  message("TFs scores computed")

  return(data.frame(t(sample_acts)))

}

#' Compute Weighted TF-coactivity network analysis (WTCNA)
#'
#' Construct a weighted signed or unsigned network using TFs activity to cluster protein regulators into modules that share similar activity patterns based on the expression of their target genes. Each TFs module will have a score per sample represented by the eigenvalue of the module.
#'
#' @param TFs.matrix TFs activity matrix with samples as rows and TFs as columns.
#' @param network.type Network type. Allowed values are “signed”, “unsigned”, “signed hybrid”, “distance”. Default value is “signed”.
#' @param clustering.method Character string specifying the function to be used to calculate co-expression similarity for distance networks. Defaults to the function dist. Default method is “ward.D2”.
#' @param minMod Integer indicating the minimum number of TFs allowed for each module. Default is 30.
#'
#' For more information about these parameters we invited the user to read the documentation in Langfelder, P et al, 2008. minMod must be carefully used as it will impact how many modules the user will have. Default value of 30 might hide modules composed of just a few TFs (around 10) that are giving an explanation of your dataset, so we invited the user to explore different options of values for this parameter.
#'
#' @param corr_mod Value from 0 to 1 used for merge modules in a second iteration. After the first TF-module construction, a correlation will take place between modules and modules correlated by a value higher than this threshold will be merged. Default value is 0.8.
#' @param cor_type Type of correlation to be used for calculation of adjacency matrix and merging modules. Default is pearson “p”, other alternatives are spearman “s”.
#' @param return Whether to save the plots or not in Results/ directory.
#'
#' @return A list containing two elements:
#' - Matrix of TFs-modules with samples as rows and eigenvalues of TF-modules as columns.
#' - Colors of TFs modules.
#' - Named list of TFs and their respective module color assignment.
#' - Proportion of variance explained by each eigenvalue of each TFs module.
#'
#' @export
#'
#' @references
#'
#' Langfelder, P., Horvath, S. WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 9, 559 (2008). https://doi.org/10.1186/1471-2105-9-559
#'
#' @examples
#'
#' network = compute.WTCNA(tfs, corr_mod = 0.9, clustering.method = "ward.D2", return = T)
#'
compute.WTCNA <- function(TFs.matrix, network.type = "unsigned", clustering.method = "ward.D2", minMod = 15, corr_mod = 0.9, cor_type = "p", return = T){

  cat("Creating weighted TF-coactivity network......................................................................\n\n")
  #####Choose parameter for scale-free network topology
  powers = c(c(1:10), seq(from = 12, to=20, by=1))
  sft = WGCNA::pickSoftThreshold(TFs.matrix, powerVector = powers, verbose = 0, networkType = network.type)

  if(return){
    pdf("Results/Soft Threshold")
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
         main = paste("Scale independence"))
    graphics::text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=0.9,col="red");
    graphics::abline(h=0.90,col="red")
    dev.off()
  }

  #Automatic choosing of soft-threshold
  target = 0.9 #Target SFT.R.sq value
  diff = abs(-sign(sft$fitIndices[,3])*sft$fitIndices[,2] - target) #Calculate absolute difference
  min_index = which.min(diff) #Identify the index with the minimum difference
  softPower = powers[min_index]
  cat("Choosing", softPower, "as soft-threshold......................................................................\n\n")

  #####Co-expression matrix using nodes adjacency and topological overlapping nodes
  cat("Calculating nodes adjacency and topological overlapping nodes.................................................\n\n")
  adjacency = WGCNA::adjacency(TFs.matrix, power =softPower, type=network.type, corFnc = "cor", corOptions = list(use = cor_type))
  TOM = WGCNA::TOMsimilarity(adjacency, TOMType = network.type)
  dissTOM = 1-TOM

  #####Unsupervised hierarchical clustering using dissimilarity matrix
  geneTree = stats::hclust(stats::dist(dissTOM), method = clustering.method)
  dynamicMods = dynamicTreeCut::cutreeDynamic(dendro = geneTree, distM = dissTOM,
                                              deepSplit = 2, pamRespectsDendro = FALSE,
                                              minClusterSize = minMod);
  dynamicColors = WGCNA::labels2colors(dynamicMods)

  #####Remove variables and clean garbage
  rm(adjacency, TOM, dissTOM)
  gc()

  if(return){
    pdf("Results/Gene dendrogram and module colors")
    WGCNA::plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                               dendroLabels = FALSE, hang = 0.03,
                               addGuide = TRUE, guideHang = 0.05,
                               main = "Gene dendrogram and module colors")
    dev.off()
  }

  #####Calculate eigenvectors from modules
  cat("Calculating eigenvectors from modules.................................................\n\n")
  MEList = WGCNA::moduleEigengenes(TFs.matrix, colors = dynamicColors, scale = F) #Data already scale
  MEs = MEList$eigengenes
  MEs = WGCNA::orderMEs(MEs)

  print(paste0("Merging modules significantly correlated with ", corr_mod, "........"))
  merge = mergeModules(MEs, dynamicColors, corr_mod)
  MEs = merge[[1]]
  dynamicColors = merge[[2]]
  modtfs = list()
  for(i in 1:ncol(MEs)){
    tfs = colnames(TFs.matrix)
    modules = c(substring(names(MEs), 3))
    inModule = is.finite(match(dynamicColors,modules[i]))
    modtfs[[i]] = tfs[inModule]
  }
  names(modtfs) = modules

  if(return){
    pdf("Results/Gene dendrogram and module colors after merging")
    WGCNA::plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                               dendroLabels = FALSE, hang = 0.03,
                               addGuide = TRUE, guideHang = 0.05,
                               main = "Gene dendrogram and module colors")
    dev.off()
  }

  TFspropVar = WGCNA::propVarExplained(TFs.matrix, dynamicColors, MEs, corFnc = "cor", corOptions = "use = 'p'")

  output = list(MEs, dynamicColors, modtfs, TFspropVar)
  names(output) = c("TFs module matrix", "TFs colors", "TFs per module", "Proportion of variance")

  contador = 1
  tfs_modules = data.frame(matrix(nrow = length(modtfs), ncol = 2))
  colnames(tfs_modules) = c("TFs module", "Composition")
  for (i in 1:length(modtfs)) {
    tfs_modules[contador,1] = names(modtfs)[i]
    tfs_modules[contador,2] = paste(modtfs[[i]], collapse = ",")
    contador = contador + 1
  }

  if(return){
    write.csv(tfs_modules, 'Results/TFs_modules.csv', row.names = F)
  }

  return(output)

}

#' Identify hub TFs
#'
#' Identifies hub TFs per module using values of module membership and degree.
#'
#' @param datExpr A matrix of TF activity (TFs as rows and samples as columns).
#' @param TF.network TF network obtained from compute.WTCNA().
#' @param MM_thresh Threshold for module membership.
#' @param degree_thresh Threshold for degree.
#'
#' TFs with high module membership (r>MM_thresh) and belonging to the top 10% (e.g. degree_thresh = 0.9) of genes with high degree are selected as hub TFs.
#'
#' @return A list with hub TFs per module.
#' @export
#'
#' @examples
#'
#' hub_tfs = identify_hub_TFs(t(tfs), network, MM_thresh = 0.8, degree_thresh = 0.9)
#'
identify_hub_TFs <- function(datExpr, TF.network, MM_thresh = 0.8, degree_thresh = 0.9) {
  moduleEigengenes = TF.network[[1]]
  moduleColors = TF.network[[2]]
  # Calculate Module Membership (MM)
  moduleMemberships <- sapply(unique(moduleColors), function(module) {
    genesInModule <- which(moduleColors == module)
    eigengene <- moduleEigengenes[, paste0("ME", module)]
    WGCNA::cor(t(datExpr[genesInModule, ]), eigengene)
  })

  # Calculate adjacency matrices and degrees
  adjacencyList <- lapply(unique(moduleColors), function(module) {
    genesInModule <- which(moduleColors == module)
    moduleData <- datExpr[genesInModule, ]
    adjacency <- WGCNA::cor(t(moduleData))
    adjacency[lower.tri(adjacency)] <- 0
    adjacency
  })

  names(adjacencyList) <- unique(moduleColors)

  moduleDegrees <- lapply(unique(moduleColors), function(module) {
    adjacency <- adjacencyList[[module]]
    degree <- rowSums(adjacency)
    names(degree) <- rownames(datExpr)[which(moduleColors == module)]
    degree
  })
  names(moduleDegrees) <- unique(moduleColors)

  # Identify hub genes
  hubGenesList <- list()
  hubGenesData <- data.frame()  # Initialize an empty data frame

  allDegrees <- numeric()
  allMemberships <- numeric()

  for (module in unique(moduleColors)) {
    genesInModule <- which(moduleColors == module)
    degrees <- moduleDegrees[[module]]
    memberships <- moduleMemberships[[module]]

    # Update the global lists of degrees and memberships
    allDegrees <- c(allDegrees, degrees)
    allMemberships <- c(allMemberships, memberships)

    # Calculate the cutoff for the top 10% by degree
    degreeCutoff <- stats::quantile(degrees, degree_thresh)

    # Plot distribution of Degrees for the current module
    # degree_plot_data <- data.frame(Degree = degrees)
    #
    # degree_plot <- ggplot(degree_plot_data, aes(x = Degree)) +
    #   geom_histogram(binwidth = 5, fill = module, color = "black", alpha = 0.6) +
    #   geom_vline(xintercept = degreeCutoff, linetype = "dashed", color = "red") +
    #   labs(
    #     title = paste("Distribution of Gene Degrees in Module", module),
    #     x = "Degree",
    #     y = "Frequency"
    #   ) +
    #   theme_minimal()
    #
    # print(degree_plot)

    # Plot distribution of Module Membership for the current module
    # membership_plot_data <- data.frame(ModuleMembership = memberships)
    #
    # membership_plot <- ggplot(membership_plot_data, aes(x = ModuleMembership)) +
    #   geom_histogram(binwidth = 0.05, fill = module, color = "black", alpha = 0.6) +
    #   geom_vline(xintercept = MM_thresh, linetype = "dashed", color = "red") +
    #   labs(
    #     title = paste("Distribution of Module Membership in Module", module),
    #     x = "Module Membership",
    #     y = "Frequency"
    #   ) +
    #   theme_minimal()
    #
    # print(membership_plot)

    # Identify hub genes (top 10% by degree)
    hubGenes <- names(degrees[degrees >= degreeCutoff])
    if (length(hubGenes) == 0) next

    # Identify hub genes (>0.8 MM)
    highMMGenes <- which(memberships > MM_thresh)
    if (length(highMMGenes) == 0) next

    # Filter hub genes (degree) to only include those with high MM
    finalHubGenes <- intersect(hubGenes, rownames(memberships)[highMMGenes])
    if (length(finalHubGenes) == 0) next

    hubGenesList[[module]] = finalHubGenes

    # Create a dataframe with detailed information for hub genes only
    moduleData <- data.frame(
      Module = module,
      ModuleMembership = unname(memberships[finalHubGenes,]),
      Degree = degrees[finalHubGenes],
      stringsAsFactors = FALSE)

    # Append to the combined dataframe
    hubGenesData <- rbind(hubGenesData, moduleData)
  }

  return(list(hubGenes = hubGenesList, detailedData = hubGenesData))
}

#' Identify cell groups
#'
#' Identifies cell dendrograms corresponding to each TF module group.
#'
#' @param features Significant deconvolution features associated with TF modules (pvalue < 0.05)
#' @param tfs.modules.groups TF modules groups obtained from compute.TF.network.classification()
#' @param cor_type Type of correlation. Either pearson "p" or spearman "s"
#' @param clustering.method Clustering method. Default is "ward.D2", see hclust() for all methods available.
#' @param distance.method Distance method. Default is "euclidean", see dist() for all methods available.
#' @param width Numeric value to define the width of plots saved in Results/ directory.
#' @param height Numeric value to define the height of plots saved in Results/ directory.
#' @param return Whether to save or not the plots in the Results/ folder.
#'
#' @return A list with the cell dendrograms per TF module group.
#' @export
#'
#' @examples
#'
#' cell_dendrograms = identify.cell.groups(corr_modules, tfs.modules.clusters, height = 20, return = T)
#'
identify.cell.groups = function(features, tfs.modules.groups, cor_type = "p", clustering.method = "ward.D2", distance.method = "euclidean", width = 12, height = 18, return = T){

  moduleTraitCor = features[[1]]

  names(features[[2]]) = paste0("ME", names(features[[2]])) #To match names of columns from corr matrix

  #Gather significant features across TF modules clusters
  features_vec = list()
  for (i in 1:length(tfs.modules.groups)) {
    features_vec[[i]] = unique(unlist(unname(features[[2]][tfs.modules.groups[[i]]])))
    names(features_vec)[i] = names(tfs.modules.groups)[i]
  }

  lis.dendrogram = list()

  for (i in 1:length(features_vec)){
    TFmoduleTraitcor = moduleTraitCor[,colnames(moduleTraitCor)%in%features_vec[[i]]]
    TFmoduleTraitcor = TFmoduleTraitcor[rownames(TFmoduleTraitcor)%in%tfs.modules.groups[[i]], , drop=F]
    data_scaled = scale(t(TFmoduleTraitcor)) #Scale per TF module
    ###Dendogram by Module
    d <- stats::dist(data_scaled, method = distance.method)
    d = d/sqrt(ncol(data_scaled)) #Adjust/Scale distance matrix for number of features to make dendrograms comparable
    dendrogram <- stats::hclust(d, method = clustering.method)
    if(return){
      pdf(paste0("Results/Dendogram_cell_types_", names(features_vec)[[i]]), width = width, height = height)
      par(mar = c(5, 2, 4, 35)) #bottom, left, top, right
      plot(as.dendrogram(dendrogram), horiz= T)
      dev.off()
    }
    lis.dendrogram[[i]] = dendrogram
  }

  names(lis.dendrogram) = names(features_vec)

  # for (i in 1:length(features_vec)){
  #   moduleTraitCor[i,!(colnames(moduleTraitCor)%in%features_vec[[i]])] = 0 ###Set to 0 not significant correlations
  # }

  #Add dendrogram "all" considering all TFs modules

  # data_scaled = scale(t(moduleTraitCor)) #Scale per TF module
  # d <- dist(data_scaled, method = "euclidean")
  # d = d/sqrt(ncol(data_scaled))
  # dendrogram_all <- hclust(d, method = clustering.method)
  # if(return){
  #   pdf("Results/Dendogram_cell_types_all", width = width, height = height)
  #   par(mar = c(5, 2, 4, 35)) #bottom, left, top, right
  #   plot(as.dendrogram(dendrogram_all), horiz= T)
  #   dev.off()
  # }
  #
  # lis.dendrogram[[length(lis.dendrogram)+1]] = dendrogram_all
  # names(lis.dendrogram)[[length(lis.dendrogram)]] = "all"

  return(lis.dendrogram)

}
