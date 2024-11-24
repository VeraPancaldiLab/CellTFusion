
##Functions from https://github.com/VeraPancaldiLab/CellTFusion

##Load libraries

libraries_set <- function(){
  suppressMessages(library("BiocManager"))
  suppressMessages(library("devtools"))
  suppressMessages(library("pak"))
  suppressMessages(library("remotes"))
  suppressMessages(library("decoupleR"))
  suppressMessages(library("OmnipathR"))
  suppressMessages(library("tidyr"))
  suppressMessages(library("dplyr"))
  suppressMessages(library("matrixStats"))
  suppressMessages(library("org.Hs.eg.db"))
  suppressMessages(library("ReactomePA"))
  suppressMessages(library("WGCNA"))
  suppressMessages(library("reshape2"))
  suppressMessages(library("purrr"))
  suppressMessages(library("tidygraph"))
  suppressMessages(library("stringr"))
  suppressMessages(library("tibble"))
  suppressMessages(library("gplots"))
  suppressMessages(library("ggplot2"))
  suppressMessages(library("AnnotationDbi"))
  suppressMessages(library("DESeq2"))
  suppressMessages(library("RColorBrewer"))
  suppressMessages(library("pheatmap"))
  suppressMessages(library("ggfortify"))
  suppressMessages(library("viper"))
  suppressMessages(library("msigdbr"))
  suppressMessages(library("GSVA"))
  suppressMessages(library("clusterProfiler"))
  suppressMessages(library("Hmisc"))
  suppressMessages(library("DOSE"))
  suppressMessages(library("fgsea"))
  suppressMessages(library("ggpubr"))
  suppressMessages(library("ComplexHeatmap"))
  suppressMessages(library("ggstatsplot"))
  suppressMessages(library("dendextend"))
  suppressMessages(library("PCAtools"))
  suppressMessages(library("stats"))
  suppressMessages(library("Boruta"))
  suppressMessages(library("caret"))
  suppressMessages(library("pROC"))
  suppressMessages(library("MLeval"))
  suppressMessages(library("survival"))
  suppressMessages(library("survminer"))
  suppressMessages(library("rms"))
  suppressMessages(library("BulkSignalR"))
  suppressMessages(library("igraph"))
  suppressMessages(library("enrichplot"))
  suppressMessages(library("uuid"))
  suppressMessages(library("parallel"))
  suppressMessages(library("factoextra"))
  suppressMessages(library("doParallel"))
  suppressMessages(library("foreach"))
  suppressMessages(library("omnideconv"))
}

libraries_set()
source("../src/cell_deconvolution.R") #Functions for deconvolution
source("../src/machine_learning.R") #Functions for machine learning

dir.create(file.path(getwd(), "Results"))

#' TPM normalization
#'
#' \code{TPM_normalization} Compute TPM normalization
#'
#' @param data matrix.
#'
#' @return A matrix of TPM values (genes X samples)
#'


vsd_normalization = function(raw.counts, coldata){
  dds <- DESeqDataSetFromMatrix(countData = round(raw.counts),
                                colData = coldata,
                                design = ~1)
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds)
  normalized.counts <- counts(dds, normalized=TRUE)

  return(normalized.counts)
}

remove.outliers = function(counts){
  sampleDists <- dist(t(counts))
  sampleDistMatrix <- as.matrix(sampleDists)
  colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pdf("Results/SampleDistanceMatrix.pdf", width = 16, height = 14)
  print(pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists,
           clustering_distance_cols=sampleDists,
           col=colors))
  dev.off()
  sampleTree = hclust(dist(t(counts)), method = "average");
  pdf("Results/SamplesClustering.pdf", width = 16)
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
  dev.off()
  response = readline(prompt = "Do you have outliers? (Y|N): ")
  if(response == 'Y'){
    cutheight = as.numeric(readline(prompt = "Cut height parameter: "))
    pdf("Results/SamplesClusteringOutliers.pdf", width = 16)
    plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
    abline(h = cutheight, col = "red")
    dev.off()
    clust = cutreeStatic(sampleTree, cutHeight = cutheight, minSize = 10) # Determine cluster under the line
    table(clust) # clust 1 contains the samples we want to keep.
    keepSamples = (clust==1)
    counts = counts[,keepSamples]
  }

  return(counts)
}

convert_ensembl_symbol = function(counts){
  require(org.Hs.eg.db)
  counts = data.frame(counts)
  # Map gene SYMBOLS
  entrz <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(counts), columns = "SYMBOL", keytype = "ENSEMBL") %>%
    distinct(SYMBOL, .keep_all=T) %>% #Remove duplicated symbols
    filter(ENSEMBL %in% rownames(counts)) %>% #Take only ENSEMBL ids in counts
    group_by(ENSEMBL) %>%
    slice(1) %>% #If two same ENSEMBL ids corresponds to different symbols, take only one
    ungroup()

  # Change rownames to SYMBOL
  counts = counts %>%
    rownames_to_column("Genes") %>%
    filter(Genes %in% entrz$ENSEMBL) %>%
    mutate(Symbol = entrz$SYMBOL) %>%
    na.omit() %>%
    remove_rownames() %>%
    column_to_rownames("Symbol") %>%
    dplyr::select(-Genes)

  return(counts)
}

filter_coding_genes = function(counts, annotation = "ENSEMBL"){

  require(biomaRt)

  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  gene_ids <- rownames(counts)

  if(annotation == "ENSEMBL"){

    # Get gene type annotations
    annotations <- getBM(attributes = c("ensembl_gene_id", "gene_biotype"),
                         filters = "ensembl_gene_id",
                         values = gene_ids,
                         mart = mart)

    # Filter annotations to keep only coding genes
    coding_genes <- annotations %>%
      filter(gene_biotype == "protein_coding") %>%
      pull(ensembl_gene_id)

    counts = counts[rownames(counts) %in% coding_genes,]

  }else if(annotation == "SYMBOL"){

    # Get gene type annotations
    annotations <- getBM(attributes = c("hgnc_symbol", "gene_biotype"),
                         filters = "hgnc_symbol",
                         values = gene_ids,
                         mart = mart)

    # Filter annotations to keep only coding genes
    coding_genes <- annotations %>%
      filter(gene_biotype == "protein_coding") %>%
      pull(hgnc_symbol)

    counts = counts[rownames(counts) %in% coding_genes,]
  }else{
    stop("Annotation provided not supported")
  }

  return(counts)

}

#Compute PCA analysis
compute_pca_analysis <- function(data, coldata, trait, trait2 = NULL, ncomp = 5, cortype = "pearson"){
  #data: genes as rows and samples as columns
  pca_res <- prcomp(t(data))
  coldata[,colnames(coldata)%in%trait] = factor(coldata[,colnames(coldata)%in%trait])
  print(autoplot(pca_res, data = coldata, colour = trait, shape = trait2, size  = 3))

  p <- PCAtools::pca(data, metadata = coldata, removeVar = 0.1)

  peigencor  <- eigencorplot(p,
                             components = PCAtools::getComponents(p, 1:ncomp),
                             metavars = colnames(coldata),
                             col = c('white', 'cornsilk1', 'gold', 'forestgreen', 'darkgreen'),
                             cexCorval = 1.2,
                             fontCorval = 2,
                             posLab = 'all',
                             rotLabX = 45,
                             scale = TRUE,
                             main = paste("Principal component", cortype, "r^2 clinical correlates"),
                             plotRsquared = TRUE,
                             corFUN = cortype,
                             corUSE = 'pairwise.complete.obs',
                             corMultipleTestCorrection = 'BH',
                             signifSymbols = c('****', '***', '**', '*', ''),
                             signifCutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1))

  print(peigencor)
}

#' Inference of TFs activity based on gene expression
#'
#' \code{compute.TFs.activity} computes TFs activty based on a gene expression matrix using VIPER algorithm and a collection of GRN (collectri) from the OmnipathR package
#'
#' @param RNA.tpm Gene expression matrix normalized by TPM (genes X samples).
#' @return A matrix of protein activity (samples X tfs).
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
    net_regulons <- aracne2regulon(network_file, as.matrix(RNA.counts), format = "3col")
  }

  if(tfs.pruned == TRUE){
    net_regulons = pruneRegulon(net_regulons, cutoff = max_size_targets)
  }

  sample_acts <- viper(as.matrix(RNA.counts), net_regulons, minsize = min_targets_size, verbose=F, method = "scale")
  message("TFs scores computed")

  return(data.frame(t(sample_acts)))

}

#' Pathways activity computation using a curated collection of pathways and their target genes (PROGENy)
#'
#' \code{compute.pathway.activity} computes pathways activty based on a protein activity matrix using a multivariate linear model (mlm)
#'
#' @param RNA.tpm Gene expression matrix normalized by TPM (genes X samples).
#' @return A matrix of pathway activity (samples X pathways).
#'
compute.pathway.activity <- function(RNA.tpm, gene_sets = NULL, paths = NULL){

  RNA.tpm = as.matrix(RNA.tpm)
  #Get universe
  if(is.null(paths)){
    paths <- get_progeny(organism = 'human', top = 500)
  }

  # Run mlm
  progeny <- run_mlm(mat=RNA.tpm, net=paths, .source='source', .target='target', .mor='weight', minsize = 5)

  #Remove variable
  rm(paths)
  gc()

  # Transform to wide matrix
  sample_acts_progeny <- progeny %>%
    pivot_wider(id_cols = 'condition', names_from = 'source',
                values_from = 'score') %>%
    column_to_rownames('condition') %>%
    as.matrix()

  if(is.null(gene_sets)==F){

    cat("Computing GSVA analysis using provided gene sets.....................................................\n")

    gsva_results <- gsva(
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

#' Enrichment of modules using Reactome database
#'
#' \code{compute.modules.enrichment} computes TFs enrichment using directed targets for each module
#'
#' @param RNA.tpm Gene expression matrix normalized by TPM (genes X samples).
#' @param TFs.matrix TFs activity matrix (TFs X samples).
#' @param ME.matrix Matrix of modules eigenvectors (samples X modules).
#' @param module_colors Named list with TFs and their respective module color assignation.
#'
#' @return Enrichment plots for each module.
#'
compute.modules.enrichment <- function(RNA.tpm, hub_tfs){
  net = get_collectri(organism = 'human', split_complexes = F) #Get universe
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

  ItemsList <- venn(pathways, show.plot = FALSE)
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
        print(dotplot(res[[i]],  title=paste0("Enrichment by Reactome\nModule ", color)))
        dev.off()
      }
    }
  }

}

module_enrich = function(tpm.counts, module_color, hub_genes, tfs_universe){
  # genes = colnames(TFs.matrix)
  # inModule = is.finite(match(module_colors,module))
  # modGenes = genes[inModule]
  targets = tfs_universe[tfs_universe$source %in% hub_genes[[1]][[module_color]],] #Extract targets from TFs
  targets = unique(targets$target) #Keep only unique targets from TFs

  targets_genes = tpm.counts[rownames(tpm.counts)%in%targets,] #Extract gene expression from targets
  targets_genes = targets_genes[order(rowVars(targets_genes), decreasing = T),][1:round(0.2*nrow(targets_genes)),] #Keep only highly variable targets (20%)

  entrz <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(targets_genes), columns = "ENTREZID", keytype = "SYMBOL") #Change to EntrezID
  universe = AnnotationDbi::select(org.Hs.eg.db, keys = rownames(tpm.counts), columns = "ENTREZID", keytype = "SYMBOL") #Change to EntrezID

  reac <- enrichPathway(gene    = entrz$ENTREZID,
                        organism     = 'human',
                        universe = universe$ENTREZID,
                        pvalueCutoff = 0.05)

  reac@result = reac@result[reac@result$p.adjust<0.05,]

  if(nrow(reac@result)!=0){
    return(reac)
  }

}

identify_hub_TFs <- function(datExpr, TF.network, MM_thresh = 0.8, degree_thresh = 0.9) {
  moduleEigengenes = TF.network[[1]]
  moduleColors = TF.network[[2]]
  # Calculate Module Membership (MM)
  moduleMemberships <- sapply(unique(moduleColors), function(module) {
    genesInModule <- which(moduleColors == module)
    eigengene <- moduleEigengenes[, paste0("ME", module)]
    cor(t(datExpr[genesInModule, ]), eigengene)
  })

  # Calculate adjacency matrices and degrees
  adjacencyList <- lapply(unique(moduleColors), function(module) {
    genesInModule <- which(moduleColors == module)
    moduleData <- datExpr[genesInModule, ]
    adjacency <- cor(t(moduleData))
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
    degreeCutoff <- quantile(degrees, degree_thresh)

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


#' Construction of Weighted TFs co-activity network
#'
#' \code{compute.WTCNA} Construct a weighted signed or unsigned network using TFs activity to cluster proteins regulators into modules that share similar activity patterns
#'
#' @param TFs.matrix TFs activity matrix (samples x TFs).
#' @param network.type Allowed values are (unique abbreviations of) "unsigned", "signed", "signed hybrid", "distance". See WGCNA documentation.
#' @param clustering.method The agglomeration method to be used. See hclust documentation.
#' @param minMod Minimun number of TFs allowed for each module.

#' @return A list including matrix of modules eigenvectors (samples X modules) and named list with TFs and their respective module color assignation.
#'
compute.WTCNA <- function(TFs.matrix, network.type = "unsigned", clustering.method = "ward.D2", minMod = 15, corr_mod = 0.9, cor_type = "p", return = T){

  cat("Creating weighted TF-coactivity network......................................................................\n\n")
  #####Choose parameter for scale-free network topology
  powers = c(c(1:10), seq(from = 12, to=20, by=1))
  sft = pickSoftThreshold(TFs.matrix, powerVector = powers, verbose = 0, networkType = network.type)

  if(return){
    pdf("Results/Soft Threshold")
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
         main = paste("Scale independence"))
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=0.9,col="red");
    abline(h=0.90,col="red")
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
  adjacency = adjacency(TFs.matrix, power =softPower, type=network.type, corFnc = "cor", corOptions = list(use = cor_type))
  TOM = TOMsimilarity(adjacency, TOMType = network.type)
  dissTOM = 1-TOM

  #####Unsupervised hierarchical clustering using dissimilarity matrix
  geneTree = hclust(dist(dissTOM), method = clustering.method)
  dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                              deepSplit = 2, pamRespectsDendro = FALSE,
                              minClusterSize = minMod);
  dynamicColors = labels2colors(dynamicMods)

  #####Remove variables and clean garbage
  rm(adjacency, TOM, dissTOM)
  gc()

  if(return){
    pdf("Results/Gene dendrogram and module colors")
    plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05,
                        main = "Gene dendrogram and module colors")
    dev.off()
  }

  #####Calculate eigenvectors from modules
  cat("Calculating eigenvectors from modules.................................................\n\n")
  MEList = moduleEigengenes(TFs.matrix, colors = dynamicColors, scale = F) #Data already scale
  MEs = MEList$eigengenes
  MEs = orderMEs(MEs)

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
    plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                        dendroLabels = FALSE, hang = 0.03,
                        addGuide = TRUE, guideHang = 0.05,
                        main = "Gene dendrogram and module colors")
    dev.off()
  }

  TFspropVar = propVarExplained(TFs.matrix, dynamicColors, MEs, corFnc = "cor", corOptions = "use = 'p'")

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

#####Merge modules that are significantly correlated
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
          mutate(new_column = rowMeans(dplyr::select(., module1, module2))) %>%
          dplyr::rename(module1 = new_column) %>%
          dplyr::select(., -module1, -module2)
      }
    }
  }

  return(list(data, colors))
}

#' Modules relationship across features
#'
#' \code{compute.modules.relationship} computes pearson correlation between two indicate matrices
#'
#' @param tfs_network TFs modules to correlate
#' @param matB Dataframe B (samples X features) to correlate.
#' @param file_name name file for plots.
#' @param width Width pdf file for plot (default = 8).
#' @param height Height pdf file for plot (default = 8).
#'
#' @return Module relationship plot between matrices.
#'
compute.modules.relationship <- function(tfs_network, matB, file_name, width = 8, height = 8, pval=0.05, padj = F, cor_type = "p", return = F, vertical=F, plot = T){

  tfs_network = data.frame(tfs_network)
  matB = data.frame(matB)

  ##check if names from both features are the same
  if(all(rownames(tfs_network)==rownames(matB)) == F){
    stop("No equal names, verify the input objects")
  }


  moduleTraitCor = cor(tfs_network, matB, method = cor_type)
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(tfs_network))

  rev = which(colSums(moduleTraitPvalue > pval)==nrow(moduleTraitPvalue)) #check if there are features no significant with any module

  if(length(rev)>0){
    moduleTraitCor = moduleTraitCor[,-rev]
    moduleTraitPvalue = moduleTraitPvalue[,-rev]
  }

  if(padj == T){
    for (i in 1:ncol(moduleTraitPvalue)) {
      moduleTraitPvalue[,i] = p.adjust(moduleTraitPvalue[,i], method = 'bonferroni')
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
        d <- dist(t(moduleTraitCor), method = "manhattan")
        hc1 <- hclust(d, method = "ward.D2")
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
          labeledHeatmap(Matrix = moduleTraitCor[vec,],
                         xLabels = colnames(moduleTraitCor),
                         yLabels = rownames(moduleTraitCor[vec,]),
                         xLabelsPosition = "top",
                         colors = blueWhiteRed(50),
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
            labeledHeatmap(Matrix = moduleTraitCor,
                           xLabels = colnames(moduleTraitCor),
                           yLabels = rownames(moduleTraitCor),
                           xLabelsPosition = "top",
                           colors = blueWhiteRed(50),
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
        d <- dist(t(moduleTraitCor), method = "manhattan")
        hc1 <- hclust(d, method = "ward.D2")
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
          labeledHeatmap(Matrix = moduleTraitCor[,vec],
                         xLabels = names(moduleTraitCor[,vec]),
                         yLabels = rownames(moduleTraitCor),
                         ySymbols = rownames(moduleTraitCor),
                         colorLabels = FALSE,
                         colors = blueWhiteRed(50),
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
            labeledHeatmap(Matrix = moduleTraitCor,
                           xLabels = names(moduleTraitCor),
                           yLabels = rownames(moduleTraitCor),
                           ySymbols = rownames(moduleTraitCor),
                           colorLabels = FALSE,
                           colors = blueWhiteRed(50),
                           textMatrix = textMatrix,
                           setStdMargins = FALSE,
                           cex.text = 0.5,
                           zlim = c(-1,1),
                           main = paste("Module-trait relationships"))
            dev.off()
          }
        }}}
}

classify.tfs = function(tfs.modules, col.data, color){
  col.data = col.data %>%
    mutate(Module = tfs.modules[,paste0("ME", color)],
           Module_level = ifelse(Module > summary(Module)[5], 'High',
                                 ifelse(Module < summary(Module)[2], "Low", "na")))
  return(col.data)
}

diff_analysis = function(RNA.counts.normalized, feature, file.name){

  ####TFs differential active
  net = get_collectri(organism = 'human', split_complexes = F) #Get universe

  collectri2viper_regulons <- function(df) {
    regulon_list <- split(df, df$source)
    regulons <- lapply(regulon_list, function(regulon) {
      tfmode <- stats::setNames(regulon$mor, regulon$target)
      list(tfmode = tfmode, likelihood = rep(1, length(tfmode)))
    })

    return(regulons)
  }

  net_regulons = collectri2viper_regulons(net)
  #net_regulons <- pruneRegulon(net_regulons, 50, adaptive = FALSE, eliminate = TRUE) #Pruned regulons to a maximum of 50 targets to avoid statistic bias

  # Generating test and ref data
  test_i <- which(feature == 1)
  ref_i <- which(feature == 2)

  mat_test <- as.matrix(RNA.counts.normalized[,test_i])
  mat_ref <- as.matrix(RNA.counts.normalized[,ref_i])

  # Generating NULL model (test, reference)
  dnull <- ttestNull(mat_test, mat_ref, per=1000, verbose = F)

  # Generating signature
  signature <- rowTtest(mat_test, mat_ref)
  signature <- (qnorm(signature$p.value/2, lower.tail = FALSE) * sign(signature$statistic))
  signature <- na.omit(signature)
  signature <- signature[,1]

  # Running msVIPER
  mra <- msviper(signature, net_regulons, dnull, verbose = F, minsize = 5)

  # Plot DiffActive TFs
  pdf(paste0("Results/Differential_TFs_", file.name))
  print(plot(mra, mrs=15, cex=1, include = c("expression","activity")))
  dev.off()

  tfs_da = mra$es$p.value
  tfs_da = tfs_da[tfs_da < 0.05]
  tfs_names = names(tfs_da)

  nes = mra$es$nes
  nes = nes[tfs_names] #only significant
  up = names(nes[nes > 0])
  down = names(nes[nes<0])

  #Map TFs targets
  targets_up = net$target[which(net$source%in%up)]
  targets_down = net$target[which(net$source%in%down)]
  common = intersect(targets_up, targets_down)
  if(length(common)!=0){
    targets_up = targets_up[-which(targets_up %in% common)]
    targets_down = targets_down[-which(targets_down %in% common)]
  }

  #Enrichment upregulated TFs
  entrz <- AnnotationDbi::select(org.Hs.eg.db, keys = targets_up, columns = "ENTREZID", keytype = "SYMBOL") #Change to EntrezID
  universe <- AnnotationDbi::select(org.Hs.eg.db, keys = rownames(RNA.counts.normalized), columns = "ENTREZID", keytype = "SYMBOL") #Change to EntrezID

  reac <- enrichPathway(gene    = entrz$ENTREZID,
                        organism     = 'human',
                        universe = universe$ENTREZID,
                        pvalueCutoff = 0.05)

  kegg <- enrichKEGG(gene    = entrz$ENTREZID,
                     organism     = 'human',
                     universe = universe$ENTREZID,
                     pvalueCutoff = 0.05)

  if(nrow(data.frame(kegg))!=0){
    pdf(paste0("Results/ORA_KEGG_UP_TFs_", file.name))
    print(dotplot(kegg))
    dev.off()
  }

  if(nrow(data.frame(reac))!=0){
    pdf(paste0("Results/ORA_Reactome_UP_TFs_", file.name))
    print(dotplot(reac))
    dev.off()
  }

  return(up)

}


anova_test = function(data, trait_name, y_name, pval = 0.05, file_name){
  print("Performing one way ANOVA test...................................")
  data$y = data[,y_name]
  data$trait = as.factor(data[,trait_name])
  model  <- lm(y ~ trait, data = data)
  res.aov <- data %>% rstatix::anova_test(y ~ trait)
  if(round(res.aov$p, 5) <= pval){
    print("Significant test found!")
    pdf(paste0("Results/ANOVA_", file_name))
    print(ggplot(data, aes(x=trait, y=y, fill=trait)) +
            geom_violin(width=0.6) +
            geom_boxplot(width=0.07, color="black", alpha=0.2) +
            scale_fill_brewer() +
            geom_smooth(aes(x=trait, y=y), method = "loess") +
            ylab(paste0("Values for ", y_name)) +
            xlab(paste0("Clinical trait: ", trait_name)) +
            labs(title="One way ANOVA test",
                 subtitle=rstatix::get_test_label(res.aov, detailed = TRUE)) +
            theme(axis.text.x = element_text(angle = 0),
                  axis.title.y = element_text(size = 8, angle = 90)))
    dev.off()
  }
}


minMax <- function(x) {
  #columns: features
  x = data.matrix(x)
  for(i in 1:ncol(x)){
    x[,i] = (x[,i] - min(x[,i], na.rm = T)) / (max(x[,i], na.rm = T) - min(x[,i], na.rm = T))
  }

  return(x)

}

check_normal_distribution = function(data_df){
  data_df = data.frame(data_df)
  # Visualize the data
  plot_histograms <- function(column, colname) {
    if(is.numeric(column)) {
      ggplot(data_df, aes_string(x = colname)) +
        geom_histogram(aes(y = ..density..), bins = 30, fill = "blue", alpha = 0.5) +
        geom_density(color = "red") +
        ggtitle(paste("Histogram and Density Plot for", colname))
    }
  }

  # Plot histograms
  for (colname in names(data_df)) {
    print(plot_histograms(data_df[[colname]], colname))
  }
}

is_scaled <- function(x) {
  mean_x <- mean(x)
  sd_x <- sd(x)
  return(abs(mean_x) < 1e-10 && abs(sd_x - 1) < 1e-10)
}

get_best_mfrow <- function(n_plots) {
  rows <- floor(sqrt(n_plots))
  cols <- ceiling(n_plots / rows)
  return(c(rows, cols))
}

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
    d <- dist(data_scaled, method = distance.method)
    d = d/sqrt(ncol(data_scaled)) #Adjust/Scale distance matrix for number of features to make dendrograms comparable
    dendrogram <- hclust(d, method = clustering.method)
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

classify.deconvolution = function(coldata, deconvolution, group){
  deconv = deconvolution[,colnames(deconvolution)%in%group]

  #Patients high in group 1 of cells
  vec = c()
  if(is.null(ncol(deconv))==T){
    idx = which(deconv > median(deconv))
    vec = c(vec,idx)
  }else{
    for (i in 1:ncol(deconv)) {
      idx = which(deconv[,i] > median(deconv[,i]))
      vec = c(vec, idx)
    }
  }

  #High in all deconv features from group
  pos = which(table(vec) == length(group))

  coldata = coldata %>%
    mutate(Cells_level = "Low")

  coldata$Cells_level[pos] = "High"
  coldata$Cells_level = factor(coldata$Cells_level)

  return(coldata)

}

compute.fisher.test = function(coldata, trait){

  contingency = table(coldata[,"Cells_level"], coldata[,trait])
  test <- fisher.test(contingency)

  df = data.frame("Cells_level" = coldata[,"Cells_level"],
                  "Trait" = coldata[,trait])

  pdf(paste0("Results/Fisher test for ", trait), width = 12, height = 9)
  ggbarstats(
    df, Cells_level, Trait,
    results.subtitle = FALSE,
    title = paste0("Analysis for ", trait),
    subtitle = paste0(
      "Fisher's exact test", ", p-value = ",
      ifelse(test$p.value < 0.001, "< 0.001", round(test$p.value, 3))
    )
  )
}

cell.groups.anova.test = function(cell.groups, coldata, trait, pval = 0.05){

  sig = c()
  coldata[,trait] = as.factor(coldata[,trait])

  for (j in 1:ncol(cell.groups[[1]])) {
    data = data.frame("Value" = cell.groups[[1]][,j], "Trait" = coldata[,trait])
    model  <- lm(Value ~ Trait, data = data)
    res.aov <- data %>% rstatix::anova_test(Value ~ Trait)

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

cell.groups.fisher.test = function(cell.groups, coldata, trait, pval = 0.05){

  sig = c()
  coldata[,trait] = as.factor(coldata[,trait])

  for (j in 1:ncol(cell.groups[[1]])) {
    coldata = coldata %>%
      mutate(level = cell.groups[[1]][,j],
             Cells_level = ifelse(level > summary(level)[3], 'High', 'Low'))

    contingency = table(coldata[,"Cells_level"], coldata[,trait])
    test = fisher.test(contingency)

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


compute.TF.network.classification = function(tf.network, pathways.features, return = T){

  tf.network = data.frame(tf.network[[1]])
  pathways.features = data.frame(pathways.features)

  moduleTraitCor = cor(tf.network, pathways.features, method = "p")
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
  silhouette = fviz_nbclust(moduleTraitCor, hcut, method = "silhouette", k.max = nrow(moduleTraitCor)-1)
  k_cluster = as.numeric(silhouette$data$clusters[which.max(silhouette$data$y)])

  if(return){
    pdf(paste0("Results/TFs_modules_Silhouette_scores"))
    print(silhouette)
    dev.off()
  }

  hc_modules = hclust(dist(moduleTraitCor), method = "ward.D2")
  dend_pathways = as.dendrogram(hc_modules)

  if(return){
    pdf(paste0("Results/TFs_modules_clusters"))
    plot(dend_pathways, cex = 0.6)
    rect.hclust(hc_modules, k = k_cluster, border = 2:5)
    dev.off()
  }

  ### Extract clusters
  sub_grp <- cutree(hc_modules, k = k_cluster)

  ### Plot PCA and biplot
  p = fviz_cluster(list(data = moduleTraitCor, cluster = sub_grp))

  if(return){
    pdf(paste0("Results/PCA_TFs_modules_clusters"))
    print(p)
    dev.off()
  }

  res.pca <- prcomp(moduleTraitCor,  scale = F)
  p = fviz_pca_biplot(res.pca, label="all", select.var = list(contrib = 6), addEllipses=TRUE, ellipse.level=0.75)

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
    rownames_to_column("Features") %>%
    arrange(desc(PC1))

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

calculate_dendrogram_cuts = function(cell.group.dendrogram, n_cuts = NULL){

  cuts = list()
  for (i in 1:length(cell.group.dendrogram)) {
    dend_heights <- dendextend::heights_per_k.dendrogram(as.dendrogram(cell.group.dendrogram[[i]])) #Calculate dendrogram heights

    sorted_heights <- sort(dend_heights) # Sort heights
    height_diffs <- diff(sorted_heights) # Calculate differences between consecutive heights

    buffer <- max(median(height_diffs[height_diffs > 0]), 1) # Buffer: take the median of the non-zero differences

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


plot.modules.categorical = function(tfs.modules, coldata){
  data = cbind(tfs.modules, coldata)
  for(i in 1:ncol(tfs.modules)){
    for (j in (ncol(tfs.modules)+1):ncol(data)) {
      module <- names(data[i])
      trait <- names(data[j])
      avz <- broom::tidy(aov(data[,i] ~ data[,j], data = data))
      if(avz$p.value[1] < 0.05) {
        pdf(paste0("Results/ANOVA_", module, "-", trait))
        print(ggplot(data, aes(x=data[,j], y=data[,i], fill=data[,j])) +
                geom_violin(width=0.6) +
                geom_boxplot(width=0.07, color="black", alpha=0.2) +
                scale_fill_brewer() +
                geom_smooth(aes(x=data[,j], y=data[,i]), method = "loess") +
                ylab(paste0("Values for ", module)) +
                xlab(paste0("Clinical trait: ", trait)) +
                labs(title="One way ANOVA test",
                     subtitle=paste0("pvalue: ", avz$p.value[1])) +
                theme(axis.text.x = element_text(angle = 0),
                      axis.title.y = element_text(size = 8, angle = 90))+
                scale_fill_discrete(name = trait))
        dev.off()
      }
    }
  }
}

compute.metada.association = function(tfs.modules, coldata, pval = 0.05, width = 20, height = 8){
  ###Association with categorical variables
  coldata_categorical = coldata %>%
    dplyr::select(where(is.character)|where(is.factor))

  if(ncol(coldata_categorical)!=0){
    data = cbind(tfs.modules, coldata_categorical)
    pvals = data.frame()
    fvals = data.frame()
    for(i in 1:ncol(tfs.modules)){
      contador = 1
      for (j in (ncol(tfs.modules)+1):ncol(data)) {
        module <- names(data[i])
        trait <- names(data[j])
        avz <- broom::tidy(aov(data[,i] ~ data[,j], data = data))
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
    dplyr::select(where(is.numeric))

  if(ncol(coldata_quantitative)!=0){
    moduleTraitCor = cor(tfs.modules, coldata_quantitative, method = "p");
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(tfs.modules))

    textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 2), ")", sep = "")
    dim(textMatrix) = dim(moduleTraitCor)

    if(ncol(coldata_categorical)!=0){
      textMatrix = cbind(textMatrix, textMatrix2)
      moduleTraitPvalue = cbind(moduleTraitPvalue, pvals)
      simulated_corr = matrix(runif(n=nrow(pvals)*ncol(pvals), min=-0.1, max=0.1), nrow = nrow(pvals), ncol = ncol(pvals))
      colnames(simulated_corr) = colnames(pvals)
      moduleTraitCor = cbind(moduleTraitCor, simulated_corr)
    }

    idx = which(round(moduleTraitPvalue,2)>pval)
    for (i in idx) {
      textMatrix[i] = NA
    }

    pdf("Results/TF.modules_metadata", width = width, height = height)
    par(mar = c(25, 15, 3, 3))
    labeledHeatmap(Matrix = moduleTraitCor,
                   xLabels = colnames(moduleTraitCor),
                   yLabels = rownames(moduleTraitCor),
                   ySymbols = rownames(moduleTraitCor),
                   colorLabels = FALSE,
                   colors = blueWhiteRed(50),
                   textMatrix = textMatrix,
                   setStdMargins = FALSE,
                   cex.text = 0.5,
                   zlim = c(-1,1),
                   main = paste0("Clinical associations with TFs modules\nOnly showing significant associations (pvalue < ", pval, ")"))
    dev.off()
  }

}





compute.survival.analysis = function(features, survival.data, time_unit, p.value = 0.05, thres = 0.5, max_factors = Inf) {
  n_features <- ncol(features)
  significant_combinations <- list() # To store significant feature combinations

  # Generate all possible combinations of the features
  contador = 1
  for (n in 1:min(n_features, max_factors)) {
    combinations <- combn(1:n_features, n, simplify = FALSE)

    for (comb in combinations) {
      # Create a formula dynamically based on the combination
      formula <- as.formula(paste("Surv(time, status) ~", paste(colnames(features)[comb], collapse = " + ")))

      # Prepare the data for survival analysis
      data_for_model <- data.frame("time" = survival.data$PFS,
                                   "status" = survival.data$DRP_st)

      data_for_model = cbind(data_for_model, features[,comb, drop=F])

      # Fit the Cox PH model with the combination of features (cox PH take into account covariates and measure the impact of each variable in the survival time)
      cox <- cph(formula, data = data_for_model)
      data_for_model$CoxPredictors <- cox$linear.predictors #linear predictors is the risk score for each individual in the dataset

      # Check that the model is significant as a predictor (maybe not useful?, it gives the same linear.predictos - to be check)
      cphmodel <- coxph(Surv(time, status) ~ CoxPredictors, data = data_for_model)
      data_for_model$CoxPredictors <- cphmodel$linear.predictors

      quantiles <- quantile(data_for_model$CoxPredictors, thres)

      # Binarize the Cox model output to draw two KM lines (linear predictors are used to stratify between high-risk and low-risk groups)
      data_for_model$coxHL <- ifelse(cphmodel$linear.predictors >= quantiles, 'High', "Low")

      # Perform Kaplan-Meier based on coxHL
      km_fit <- survfit(Surv(time, status) ~ coxHL, data = data_for_model)

      pval <- surv_pvalue(km_fit, data = data_for_model)$pval #Performs log-rank test to see whether both survival curves are significantly different

      if (!is.na(pval) && pval < p.value) {
        significant_combinations[[contador]] <- formula
        names(significant_combinations)[contador] = paste0("Formula_", contador)

        pdf(paste0("Results/SurvPlot_", names(significant_combinations)[contador]), width = 10, height = 5, onefile = FALSE)
        print(ggsurvplot(km_fit,
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

compute.cell.communication = function(counts, r2.thres=0.9, file_name){
  bsrdm <- prepareDataset(counts = counts, normalize = F, method = "log2transformed", log.transformed = T)
  #L-R interactions
  bsrdm <- learnParameters(bsrdm, plot.folder = "Results/", filename = "sdc", verbose = T)
  bsrinf <- initialInference(bsrdm)
  #Cell scores
  data(immune.signatures, package="BulkSignalR")
  data(tme.signatures, package="BulkSignalR")
  immune.signatures <- immune.signatures[!(immune.signatures$signature %in% c("Red pulp macrophages","Nuocytes","Megakaryocytes",
                                                                              "Eosinophils","Plasmacytoid dendritic cells","T follicular helper cells")),]
  signatures <- rbind(immune.signatures,tme.signatures[tme.signatures$signature%in%c("Endothelial cells","Fibroblasts"),])
  tme.scores <- scoreSignatures(bsrdm, signatures)
  #Assign cell types to interactions
  lr2ct <- assignCellTypesToInteractions(bsrdm, bsrinf, tme.scores)
  #Cellular network computation and plot
  g.table <- cellularNetworkTable(lr2ct)
  g.table.filt = g.table[g.table$r2>r2.thres&g.table$score>summary(g.table$score)[5],]
  gSummary <- summarizedCellularNetwork(g.table.filt)
  png(paste0("Results/CC_network_", file_name, ".png"), width = 800, height = 600)
  plot(gSummary, edge.width=1+30*E(gSummary)$score)
  dev.off()

  return(g.table.filt)
}

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
  extracted_names <- unique(na.omit(extracted_names))
  return(extracted_names)
}

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

remove.cell.groups.corr <- function(data, colors, threshold = 0.9) {

  features_high_corr = c()
  # Compute correlation matrix
  corr_matrix <- cor(data[[1]])
  # Find highly correlated features
  contador = 1
  while(nrow(corr_matrix)>0){
    color_features = c()
    feature = data.frame(corr_matrix[1, , drop = FALSE]) #Extract first row feature
    feature = feature %>%                                #Take only high corr above threshold
      mutate_all(~ifelse(. > threshold, ., NA)) %>%
      select_if(~all(!is.na(.)))

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
    column_to_rownames("Samples")

  cell.groups = remove.cell.groups.corr(cell.groups, module_colors, threshold = 0.9)

  cat("Removing low variance features (if present)........................\n")
  zero = nearZeroVar(cell.groups[[1]], saveMetrics = TRUE)
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

cell.groups.analysis = function(deconvolution, tfs.module.network, cell.dendrograms, cut.height, return = T, width = 12, height = 14){

  #If no different cuts have been defined
  if(length(cut.height)==1){
    cut.height = rep(cut.height, length(cell.dendrograms))
  }else if(length(cut.height)!=length(cell.dendrograms)){
    stop("Number of cuts specified don't match with number of dendrograms")
  }


  tfs.module.matrix = tfs.module.network[[1]]
  module_colors = unique(tfs.module.network[[2]])

  if(return){
    for (i in 1:length(cell.dendrograms)) {
      pdf(paste0("Results/Dendrogram_cell_types_", names(cell.dendrograms)[i], "_cut_", cut.height[i]), height = 20)
      par(mar = c(5, 2, 4, 35)) #bottom, left, top, right
      plot(as.dendrogram(cell.dendrograms[[i]]), horiz= T)
      dendextend::rect.dendrogram(as.dendrogram(cell.dendrograms[[i]]), h=cut.height[i], horiz=TRUE)
      dev.off()
    }
  }

  ###Identify cell groups
  cell.groups.dendrograms = list()
  cell.groups = list()
  exists = FALSE
  contador = 1
  for (i in 1:length(cell.dendrograms)) {
    x = list() #List for saving cell groups scores
    y = list() #List for saving cell names
    clusters <- dendextend::cutree(cell.dendrograms[[i]], h = cut.height[i], order_clusters_as_data=FALSE)
    for (j in 1:length(table(clusters))) {
      cells = names(clusters)[clusters==j]
      y[[j]] = cells
      ###Compute score for each cell group
      if(length(cells)>1){
        pca_group = deconvolution[,colnames(deconvolution) %in% cells, drop = F]
        color = extract_colors(module_colors, names(cell.dendrograms)[i])
        x[[j]] <- compute_composite_score(pca_group, color, tfs.module.matrix)
      }else{
        x[[j]] = deconvolution[,colnames(deconvolution) %in% cells]
      }
      names(x)[j] = paste0("group_", j)
      names(y)[j] = paste0("group_", j)
    }
    cell.groups[[contador]] = y
    cell.groups.dendrograms[[contador]] = data.frame(x)
    names(cell.groups.dendrograms)[contador] = names(cell.dendrograms)[contador]
    names(cell.groups)[contador] = names(cell.dendrograms)[contador]
    contador = contador + 1
  }

  ###################################################Remove groups composed of 1 feature
  for (i in 1:length(cell.groups)) {
    vec = c()
    for (j in 1:length(cell.groups[[i]])) {
      if(length(cell.groups[[i]][[j]])==1){
        vec = c(vec, j)
      }
    }
    if(length(vec)>0){
      cell.groups[[i]] = cell.groups[[i]][-vec]
      cell.groups.dendrograms[[i]] = cell.groups.dendrograms[[i]][-vec]
    }
  }

  ###################################################Naming cell groups
  cell.groups.dendrograms_all = list()
  cell.groups_all = list()
  contador = 1

  for (i in 1:length(cell.groups.dendrograms)) {
    for (j in 1:ncol(cell.groups.dendrograms[[i]])) {
      cell.groups.dendrograms_all[[contador]] = cell.groups.dendrograms[[i]][[j]]
      names(cell.groups.dendrograms_all)[contador] = paste0("Dendrogram_", names(cell.groups.dendrograms)[[i]], ".", names(cell.groups.dendrograms[[i]])[[j]])

      cell.groups_all[[contador]] = cell.groups[[i]][[j]]
      names(cell.groups_all)[contador] = paste0("Dendrogram_", names(cell.groups)[[i]], ".", names(cell.groups[[i]])[[j]])

      contador = contador + 1
    }
  }

  cell.groups.dendrograms_all = data.frame(cell.groups.dendrograms_all)
  rownames(cell.groups.dendrograms_all) = rownames(deconvolution)

  ###################################################Remove high correlated cell groups between each TF module
  groups = list(cell.groups.dendrograms_all, cell.groups_all)
  output = remove.cell.groups.corr(groups, unique(tfs.module.network[[2]]), threshold = 0.9)

  cell.groups.dendrograms_all = output[[1]]
  cell.groups_all = output[[2]]

  ###################################################Export clusters in a table and save
  clusters = data.frame(matrix(nrow = length(cell.groups_all), ncol = 2))
  for (i in 1:length(cell.groups_all)) {
    clusters[i,1] = names(cell.groups_all)[[i]]
    clusters[i,2] = paste(cell.groups_all[[i]], collapse ="\n")
  }
  colnames(clusters) = c("Cell groups", "Methods-signatures")

  write.csv(clusters, "Results/Cell.groups.names.csv", row.names = F)
  write.csv(cell.groups.dendrograms_all, "Results/Cell.groups.values.csv")


  return(list(data.frame(scale(cell.groups.dendrograms_all)), cell.groups_all))

}

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

find.communities = function(network, file_name){
  g <- graph_from_data_frame(network, directed=F, vertices=NULL)
  #Compute algorithms for communities finding
  k1 = cluster_walktrap(g)
  k2 = cluster_infomap(g)
  k3 = cluster_edge_betweenness(g)
  k4 = cluster_louvain(g)

  png(paste0("Results/CC_network_walktrap_", file_name, ".png"), width = 800, height = 600)
  plot(k1,g)
  dev.off()

  png(paste0("Results/CC_network_infomap_", file_name, ".png"), width = 800, height = 600)
  plot(k2,g)
  dev.off()

  png(paste0("Results/CC_network_edge_betweenness_", file_name, ".png"), width = 800, height = 600)
  plot(k3,g)
  dev.off()

  png(paste0("Results/CC_network_louvain_", file_name, ".png"), width = 800, height = 600)
  plot(k4,g)
  dev.off()
}

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
    cor_group = cor(aveg_group, PCs[,i], use = "p")

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

create_tfs_modules = function(TF.matrix, network_tfs){

  tfs.modules = TF.matrix %>%
    t() %>%
    data.frame() %>%
    mutate(Module = "na")

  for (i in 1:length(network_tfs[[3]])) {
    tfs.modules$Module[which(rownames(tfs.modules) %in% network_tfs[[3]][[i]])] = names(network_tfs[[3]])[i]
  }

  tfs_colors = tfs.modules %>%
    pull(Module)

  MEList = moduleEigengenes(TF.matrix, colors = tfs_colors, scale = F) #Data already scale
  MEs = MEList$eigengenes
  MEs = orderMEs(MEs)

  return(MEs)
}

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
      deconv_subgroups_values = cbind(deconv_subgroups_values, rowMedians(as.matrix(deconvolution_test[,base_groups[[i]]]))) #Compute median using base groups
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

CellTFusion = function(raw.counts, coldata = NULL, cbsx.mail = NULL, cbsx.token = NULL, file_name = NULL){

  #Normalize counts
  counts.norm = data.frame(ADImpute::NormalizeTPM(raw.counts, log = T))
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

  cat("Everything done! Results are saved in Results/ folder............................................................\n")

  res = list("Deconvolution" = deconv, "TFs_matrix" = tfs, "TF_network" = network, "Pathways_scores" = pathways,
             "Processed_deconvolution" = dt, "Cell_dendrograms" = cell_dendrograms, "Cell_groups" = cell.groups)

  return(res)

}
