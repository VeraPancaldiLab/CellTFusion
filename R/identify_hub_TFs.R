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
