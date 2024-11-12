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
