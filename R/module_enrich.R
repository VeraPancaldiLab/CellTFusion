#' Module enrichment
#'
#' @param tpm.counts A matrix with normalized counts (genes as rows and samples as columns)
#' @param module_color A character vector with TF module colors.
#' @param hub_genes List of hub TFs per module.
#' @param tfs_universe A matrix with TF-gene interactions
#'
#' @return Reactome results
#' @export
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
