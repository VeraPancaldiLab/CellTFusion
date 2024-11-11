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

  res = list("Deconvolution" = deconv, "TFs_matrix" = tfs, "TF_network" = network, "Pathways_scores" = pathways,
             "Processed_deconvolution" = dt, "Cell_dendrograms" = cell_dendrograms, "Cell_groups" = cell.groups)

  return(res)

}
