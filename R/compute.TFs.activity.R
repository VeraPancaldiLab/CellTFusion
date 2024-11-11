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
