#' Compute TF activity
#'
#' The function infers TFs activity based on a gene expression matrix using the VIPER algorithm from [1] and a collection of TF-genes interactions (CollecTRI or Dorothea) from the OmnipathR package [2]. Users can also input TFs-gene interactions from ARACNE [4,5].
#'
#' @param RNA.counts Gene expression matrix with genes as rows and samples as columns, normalized either by TPM or other type of normalization.
#' @param TF.collection TF-gene interactions collection matrix used to infer the TFs activity. Default argument is “CollecTRI”, other arguments accepted are “Dorothea” [3] or “ARACNE” [4,5] for input your own TF-gene collection.
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
#' @examples
#'
#' tfs = compute.TFs.activity(counts.normalized)
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
