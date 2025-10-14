

#' Raw counts
#'
#' Raw gene expression matrix from bulk RNAseq.
#'
#' @format Matrix with genes as rows and samples as columns
#'
#' @source Mariathasan et al. (2018), doi: https://doi.org/10.1038/nature25501
#'
#' @examples
#' data(raw.counts.tuto)
#' head(raw.counts.tuto)
"raw.counts.tuto"

#' Log(TPM+1) normalized counts
#'
#' Normalized gene expression matrix from bulk RNAseq.
#'
#' @format Matrix with genes as rows and samples as columns
#'
#' @examples
#' data(counts.norm.tuto)
#' head(counts.norm.tuto)
"counts.norm.tuto"

#' Clinical data
#'
#' Data frame with the clinical data across samples
#'
#' @format Matrix with samples as rows and traits as columns
#'
#' @source Mariathasan et al. (2018), doi: https://doi.org/10.1038/nature25501
#'
#' @examples
#' data(traitdata_tuto)
#' head(traitdata_tuto)
"traitdata_tuto"

#' TFs data
#'
#' Data frame with the inferred TFs activity across samples
#'
#' @format Matrix with samples as rows and TFs as columns
#'
#' @examples
#' data(tfs.tuto)
#' head(tfs.tuto)
"tfs.tuto"

#' Cell subgroups
#'
#' Cell subgroups composition
#'
#' @format A list with the cell subgroups
#'
#' @examples
#' data(deconv_subgroups.tuto)
#' deconv_subgroups.tuto[[1]]
"deconv_subgroups.tuto"

#' TF Network
#'
#' Network file obtained from compute.WTCNA()
#'
#' @format List where first element corresponds to the TF modules scores per sample
#'
#' @examples
#' data(network.tuto)
#' head(network.tuto)
"network.tuto"


#' Example Deconvolution Results
#'
#' A toy dataset included with the package to illustrate usage of
#' functions in **CellTFusion**. This dataset contains the output
#' of a deconvolution run for a small subset of genes.
#'
#' @format A data frame with X rows and Y variables:
#'
#' @examples
#' data(deconv.tuto)
"deconv.tuto"
