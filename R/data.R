

#' Raw counts
#'
#' Raw gene expression matrix from bulk RNAseq.
#'
#' @format Matrix with genes as rows and samples as columns
#'
#' @source Mariathasan et al. (2018), doi: https://doi.org/10.1038/nature25501
#'
#' @examples
#' data(raw.counts)
#' head(raw.counts)
"raw.counts"

#' Log(TPM+1) normalized counts
#'
#' Normalized gene expression matrix from bulk RNAseq.
#'
#' @format Matrix with genes as rows and samples as columns
#'
#' @examples
#' data(counts.norm)
#' head(counts.norm)
"counts.norm"

#' Clinical data
#'
#' Data frame with the clinical data across samples
#'
#' @format Matrix with samples as rows and traits as columns
#'
#' @source Mariathasan et al. (2018), doi: https://doi.org/10.1038/nature25501
#'
#' @examples
#' data(traitdata)
#' head(traitdata)
"traitdata"

#' TFs data
#'
#' Data frame with the inferred TFs activity across samples
#'
#' @format Matrix with samples as rows and TFs as columns
#'
#' @examples
#' data(tfs)
#' head(tfs)
"tfs"

#' Cell subgroups
#'
#' Cell subgroups composition
#'
#' @format A list with the cell subgroups
#'
#' @examples
#' data(deconv_subgroups)
#' deconv_subgroups[[1]]
"deconv_subgroups"

#' TF Network
#'
#' Network file obtained from compute.WTCNA()
#'
#' @format List where first element corresponds to the TF modules scores per sample
#'
#' @examples
#' data(network)
#' head(network)
"network"
