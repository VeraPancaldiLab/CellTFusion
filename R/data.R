

#' Raw counts
#'
#' Raw gene expression matrix from bulk RNAseq.
#'
#' @format A data frame with 31086 genes as rows and 192 samples as columns
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
#' @format A matrix with 23205 genes as rows and 192 samples as columns
#'
#' @examples
#' data(counts.norm.tuto)
#' head(counts.norm.tuto)
"counts.norm.tuto"

#' Clinical data
#'
#' Data frame with the clinical data across samples
#'
#' @format A data frame with 192 samples as rows and 3 traits as columns
#'   (\code{Best.Confirmed.Overall.Response}, \code{binaryResponse} and \code{Enrollment.IC})
#'
#' @source Mariathasan et al. (2018), doi: https://doi.org/10.1038/nature25501
#'
#' @examples
#' data(traitdata.tuto)
#' head(traitdata.tuto)
"traitdata.tuto"

#' TFs data
#'
#' Data frame with the inferred TFs activity across samples
#'
#' @format A data frame with 192 samples as rows and 770 TFs as columns
#'
#' @examples
#' data(tfs.tuto)
#' head(tfs.tuto)
"tfs.tuto"

#' Cell subgroups
#'
#' Cell subgroups composition, as returned by \code{multideconv::compute.deconvolution.analysis()}
#' (multideconv >= 0.2.0, \code{corr = 0.7}) on \code{deconv.tuto}
#'
#' @format A list of 6 elements; the first one (\code{Deconvolution matrix}) is the processed
#'   deconvolution matrix (192 samples x 184 features) and the third one
#'   (\code{Deconvolution subgroups composition}) the composition of the cell subgroups
#'
#' @examples
#' data(deconv_subgroups.tuto)
#' deconv_subgroups.tuto[[1]]
"deconv_subgroups.tuto"

#' TF Network
#'
#' Network file obtained from compute.WTCNA()
#'
#' @format A list of 5 elements (\code{TFs module matrix}, \code{TFs colors}, \code{TFs per module},
#'   \code{Proportion of variance} and \code{TFs_matrix}), where the first element corresponds to the
#'   TF modules scores per sample (192 samples x 10 modules)
#'
#' @examples
#' data(network.tuto)
#' head(network.tuto)
"network.tuto"


#' Example Deconvolution Results
#'
#' A toy dataset included with the package to illustrate usage of
#' functions in **CellTFusion**. This dataset contains the output
#' of a deconvolution run (cell-type proportions estimated with different
#' methods and signatures) on the tutorial samples.
#'
#' @format A data frame with 192 samples as rows and 399 deconvolution features as columns
#'   (named \code{<method>_<signature>_<cell type>}, or \code{<method>_<cell type>} for methods
#'   with a single signature such as Quantiseq)
#'
#' @examples
#' data(deconv.tuto)
"deconv.tuto"
