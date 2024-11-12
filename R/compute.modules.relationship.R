#' Compute modules relationship
#'
#' Performs linear correlation between two matrices of features across samples
#'
#' @param tfs_network TFs module matrix (samples x modules).
#' @param matB matrix B to correlate (samples x features).
#' @param file_name string indicating the name of the figure to save.
#' @param width integer indicating the width of the figure, default is 8.
#' @param height integer indicating the height of the figure, default is 8.
#' @param pval p value to use as threshold to differentiate between significant and no significant variables. Default is 0.05.
#' @param padj logical value indicating correction of p-value by Bonferroni method has to be applied. Default is false.
#' @param cor_type type of correlation to be used. Default is pearson “p”.
#' @param return logical value indicating if significant features have to be returned. Default is False.
#' @param vertical logical value indicating if function should return a horizontal or vertical plot (this can be useful when looking at several deconvolution features).
#' @param plot whether to saved or not the plots.
#'
#' @return
#'
#' A heatmap showing the level of correlation between TFs modules and corresponding features. Only significant features are being shown. If return = TRUE it will return a list with two elements: the correlation matrix between the TFs modules and the other features and a character vector containing the names of the significant associated features. Note that features not significantly associated with any module are not returned.
#' @export
#'
#' @examples
#'
#' pathways = compute.pathway.activity(counts.norm)
#' compute.modules.relationship(tfs_modules, pathways, "Pathways_Progeny-TFs_Modules", width = 15)
#' dt = compute.deconvolution.analysis(deconv, corr = 0.7, seed = 123)
#' compute.modules.relationship(tfs_modules, deconvolution_matrix, "Deconvolution-TFs_Modules", vertical = T, height = 30, width = 10, pval = 0.05)
#'
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
