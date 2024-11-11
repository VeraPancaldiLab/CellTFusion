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
