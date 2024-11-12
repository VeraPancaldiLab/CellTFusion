#' Compute metadata association between TF-modules and clinical data
#'
#' It tests if there is an association between the TF modules scores and the available clinical traits. It does a pearson correlation for quantitative traits and an anova test for the qualitative ones.
#'
#' @param tfs.modules A matrix of TF modules scores across samples obtained from compute.WTCNA()
#' @param coldata A matrix with the clinical data to test
#' @param pval A numeric value to set if a test it is significant or not (default is 0.05)
#' @param width A numeric value to set the width of the labeled heatmap to plot.
#' @param height A numeric value to set the height of the labeled heatmap to plot.
#'
#' @return A labeled heatmap with the pearson correlations and violin plots for the Anova tests are saved in the Results/ directory.
#' @export
#'
#' @examples
#'
#' compute.metada.association(tf_modules, traitData, pval = 0.05)
#'
compute.metada.association = function(tfs.modules, coldata, pval = 0.05, width = 20, height = 8){
  ###Association with categorical variables
  coldata_categorical = coldata %>%
    dplyr::select(where(is.character)|where(is.factor))

  if(ncol(coldata_categorical)!=0){
    data = cbind(tfs.modules, coldata_categorical)
    pvals = data.frame()
    fvals = data.frame()
    for(i in 1:ncol(tfs.modules)){
      contador = 1
      for (j in (ncol(tfs.modules)+1):ncol(data)) {
        module <- names(data[i])
        trait <- names(data[j])
        avz <- broom::tidy(aov(data[,i] ~ data[,j], data = data))
        pvals[i,contador] = avz$p.value[1]
        fvals[i,contador] = avz$statistic[1]
        contador = contador + 1
        if(avz$p.value[1] < pval) {
          pdf(paste0("Results/ANOVA_", module, "-", trait))
          print(ggplot(data, aes(x=data[,j], y=data[,i], fill=data[,j])) +
                  geom_violin(width=0.6) +
                  geom_boxplot(width=0.07, color="black", alpha=0.2) +
                  scale_fill_brewer() +
                  geom_smooth(aes(x=data[,j], y=data[,i]), method = "loess") +
                  ylab(paste0("Values for ", module)) +
                  xlab(paste0("Clinical trait: ", trait)) +
                  labs(title="One way ANOVA test",
                       subtitle=paste0("F statistic: ", round(avz$statistic[1],3), "\npvalue: ", round(avz$p.value[1], 3))) +
                  theme(axis.text.x = element_text(angle = 0),
                        axis.title.y = element_text(size = 8, angle = 90))+
                  scale_fill_discrete(name = trait))
          dev.off()
        }
      }
    }
    rownames(pvals) = colnames(tfs.modules)
    colnames(pvals) = colnames(coldata_categorical)

    rownames(fvals) = colnames(tfs.modules)
    colnames(fvals) = colnames(coldata_categorical)

    pvals = as.matrix(pvals)
    textMatrix2 = paste("ANOVA\n(", signif(pvals, 2), ")", sep = "")
    dim(textMatrix2) = dim(pvals)
  }

  ###Association with quantitative variables
  coldata_quantitative = coldata %>%
    dplyr::select(where(is.numeric))

  if(ncol(coldata_quantitative)!=0){
    moduleTraitCor = cor(tfs.modules, coldata_quantitative, method = "p");
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(tfs.modules))

    textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 2), ")", sep = "")
    dim(textMatrix) = dim(moduleTraitCor)

    if(ncol(coldata_categorical)!=0){
      textMatrix = cbind(textMatrix, textMatrix2)
      moduleTraitPvalue = cbind(moduleTraitPvalue, pvals)
      simulated_corr = matrix(runif(n=nrow(pvals)*ncol(pvals), min=-0.1, max=0.1), nrow = nrow(pvals), ncol = ncol(pvals))
      colnames(simulated_corr) = colnames(pvals)
      moduleTraitCor = cbind(moduleTraitCor, simulated_corr)
    }

    idx = which(round(moduleTraitPvalue,2)>pval)
    for (i in idx) {
      textMatrix[i] = NA
    }

    pdf("Results/TF.modules_metadata", width = width, height = height)
    par(mar = c(25, 15, 3, 3))
    labeledHeatmap(Matrix = moduleTraitCor,
                   xLabels = colnames(moduleTraitCor),
                   yLabels = rownames(moduleTraitCor),
                   ySymbols = rownames(moduleTraitCor),
                   colorLabels = FALSE,
                   colors = blueWhiteRed(50),
                   textMatrix = textMatrix,
                   setStdMargins = FALSE,
                   cex.text = 0.5,
                   zlim = c(-1,1),
                   main = paste0("Clinical associations with TFs modules\nOnly showing significant associations (pvalue < ", pval, ")"))
    dev.off()
  }

}
