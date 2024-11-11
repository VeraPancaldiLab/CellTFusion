#' Fisher test using cell groups scores
#'
#' @param cell.groups A list of cell groups obtained from cell.groups.computation()
#' @param coldata A data frame containing the trait column to test with Fisher
#' @param trait Name of the column to test from coldata
#' @param pval pvalue to consider significant (default 0.05)
#'
#' @return A list containing the significant cell groups after Fisher test. Additionally it saves in the Results/ folder its respective plots
#' @export
#'
#' @examples
#'
#' res_fisher = cell.groups.fisher.test(cell.groups, clinical.data, "Response")
#'
cell.groups.fisher.test = function(cell.groups, coldata, trait, pval = 0.05){

  sig = c()
  coldata[,trait] = as.factor(coldata[,trait])

  for (j in 1:ncol(cell.groups[[1]])) {
    coldata = coldata %>%
      mutate(level = cell.groups[[1]][,j],
             Cells_level = ifelse(level > summary(level)[3], 'High', 'Low'))

    contingency = table(coldata[,"Cells_level"], coldata[,trait])
    test = fisher.test(contingency)

    ##Extract only significant features
    if(round(test$p.value, 5) <= pval){
      cat("Significant pval after doing Fisher test for", colnames(cell.groups[[1]])[j], "\n")
      df = data.frame("Cells_level" = coldata[,"Cells_level"], "Trait" = coldata[,trait])
      pdf(paste0("Results/Fisher_", trait, "_", colnames(cell.groups[[1]])[j]), width = 12, height = 9)
      print(ggbarstats(df, Cells_level, Trait, results.subtitle = F,
                       title= paste0("Dendrogram_", colnames(cell.groups[[1]])[j]),
                       subtitle = paste0("Fisher's exact test, p-value = ", ifelse(test$p.value < 0.001, "< 0.001", round(test$p.value, 5))))+
              ggplot2::theme(plot.title = ggplot2::element_text(size=15), axis.text = ggplot2::element_text(size=14), legend.title = ggplot2::element_text(size=14)))
      dev.off()

      sig = c(j, sig)
    }
  }

  cell.groups.sig = list()
  cell.groups.sig[[1]] = cell.groups[[1]][,sig]
  cell.groups.sig[[2]] = cell.groups[[2]][sig]

  if(length(sig)==0){
    message("No significant cell groups (pvalue < ", pval, ") after Fisher test")
  }else{
    return(cell.groups.sig)
  }

}
