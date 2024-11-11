cell.groups.anova.test = function(cell.groups, coldata, trait, pval = 0.05){

  sig = c()
  coldata[,trait] = as.factor(coldata[,trait])

  for (j in 1:ncol(cell.groups[[1]])) {
    data = data.frame("Value" = cell.groups[[1]][,j], "Trait" = coldata[,trait])
    model  <- lm(Value ~ Trait, data = data)
    res.aov <- data %>% rstatix::anova_test(Value ~ Trait)

    ##Extract only significant features
    if(round(res.aov$p, 5) <= pval){
      cat("Significant pval after doing Anova test for", colnames(cell.groups[[1]])[j], "\n")
      pdf(paste0("Results/Anova_", trait, "_", colnames(cell.groups[[1]])[j]), width = 12, height = 9)
      print(ggplot(data, aes(x=Trait, y=Value, fill=Trait)) +
              geom_violin(width=0.6) +
              geom_boxplot(width=0.07, color="black", alpha=0.2) +
              scale_fill_brewer() +
              geom_smooth(aes(x=Trait, y=Value), method = "loess") +
              xlab(paste0("Clinical trait: ", trait)) +
              labs(title= paste0("Dendrogram_", colnames(cell.groups[[1]])[j]),
                   subtitle = rstatix::get_test_label(res.aov, detailed = TRUE)) +
              theme(axis.text.x = element_text(angle = 0),
                    axis.title.y = element_text(size = 8, angle = 90)))
      dev.off()

      sig = c(sig, j)
    }
  }

  cell.groups.sig = list()
  cell.groups.sig[[1]] = cell.groups[[1]][,sig]
  cell.groups.sig[[2]] = cell.groups[[2]][sig]

  if(length(sig)==0){
    message("No significant cell groups (pvalue < ", pval, ") after Anova test")
  }else{
    return(cell.groups.sig)
  }

}
