remove_low_variance <- function(data, plot = T) {
  vars <- apply(data, 2, var)
  threshold = summary(vars)[[2]]
  cat("Checking distribution of variances...............................................................\n\n")
  cat("Chosen threshold is:", threshold, "\n\n")
  cat("Saving variance distribution plot in Results/ folder\n\n")
  low_variance <- which(vars < threshold)
  cat("Removing", length(low_variance), "features with variance across samples below this threshold...............................................................\n\n")

  if(plot){
    pdf("Results/Distribution_variances_deconvolution.pdf")
    hist(vars, main = "Distribution of deconvolution variances across samples\nRemoving features below threshold (low variance)", xlab = "Variance", col = "skyblue", border = "white", xlim = range(vars))
    lines(density(vars), col = "red", lwd = 2)
    legend("topright", legend = c("Density", paste("Threshold =", round(threshold, 5))), col = c("red", "orange"), lty = c(1, 2), lwd = c(2, 2))

    # Shade region below threshold
    abline(v = threshold, col = "orange", lwd = 2, lty = 2)
    x <- density(vars)$x
    y <- density(vars)$y
    polygon(c(min(x[vars < threshold]), x[vars < threshold], max(x[vars < threshold])),
            c(0, y[vars < threshold], 0), col = adjustcolor("orange", alpha.f = 0.3), border = NA)
    dev.off()
  }

  data_filt = data[,-low_variance, drop = F]
  low_var_features = data[,low_variance, drop = F]

  res = list(data_filt, low_var_features)
  return(res)
}
