#' Perform pairwise correlation across all features
#'
#' @param data Matrix with features to correlate
#'
#' @return Dataframe containing all significant correlations (pvalue < 0.05)
#' @export
#'
#' @examples
#'
#' df = correlation(data)
#'
correlation <- function(data) {

  M <- Hmisc::rcorr(as.matrix(data), type = "pearson")
  Mdf <- map(M, ~data.frame(.x))

  corr_df = Mdf %>%
    map(~rownames_to_column(.x, var="measure1")) %>%
    map(~pivot_longer(.x, -measure1, names_to = "measure2")) %>%
    bind_rows(.id = "id") %>%
    pivot_wider(names_from = id, values_from = value) %>%
    dplyr::rename(p = P) %>%
    mutate(sig_p = ifelse(p < .05, T, F),
           p_if_sig = ifelse(sig_p, p, NA),
           r_if_sig = ifelse(sig_p, r, NA))

  corr_df = na.omit(corr_df)  #remove the ones that are the same TFs (pval = NA)
  corr_df <- corr_df[which(corr_df$sig_p==T),]  #remove not significant
  corr_df <- corr_df[order(corr_df$r, decreasing = T),]
  corr_df$AbsR =  abs(corr_df$r)

  return(corr_df)

}
