compute.boruta <- function(data, seed, fix = TRUE) {

  set.seed(seed)
  boruta_output <- Boruta(target ~ ., data = data, doTrace = 0)

  if (fix) {
    roughFixMod <- TentativeRoughFix(boruta_output)
    boruta_output <- roughFixMod
  }

  imps <- attStats(boruta_output)
  decision <- as.character(imps$decision)

  res <- imps %>%
    data.frame() %>%
    rownames_to_column("Variable") %>%
    dplyr::select(-decision)


  return(list(res, decision))
}
