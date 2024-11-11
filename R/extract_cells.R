extract_cells = function(groups){
  names_cells = c("B.cells", "B.naive", "B.memory", "Macrophages.cells", "Macrophages.M0", "Macrophages.M1", "Macrophages.M2", "Monocytes", "Neutrophils", "NK.cells", "NK.activated",
                  "NK.resting", "NKT.cells", "CD4.cells", "CD4.memory.activated", "CD4.memory.resting", "CD4.naive", "CD8.cells", "T.cells.regulatory", "T.cells.non.regulatory","T.cells.helper",
                  "T.cells.gamma.delta", "Dendritic.cells", "Dendritic.activated", "Dendritic.resting", "Cancer", "Endothelial", "Eosinophils", "Plasma.cells", "Myocytes", "Fibroblasts",
                  "Mast.cells", "Mast.activated", "Mast.resting", "CAF")

  regex_pattern <-  paste0("(", paste(names_cells, collapse = "|"), ")")

  extracted_names <- sapply(groups, function(x) {
    match <- regexpr(regex_pattern, x)
    if (match != -1) {
      return(regmatches(x, match))
    } else {
      return(NA)
    }
  })

  extracted_names <- unname(extracted_names)
  extracted_names <- unique(na.omit(extracted_names))
  return(extracted_names)
}
