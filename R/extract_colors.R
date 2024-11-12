#' Extract colors
#'
#' Extract TF module colors from cell type groups
#'
#' @param module_colors A character vector with the TF module colors. This can be found as an element from the output of compute.WTCNA()
#' @param cell_group_name Cell type group name
#'
#' @return A character vector with the module colors
#' @export
#'
#' @examples
#'
#' color = extract_colors(module_colors, cell.dendrogram.group.1)
#'
extract_colors <- function(module_colors, cell_group_name) {
  module_colors = c(module_colors, "all")
  matches <- c() # For storing the matches
  for (color in module_colors) {
    match <- regexpr(color, cell_group_name) # Find the position of the match
    if (match != -1) {
      matches <- c(matches, regmatches(cell_group_name, match))
    }
  }


  if (length(matches) > 0) {
    order <- sapply(matches, function(m) regexpr(m, cell_group_name)) # Sort matches based on their position in the original string to ensure names are the same
    matches <- matches[order(order)]  # Order the matches based on their position
    return(matches)  # Return the ordered matches
  } else {
    return(NA)  # If no matches are found, return NA
  }

}
