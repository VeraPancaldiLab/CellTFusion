#' Find maximum iteration from subgroups
#'
#' @param cells.groups Cell groups corresponding to a specific cell type.
#'
#' @return Maximum subgroupping iteration
#' @export
#'
#' @examples
#'
#' iterations = find.maximum.iteration(deconv_subgroups)
#'
find.maximum.iteration = function(cells.groups){
  max_iteration = c()
  for (i in 1:length(cells.groups)){
    if(is.null(names(cells.groups[[i]]))==F){
      iterations <- sapply(names(cells.groups[[i]]), function(x) {
        as.numeric(sub(".*\\.Iteration\\.(\\d+)", "\\1", x))
      })
      local_max = max(iterations)
      max_iteration = c(max_iteration, local_max)
    }
  }

  return(max(max_iteration))
}
