#' Remove cell groups with only one feature
#'
#' @param cell.values Cell groups scores
#' @param cell.composition Cell groups composition
#'
#' @return A list containing
#'
#' - Cell groups scores after removal of single cell groups
#' - Cell groups composition after removal of single cell groups
#'
#' @export
#'
#' @examples
#'
#' cell.groups = remove_single_groups(cell.values, cell.composition)
#'
remove_single_groups = function(cell.values, cell.composition){

  message("Removing cell groups composed of one single feature..............................................................................")
  vec = c()
  for (i in 1:length(cell.composition)) {
    if(length(cell.composition[[i]])==1){
      vec = c(vec, i)
    }
  }

  if(length(vec)>0){
    cell.composition = cell.composition[-vec]
    cell.values = cell.values[-vec]
  }

  if(length(cell.composition)==0){
    return(NULL)
  }else{
    return(list(cell.values, cell.composition))
  }

}
