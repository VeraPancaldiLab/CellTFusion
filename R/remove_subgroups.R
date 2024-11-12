#' Remove subgroups that have the same method across different signatures
#'
#' @param groups Cell groups of features within cell types.
#'
#' @return List of position of groups which have features of same method.
#' @export
#'
#' @examples
#'
#' lis = remove_subgroups(subgroup)
#'
remove_subgroups = function(groups){
  lis = c()
  for (pos in 1:length(groups)){
    x = c()
    if(length(groups[[pos]])!=0){
      for (i in 1:length(groups[[pos]])) {
        x =  c(x,str_split(groups[[pos]][[i]], "_")[[1]][[1]])
      }
      if(length(unique(x)) == 1){
        lis = c(lis, pos)
      }
    }
  }

  return(lis)
}
