#' Unregister workers
#'
#' @return Clean parallelization
#' @export
#'
#' @examples
#'
#' unregister_dopar()
#'
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
  gc()
}
