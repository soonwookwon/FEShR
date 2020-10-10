##' ...
##'
##' ..
##' @param L 
##' @param T 
make_from_lowertri <- function(L, T) {
  L_mat <- matrix(0, nrow = T, ncol = T)
  L_mat[lower.tri(L_mat, diag = TRUE)] <- L
  return(tcrossprod(L_mat))
}

##' ...
##'
##' ..
##' @param Lambda 
extract_lowertri <- function(Lambda) {
   return(as.vector(Lambda[lower.tri(Lambda, diag = TRUE)]))
}
