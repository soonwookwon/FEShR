##' @title Negative log likelihood of the marginal distribution of the data
##' 
##' @param mu 
##' @param Lambda 
##' @param y a T-by-J data matrix
##' @param M a length J list with the corresponding covariance matrices 
nll <- function(mu, Lambda, y, M) {
  if (!is.matrix(y)) {
    stop("y must be a matrix")
  }

  T <- nrow(y)
  J <- ncol(y)
 
  nll <- 0
  
  for (j in 1:J) {
    y_j <- y[, j]
    M_j <- M[[j]]
    nll <- nll + log(det(Lambda + M_j)) +
      t(y_j - mu) %*% solve(Lambda + M_j) %*% (y_j - mu)
  }
  
  return((1/J) * nll)
}
