##' @title URE
##' @description Calculates URE(mu, Lambda) as defined in Kwon (2020)
##' 
##' @param mu 
##' @param Lambda 
##' @param y a T-by-J data matrix
##' @param M a length J list with the corresponding covariance matrices 
URE <- function(mu, Lambda, y, M) {

  if (!is.matrix(y)) {
    stop("y must be a matrix")
  }

  if (mu == 0) {
    mu <- matrix(0, nrow = nrow(y), ncol = ncol(y))
  } else if (length(mu) == ncol(y)) {
    mu <- matrix(rep(mu, nrow(y)), nrow = nrow(y), ncol = ncol(y))
  }
  
  T <- nrow(y)
  J <- ncol(y)
 
  URE <- 0
  for (j in 1:J) {
    mu_j <- mu[, j]
    y_j <- y[, j]
    M_j <- M[[j]]
    mu_j <- mu[, j]
    URE <- URE + URE_j(mu_j, Lambda, y_j, M_j)
  }
  return((1/J) * URE)
}


##' @title URE_j
##' @description Makes the function that calculate the component of URE(mu,
##'   Lambda) that cooresponds to the jth observation ..
##' @param mu 
##' @param Lambda 
##' @param y_j 
##' @param M_j 
URE_j <- function(mu, Lambda, y_j, M_j) {

  inv_Lam_Mj <- chol2inv(chol(Lambda + M_j))
  URE_j <- -2 * sum(diag(inv_Lam_Mj %*% M_j %*% M_j))
  URE_j <- URE_j + sum((M_j %*% inv_Lam_Mj %*% (y_j - mu))^2)

  return(URE_j)
}
