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

  T <- nrow(y)
  J <- ncol(y)
 
  URE <- 0
  for (j in 1:J) {
    y_j <- y[, j]
    M_j <- M[[j]]
    URE <- URE + URE_j(mu, Lambda, y_j, M_j)
  }
  return((1/J) * URE)
}


##' @title URE with covariates
##' @description Calculates URE(Z gamma, Lambda) as defined in Kwon (2020)
##' 
##' @param gamma 
##' @param Lambda 
##' @param y a T-by-J data matrix
##' @param M a length J list with the corresponding covariance matrices 
##' @param Z 
URE_cov <- function(gamma, Lambda, y, M, Z) {

  if (!is.matrix(y)) {
    stop("y must be a matrix")
  }

  T <- nrow(y)
  J <- ncol(y)
 
  URE <- 0
  for (j in 1:J) {
    y_j <- y[, j]
    M_j <- M[[j]]
    URE <- URE + URE_j(Z[[j]] %*% gamma, Lambda, y_j, M_j)
  }
  return((1/J) * URE)
}

##' @title make_URE_j
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


##' @title URE_diag
##' @description Calculates URE(mu, Lambda) as defined in Xie, Kou, and Brown
##'   (2012), but allowed to have different lambdas periodwise
##' 
##' @param mu 
##' @param Lambda 
##' @param y 
##' @param M 
URE_diag <- function(mu, Lambda, y, M) {
  inv_Lam_M <- 1 / (Lambda + M)
  URE <- -2 * sum(inv_Lam_M * M^2) + sum(M * inv_Lam_M * (y - mu)^2)
  return((1/J) * URE)
}
