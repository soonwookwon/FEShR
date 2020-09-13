##' @title Shrinkage estimator for jth observation
##' @description Computes the shrinkage estimator for theta_j
##'
##' ..
##' @param mu 
##' @param Lambda 
##' @param Mj 
##' @param yj 
thetahat_j <- function(mu, Lambda, Mj, yj) {

  if(mu == 0) mu <- rep(0, nrow(Lambda))
  shr_mat <- Lambda %*% solve (Lambda + Mj)
  thetahat <- (diag(ncol(Lambda)) - shr_mat) %*% mu + shr_mat %*% yj
  
  return(thetahat)
}

##' @title Shrinkage estimator for theta in the independent case
##'
##' ..
##' @param mu 
##' @param Lambda  
##' @param M JT vector of variances
##' @param y JT vector of observations
thetahat_diag <- function(mu = 0, Lambda, M, y) {

  shr_vec <- Lambda / (Lambda + M)
  thetahat <- (1 -shr_vec) * mu + shr_vec * y
  
  return(thetahat)
}
