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
  shr_mat <- Lambda %*% chol2inv(chol(Lambda + Mj))
  thetahat <- (diag(ncol(Lambda)) - shr_mat) %*% mu + shr_mat %*% yj
  
  return(thetahat)
}
