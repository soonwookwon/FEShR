##' @title Loss function of thetahat(mu, Lambda) at theta
##'
##' ..
##' @param mu 
##' @param Lambda 
##' @param y 
##' @param M 
##' @param T 
##' @param theta T-by-J matrix of true mean vector  
loss <- function(mu, Lambda, y, M, T, theta) {

  if (!is.matrix(y)) {
    y <- matrix(y, nrow = T)
  }

  J <- ncol(y_mat)
  
  loss <- 0
   
  for(j in 1:J) {
    diff <- thetahat_j(mu, Lambda, M[[j]], y[, j]) - theta[, j]
    loss <- loss + sum(diff^2)
  }
  
  return((loss / J))
}

##' @title Loss fucntion for the independent case
##'
##' ..
##' @param mu 
##' @param Lambda 
##' @param y 
##' @param M 
##' @param T 
##' @param theta 
loss_diag <- function(mu, Lambda, y, M, T, theta) {

  loss <- 0
  diff <- thetahat_diag(mu, Lambda, M, y) - theta
  
  return(sum(diff^2))
}

##' @title Loss function for general estimator
##'
##' ..
##' @param theta 
##' @param thetahat 
##' @param T 
loss_generic <- function(theta, thetahat, T) {

  if (!is.matrix(thetahat)) {
    y <- matrix(thetahat, nrow = T)
  }
  
  J <- ncol(thetahat)

  return(sum((theta - thetahat)^2) / J)
}
