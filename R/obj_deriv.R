##' @title URE
##' @description Calculates URE(mu, Lambda) as defined in Kwon (2020)
##' 
##' @param y a T-by-J data matrix
##' @param M a length J list with the corresponding covariance matrices 
make_URE_deriv <- function(y, M) {

  T <- nrow(y)
  J <- ncol(y)
  
  URE_deriv <- function(Lambda_entries) {

    Lambda <- matrix(Lambda_entries, nrow = T)
    URE_deriv <- 0

    for (j in 1:J) {
      y_j <- y[, j]
      M_j <- M[[j]]
      URE_deriv <- URE_deriv + URE_deriv_j(Lambda, y_j, M_j)
    }
    
    return((1/J) * as.vector(t(URE_deriv))) 
  }
}

##' @title make_URE_j
##' @description Makes the function that calculate the component of URE(mu,
##'   Lambda) that cooresponds to the jth observation ..
##' @param mu 
##' @param Lambda 
##' @param y_j 
##' @param M_j

URE_deriv_j <- function(Lambda, y_j, M_j) {
  
  inv_Lam_Mj <- solve(Lambda + M_j)
  Mj_inv_Lam_Mj <- M_j %*% inv_Lam_Mj
  URE_deriv_j <- 2 * crossprod(Mj_inv_Lam_Mj)
  URE_deriv_j <- URE_deriv_j - inv_Lam_Mj %*% tcrossprod(y_j) %*% URE_deriv_j
  
  return(URE_deriv_j)
}

## @export
check_deriv <- function(y, M, Lambda_entries) {
  T <- nrow(y)
  J <- ncol(y)

  obj <- make_URE_obj(y, M, centering = "0", spg = TRUE)
  numer_gr <- numDeriv::grad(obj, Lambda_entries)
  obj_deriv <- make_URE_deriv(y, M)
  an_gr <- obj_deriv(Lambda_entries)

  return(list(matrix(numer_gr, nrow = T), matrix(an_gr, nrow = T)))
}
