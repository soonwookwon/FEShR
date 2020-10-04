##' @title URE
##' @description Calculates URE(mu, Lambda) as defined in Kwon (2020)
##' 
##' @param y a T-by-J data matrix
##' @param M a length J list with the corresponding covariance matrices 
make_URE_deriv <- function(y, M, mu = 0) {

  T <- nrow(y)
  J <- ncol(y)
  
  URE_deriv <- function(Lambda_entries) {

    Lambda <- matrix(Lambda_entries, nrow = T)
    URE_deriv <- 0

    for (j in 1:J) {
      y_j <- y[, j]
      M_j <- M[[j]]
      URE_deriv <- URE_deriv + URE_deriv_j(Lambda, y_j, M_j, mu = mu)
    }
    
    return((1/J) * as.vector(t(URE_deriv))) 
  }
}

##' @title URE
##' @description Calculates URE(mu, Lambda) as defined in Kwon (2020)
##' 
##' @param y a T-by-J data matrix
##' @param M a length J list with the corresponding covariance matrices 
make_URE_deriv_lt <- function(y, M, mu = 0) {

  T <- nrow(y)
  J <- ncol(y)

  URE_deriv_lt <- function(L) {

    Lambda <- make_from_lowertri(L, T)
    URE_deriv <- 0

    for (j in 1:J) {
      y_j <- y[, j]
      M_j <- M[[j]]
      URE_deriv <- URE_deriv + URE_deriv_j(Lambda, y_j, M_j, mu = mu)
    }

    L_mat <- matrix(0, nrow = T, ncol = T)
    L_mat[lower.tri(L_mat, diag = TRUE)] <- L
    URE_deriv_lt <- (URE_deriv + t(URE_deriv)) %*% L_mat 
    URE_deriv_lt <- URE_deriv_lt[lower.tri(URE_deriv_lt, diag = TRUE)]
    return((1/J) * as.vector(URE_deriv_lt)) 
  }
}


##' @title Derivative of URE
##' @description Makes the function that calculate the component of URE(mu,
##'   Lambda) that cooresponds to the jth observation ..
##' @param Lambda 
##' @param y_j 
##' @param M_j 
URE_deriv_j <- function(Lambda, y_j, M_j, mu = 0) {
  
  inv_Lam_Mj <- chol2inv(chol(Lambda + M_j))
  Mj_inv_Lam_Mj <- M_j %*% inv_Lam_Mj
  URE_deriv_j <- 2 * crossprod(Mj_inv_Lam_Mj)
  URE_deriv_j <- URE_deriv_j - inv_Lam_Mj %*% tcrossprod(y_j - mu) %*% URE_deriv_j
  
  return(URE_deriv_j)
}


##' @title URE
##' @description Calculates URE(mu, Lambda) as defined in Kwon (2020)
##' 
##' @param y a T-by-J data matrix
##' @param M a length J list with the corresponding covariance matrices 
make_nll_deriv <- function(y, M, mu = 0) {

  T <- nrow(y)
  J <- ncol(y)
  
  nll_deriv <- function(Lambda_entries) {

    Lambda <- matrix(Lambda_entries, nrow = T)
    nll_deriv <- 0

    for (j in 1:J) {
      y_j <- y[, j]
      M_j <- M[[j]]
      nll_deriv <- nll_deriv + nll_deriv_j(Lambda, y_j, M_j, mu = mu)
    }
    
    return((1/J) * as.vector(t(nll_deriv))) 
  }
}

##' @title URE
##' @description Calculates URE(mu, Lambda) as defined in Kwon (2020)
##' 
##' @param y a T-by-J data matrix
##' @param M a length J list with the corresponding covariance matrices 
make_nll_deriv_lt <- function(y, M, mu = 0) {

  T <- nrow(y)
  J <- ncol(y)
  
  nll_deriv_lt <- function(L) {

    Lambda <- make_from_lowertri(L, T)
    nll_deriv <- 0

    for (j in 1:J) {
      y_j <- y[, j]
      M_j <- M[[j]]
      nll_deriv <- nll_deriv + nll_deriv_j(Lambda, y_j, M_j, mu = mu)
    }

    L_mat <- matrix(0, nrow = T, ncol = T)
    L_mat[lower.tri(L_mat, diag = TRUE)] <- L
    nll_deriv_lt <- (nll_deriv + t(nll_deriv)) %*% L_mat 
    nll_deriv_lt <- nll_deriv_lt[lower.tri(nll_deriv_lt, diag = TRUE)]
    return((1/J) * as.vector(nll_deriv_lt)) 
  }
}


##' @title Derivative of URE
##' @description Makes the function that calculate the component of URE(mu,
##'   Lambda) that cooresponds to the jth observation ..
##' @param Lambda 
##' @param y_j 
##' @param M_j 
nll_deriv_j <- function(Lambda, y_j, M_j, mu) {
  
  inv_Lam_Mj <- chol2inv(chol(Lambda + M_j))
  nll_deriv_j <- inv_Lam_Mj - tcrossprod(inv_Lam_Mj %*% (y_j - mu))
  
  return(nll_deriv_j)
}

##'@export
check_deriv <- function(y, M, mu = 0, Lambda_lt, thetas = NULL, method) {
  T <- nrow(y)
  J <- ncol(y)
  
  if (method == "URE") {
    obj <- make_URE_obj(y, M, centering = "0")
    numer_gr <- numDeriv::grad(obj, Lambda_lt)
    obj_deriv <- make_URE_deriv_lt(y, M)
    an_gr <- obj_deriv(Lambda_lt)
  } else if (method == "EBMLE") {
    obj <- make_nll_obj(y, M, centering = "0")
    numer_gr <- numDeriv::grad(obj, Lambda_lt)
    obj_deriv <- make_nll_deriv_lt(y, M)
    an_gr <- obj_deriv(Lambda_lt)
  } else if (method == "ol") {
    obj <- make_oracle_obj(thetas, y, M, centering = "0", diag = FALSE)
    numer_gr <- numDeriv::grad(obj, Lambda_lt)
    obj_deriv <- make_ol_deriv_lt(y, M, thetas)
    an_gr <- obj_deriv(Lambda_lt)
  }
  

  return(data.frame(numerical = numer_gr, analytical = an_gr))
}
