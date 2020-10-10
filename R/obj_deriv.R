##' @title URE
##' @description Calculates URE(mu, Lambda) as defined in Kwon (2020)
##' 
##' @param mu 
##' @param y a T-by-J data matrix
##' @param M a length J list with the corresponding covariance matrices 
gen_deriv_0 <- function(y, M, type) {

  T <- nrow(y)
  J <- ncol(y)

  if (type == "URE") {
    
    deriv_lt <- function(L) {

      Lambda <- make_from_lowertri(L = L, T = T)
      URE_deriv <- URE_deriv(matrix(0, T, J), Lambda, y, M)

      L_mat <- matrix(0, nrow = T, ncol = T)
      L_mat[lower.tri(L_mat, diag = TRUE)] <- L
      URE_deriv_lt <- (URE_deriv + t(URE_deriv)) %*% L_mat 
      URE_deriv_lt <- URE_deriv_lt[lower.tri(URE_deriv_lt, diag = TRUE)]
      
      return((1/J) * as.vector(URE_deriv_lt)) 
    }  
  } else if (type == "EBMLE") {
    
    deriv_lt <- function(L) {

      Lambda <- make_from_lowertri(L = L, T = T)
      nll_deriv <- nll_deriv(matrix(0, T, J), Lambda, y, M)

      L_mat <- matrix(0, nrow = T, ncol = T)
      L_mat[lower.tri(L_mat, diag = TRUE)] <- L
      nll_deriv_lt <- (nll_deriv + t(nll_deriv)) %*% L_mat 
      nll_deriv_lt <- nll_deriv_lt[lower.tri(nll_deriv_lt, diag = TRUE)]
      return((1/J) * as.vector(nll_deriv_lt)) 
    }
  }

  return(deriv_lt)
}


##' Check gradients are correct
##'
##' Function to make sure the numerical and anlytical gradients indeed agree
##' @param thetas 
##' @param mu 
##' @param Lambda_lt 
##' @param y a T-by-J data matrix
##' @param M a length J list with the corresponding covariance matrices 
##' @param method 
check_deriv <- function(thetas = NULL,  mu = 0, Lambda_lt, y, M, method) {
  T <- nrow(y)
  J <- ncol(y)
  
  if (method == "URE") {
    
    obj <- make_URE_cpp_obj(mu = matrix(0, T, J), y, M)
    numer_gr <- numDeriv::grad(obj, Lambda_lt)
    obj_deriv <- make_URE_cpp_deriv_lt(y, M)
    an_gr <- obj_deriv(Lambda_lt)
    
  } else if (method == "EBMLE") {
    
    obj <- make_nll_cpp_obj(mu = matrix(0, T, J), y, M)
    numer_gr <- numDeriv::grad(obj, Lambda_lt)
    obj_deriv <- make_nll_cpp_deriv_lt(y, M)
    an_gr <- obj_deriv(Lambda_lt)
    
  } else if (method == "ol") {
    
    obj <- make_oracle_cpp_obj(thetas, mu = matrix(0, T, J), y, M, diag = FALSE)
    numer_gr <- numDeriv::grad(obj, Lambda_lt)
    obj_deriv <- make_ol_cpp_deriv_lt(y, M, thetas)
    an_gr <- obj_deriv(Lambda_lt)
    
  }
  
  return(data.frame(numerical = numer_gr, analytical = an_gr))
}
