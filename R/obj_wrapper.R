##' Create objective function (URE) for optimizing Lambda 
##'
##' Create objective function for optimizing Lambda that takes the centring term
##' as given
##' @param mu 
##' @param y 
##' @param M 
gen_obj_0 <- function(y, M, type) {

  T <- nrow(y)
  J <- ncol(y)

  if (type == "URE") {

    obj <- function(L) {
      Lambda <- make_from_lowertri(L, T)
      return(URE(mu = matrix(0, T, J), Lambda, y, M))
      
    }
  } else if (type == "EBMLE") {

    obj <- function(L) {
      Lambda <- make_from_lowertri(L, T)
      return(nll(mu = matrix(0, T, J), Lambda, y, M))    
    }
  }
  
  return(obj)
}

##' ...
##'
##' ..
##' @param y 
##' @param M 
##' @param type 
gen_obj_gen <- function(y, M, type, mu_abs_bounds) {

  T <- nrow(y)

  if (type == "URE") {
    
    obj <- function(L) {
      Lambda <- make_from_lowertri(L, T)
      return(URE(opt_mu_Lambda_URE(Lambda, y, M, mu_abs_bounds), Lambda, y, M))     
    }
    
  } else if (type == "EBMLE") {

    obj <- function(L) {
      Lambda <- make_from_lowertri(L, T)
      return(URE(opt_mu_Lambda_nll(Lambda, y, M), Lambda, y, M))     
    }
  }

  return(obj)
}

##' ...
##'
##' ..
##' @param Lambda 
##' @param y 
##' @param M 
##' @param mu_abs_bounds 
opt_mu_Lambda_URE <- function(Lambda, y, M, mu_abs_bounds)  {

  T <- nrow(y)
  J <- ncol(y)
  
  Dmat_dvec <- get_Dmat_dvec(Lambda, y, M)
  mu_opt <- quadprog::solve.QP(Dmat = Dmat_dvec$Dmat,
                               dvec = as.numeric(Dmat_dvec$dvec),
                               Amat = cbind(diag(T), - diag(T)),
                               bvec = - rep(mu_abs_bounds, 2))$solution
  
  return(matrix(mu_opt, T, J))
}
