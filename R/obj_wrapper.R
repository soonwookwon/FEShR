##' Create objective function (URE) for optimizing Lambda 
##'
##' Create objective function for optimizing Lambda that takes the centring term
##' as given
##' @param mu 
##' @param y 
##' @param M 
make_URE_obj <- function(mu, y, M) {

  T <- nrow(y)
  
  URE_obj <- function(L) {
    Lambda <- make_from_lowertri(L, T)
    return(URE(mu = mu, Lambda, y, M))
  }
  
  return(URE_obj)
}

##' Create objective function (URE) for optimizing Lambda 
##'
##' Create objective function for optimizing Lambda that takes the centring term
##' as given
##' @param mu 
##' @param y 
##' @param M 
make_URE_cpp_obj <- function(mu, y, M) {

  T <- nrow(y)
  
  URE_obj <- function(L) {
    Lambda <- make_from_lowertri(L, T)
    return(URE_cpp(mu = mu, Lambda, y, M))
  }
  
  return(URE_obj)
}


##' Create objective function (EBMLE) for optimizing Lambda 
##'
##' Create objective function for optimizing Lambda that takes the centring term
##' as given
##' @param mu 
##' @param y 
##' @param M 
make_nll_obj <- function(mu, y, M) {

  T <- nrow(y)

  nll_obj <- function(L) {
    Lambda <- make_from_lowertri(L, T)
    return(nll(mu = mu, Lambda, y, M))
  }
  
  return(nll_obj)
}
