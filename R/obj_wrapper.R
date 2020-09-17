make_URE_obj <- function(y, M, centering) {

  T <- nrow(y)

  if (centering == "0") {
    URE_obj <- function(L) {
      Lambda <- make_from_lowertri(L, T)
      return(URE(mu = 0, Lambda, y, M))
    }
  } else if (centering == "gen") {
    ## insert URE for general centering case
  } else {
    ##  insert URE for covariate case
  }
   
  return(URE_obj)
}

make_nll_obj <- function(y, M, centering) {

  T <- nrow(y)

  if (centering == "0") {
    nll_obj <- function(L) {
      Lambda <- make_from_lowertri(L, T)
      return(nll(mu = 0, Lambda, y, M))
    }
  } else if (centering == "gen") {
    ## insert URE for general centering case
  } else {
    ##  insert URE for covariate case
  }
   
  return(nll_obj)
}
