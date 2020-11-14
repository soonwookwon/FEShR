##' Create objective function (URE) for optimizing Lambda 
##'
##' Create objective function for optimizing Lambda that takes the centring term
##' as given
##' @param y 
##' @param M 
##' @param type 
##' @param W 
##' @param missing_cells 
##' @param O 
##' @param centering 
##' @param mu_abs_bounds 
gen_obj <- function(y, M, type, W, missing_cells, O, diag_lam, centering,
                    bounds = NULL, Z = list()) {

  T <- nrow(y)
  J <- ncol(y)

  if (!diag_lam) {
    if (centering == "0") {
      if (type == "URE") {

        obj <- function(L) {
          Lambda <- make_from_lowertri(L, T)
          return(URE(mu = matrix(0, T, J), Lambda, y, M, W, missing_cells, O))
        }
        
      } else if (type == "EBMLE") {

        obj <- function(L) {
          Lambda <- make_from_lowertri(L, T)
          return(nll(mu = matrix(0, T, J), Lambda, y, M, missing_cells, O))    
        }
        
      }
    } else if (centering == "gen" | centering == "cov") {

      if (type == "URE") {
        
        obj <- function(L) {
          Lambda <- make_from_lowertri(L, T)
          opt_mu <- opt_mu_Lambda_URE(Lambda, y, M, bounds, W,
                                      missing_cells, O, Z)
          return(URE(opt_mu, Lambda, y, M, W, missing_cells, O))     
        }
        
      } else if (type == "EBMLE") {

        obj <- function(L) {
          Lambda <- make_from_lowertri(L, T)
          opt_mu <- opt_mu_Lambda_nll(Lambda, y, M, missing_cells, O, Z)

          if (centering == "cov") {
            opt_mu <- gam_to_mu(opt_mu, Z, T)
          }
          
          return(nll(opt_mu, Lambda, y, M, missing_cells, O))     
        }
      }
    } 
  } else {
    if (centering == "0") {
      if (type == "URE") {

        obj <- function(D) {
          Lambda <- diag(D)^2
          return(URE(mu = matrix(0, T, J), Lambda, y, M, W, missing_cells, O))
        }
        
      } else if (type == "EBMLE") {

        obj <- function(D) {
          Lambda <- diag(D)^2
          return(nll(mu = matrix(0, T, J), Lambda, y, M, missing_cells, O))    
        }
        
      }
    } else if (centering == "gen" | centering == "cov") {

      if (type == "URE") {
        
        obj <- function(D) {
          Lambda <- diag(D)^2
          opt_mu <- opt_mu_Lambda_URE(Lambda, y, M, bounds, W,
                                      missing_cells, O, Z)
          return(URE(opt_mu, Lambda, y, M, W, missing_cells, O))     
        }
        
      } else if (type == "EBMLE") {

        obj <- function(D) {
          Lambda <- diag(D)^2
          opt_mu <- opt_mu_Lambda_nll(Lambda, y, M, missing_cells, O, Z)
          return(nll(opt_mu, Lambda, y, M, missing_cells, O))     
        }
      }

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
##' @param bounds 
##' @param W 
##' @param missing_cells 
##' @param O 
##' @param Z 
opt_mu_Lambda_URE <- function(Lambda, y, M, bounds, W, missing_cells, O,
                              Z = list()) {

  T <- nrow(y)
  J <- ncol(y)

  
  Dmat_dvec <- get_Dmat_dvec(Lambda, y, M, W, missing_cells, O, Z)

  if (length(Z) == 0) {
    mu_opt <- quadprog::solve.QP(Dmat = Dmat_dvec$Dmat,
                                 dvec = as.numeric(Dmat_dvec$dvec),
                                 Amat = cbind(diag(T), - diag(T)),
                                 bvec = - rep(bounds, 2))$solution
  } else {
    mu_opt <- solve(Dmat_dvec$Dmat, Dmat_dvec$dvec)
    mu_opt <- gam_to_mu(mu_opt, Z, T)
  }
  
  return(matrix(mu_opt, T, J))
}
