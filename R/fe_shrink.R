##' Calculate the optiamal fixed effects estimator
##'
##' ..
##' @param y a T-by-J data matrix, possibly correlated
##'   within each column. If a vector is provided, it is assumed that T=1 and
##'   will be coerced to a 1-by-length(y) matrix.
##' @param M length J list of the covaraince matrices.
##' @param centering Centering
##' @param type 
##' @export
fe_shrink <- function(y, M, centering = c("0", "gen", "cov"), W = NULL,
                      type = c("URE", "EBMLE"), init_vals = 1) {
  
  if(!is.matrix(y)) {
    warning("y is not a matrix, assuming T=1.")
    y <- matrix(y, nrow = 1)
  }

  centering <- match.arg(centering)
  type <- match.arg(type)

  T <- nrow(y)
  J <- ncol(y)

  if (is.null(W)) W <- diag(T) # W is always diagonal for now
  
  thetahat <- matrix(0, nrow = T, ncol = J)
  
  if (centering == "0") {
    
    if (type == "URE") {
      obj <- make_URE_obj(y, M, centering, spg = spg)
    } else {
      obj <- make_nll_obj(y, M, centering, spg = spg)
    }
    
    opt <- NULL

    for(j in 1:init_vals) {
      if (j == 1) {
        init_val <- .5 * rep(1, T*(T+1)/2)
      } else {
        init_val <- .5 * j * rep(1, T*(T+1)/2)
      }

      opt_res <- optim(par = rep(0, T*(T+1)/2), obj)
      ## control=list(trace=TRUE))
      opt <- min(opt_res$val, opt)

      if (opt == opt_res$val) {
        L <- opt_res$par
        Lambda_opt <- make_from_lowertri(L, T)
      }
      
    for (j in 1:J) {
      thetahat[, j] <- thetahat_j(mu = 0, Lambda_opt, M[[j]], y[, j])
    }
    
  } else if (centering == "gen") {
    ## insert URE for general centering case
  } else {
    ##  insert URE for covariate case
  }

  return(list(thetahat = thetahat, Lambda_opt = Lambda_opt, obj_val = opt))
}
