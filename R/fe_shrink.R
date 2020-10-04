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
                      type = c("URE", "EBMLE"), n_init_vals = 1,
                      all_init_val = FALSE, print_by_init_val = FALSE,
                      verbose = FALSE) {
  
  if(!is.matrix(y)) {
    warning("y is not a matrix, assuming T=1.")
    y <- matrix(y, nrow = 1)
  }

  centering <- match.arg(centering)
  type <- match.arg(type)

  print(paste0("method=", type, "#init_vals=", n_init_vals))
  T <- nrow(y)
  J <- ncol(y)

  if (is.null(W)) W <- diag(T) # W is always diagonal for now
  
  thetahat <- matrix(0, nrow = T, ncol = J)
  
  if (centering == "0") {

    if (type == "URE") {
        obj <- make_URE_obj(y, M, centering)
        grad <- make_URE_deriv_lt(y, M)
      } else {
        obj <- make_nll_obj(y, M, centering)
        grad <- make_nll_deriv_lt(y, M)
      }

    all_vals <- NULL
    all_pars <- NULL
    
    
    # EBMLE opt par
    ## if (use_EBMLE_opt) {
    ##   EBMLE_obj <- make_nll_obj(y, M, centering)
    ##   EBMLE_opt_par <-
    ##     optim(par = init_val_mat[1, ],
    ##           obj,
    ##           method = "BFGS")$par
    ##   init_val_mat <- rbind(EBMLE_opt_par, init_val_mat)
    ## }
    ## solver <- ifelse(use_optimx, optimx::optimx, optim)
    ## opt_res <- solver(par = init_val, obj)

    init_val_mat <- matrix(rnorm(n_init_vals * T * (T + 1) / 2),
                           nrow = n_init_vals)
    opt <- NULL
    
    for (j in 1:nrow(init_val_mat)) {
      init_val <- init_val_mat[j, ]
      opt_res <- optim(par = init_val, gr = grad, obj, method = "BFGS")
      
      if (all_init_val) {
        all_vals <- c(all_vals, opt_res$val)
        all_pars <- rbind(all_pars, opt_res$par)

        if (print_by_init_val) {         
          print(paste0(rep("=", 80), collapse = ""))
          print(all_vals)
          print(all_pars)
        }       
      }
      
      opt <- min(opt_res$val, opt)

      if (opt == opt_res$val) {
        L <- opt_res$par
        Lambda_opt <- make_from_lowertri(L, T)
      }
    }
    
    
    for (j in 1:J) {
      thetahat[, j] <- thetahat_j(mu = 0, Lambda_opt, M[[j]], y[, j])
    }
    
  } else if (centering == "gen") {
    ## insert URE for general centering case
  } else {
    ##  insert URE for covariate case
  }
  if (verbose) {
    print(opt)
    print(Lambda_opt)
  }
  
  return(list(thetahat = thetahat, Lambda_opt = Lambda_opt, obj_val = opt,
              all_vals = all_vals, all_pars = all_pars))
}
