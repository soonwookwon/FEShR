##' Calculate the URE shrinkage estimator
##'
##' ...
##' @param y a T-by-J data matrix, possibly correlated
##'   within each column. If a vector is provided, it is assumed that T=1 and
##'   will be coerced to a 1-by-length(y) matrix.
##' @param M length J list of the covaraince matrices.
##' @param centering Centering
##' @param type 
##' @export
fe_shrink <- function(y, M, centering = c("0", "gen", "cov"), W = NULL,
                      type = c("URE", "EBMLE"), tau = .95, n_init_vals = 1,
                      all_init_vals = FALSE, optim_control = list()) {

  # Some sanity checks
  if(!is.matrix(y)) {
    warning("y is not a matrix; converting to a matrix assuming T=1.")
    y <- matrix(y, nrow = 1)
  }

  if(!is.list(M)) {
    stop("M is not a list.")
  }

  if(ncol(y) != length(M)) {
    stop("ncol(y) and length(M) must match.")
  }

  centering <- match.arg(centering)
  type <- match.arg(type)

  print(paste0("method=", type, ", centering=", centering))
  T <- nrow(y)
  J <- ncol(y)

  if (is.null(W)) W <- diag(T) # W is always diagonal for now
  
  if (centering == "0") {
    mu_opt <- matrix(0, T, J)
    obj <- gen_obj_0(y, M, type = type)
    grad <- gen_deriv_0(y, M, type = type)
    all_vals <- NULL
    all_pars <- NULL

    init_val_mat <- matrix(rnorm(n_init_vals * T * (T + 1) / 2),
                           nrow = n_init_vals)
    opt <- NULL
    for (j in 1:n_init_vals) {
      init_val <- init_val_mat[j, ]
      opt_res <- optim(par = init_val, gr = grad, obj, method = "BFGS",
                       control = optim_control)
 
      if (all_init_vals) {
        all_vals <- c(all_vals, opt_res$val)
        all_pars <- rbind(all_pars, opt_res$par)

        if (j == n_init_vals) {
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
    
    thetahat <- get_thetahat(mu_opt, Lambda_opt, y, M) 
    
  } else if (centering == "gen") {

    mu_abs_bounds <- sapply(1:T, function(t) quantile(y[t, ], tau))
    obj <- gen_obj_gen(y, M, type = type, mu_abs_bounds)
    
    all_vals <- NULL
    all_pars <- NULL
    
    init_val_mat <- matrix(rnorm(n_init_vals * T * (T + 1) / 2),
                           nrow = n_init_vals)
    opt <- NULL
    
    for (j in 1:n_init_vals) {
      init_val <- init_val_mat[j, ]
      opt_res <- optim(par = init_val, obj, method = "BFGS",
                       control = optim_control)
      
      if (all_init_vals) {
        all_vals <- c(all_vals, opt_res$val)
        all_pars <- rbind(all_pars, opt_res$par)

        if (j == n_init_vals) {
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

    if (type == "URE") {
      mu_opt <- opt_mu_Lambda_URE(Lambda_opt, y, M, mu_abs_bounds)
    } else if (type == "EBMLE") {
      mu_opt <- opt_mu_Lambda_nll(Lambda_opt, y, M)
    }
        
    thetahat <- get_thetahat(mu_opt, Lambda_opt, y, M)
    
  } else {
    ##  insert URE for covariate case
  }
  
  return(list(thetahat = thetahat, mu = mu_opt, Lambda_opt = Lambda_opt,
              obj_val = opt, all_vals = all_vals, all_pars = all_pars))
}
