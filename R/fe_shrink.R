##' Calculate URE estimator for the fixed effects
##'
##' This function calculates the URE estimator for the fixed effects, as
##' described in Kwon (2020). If \code{type="EBMLE"}, the Empirical Bayes MLE is
##' calculated instead.
##' 
##' @param y a T-by-J data matrix that consists of the least sqaures estimators
##'   of the fixed effects. If a vector is provided, it is assumed that T=1 and
##'   will be coerced to a 1-by-length(y) matrix.
##' @param M length J list of the varaince matrices
##' @param centering location to shrink toward. If \code{centering="0"} the
##'   estimator that shrinks to zero is computed (this is the same with
##'   shrinking to the grand mean if the fixed effects are demeaned for each
##'   period); \code{centering="gen"} shrinks the data to a general data-driven
##'   location; and \code{centering="cov"} shrinks to a linear combination of
##'   the covariates.
##' @param W T-by-T weight matrix.
##' @param Z length J list of the T-by-k covariate matrices
##' @param type type of shrinkage estimator; possible options are \code{URE} and
##'   \code{EBMLE}
##' @param tau \eqn{\tau} used to define the set \eqn{\mathcal{M}}
##' @param n_init_vals Number of initial values to try for optimization
##' @param all_init_vals If TRUE, print and return optimization results for all
##'   initial values tried
##' @param optim_control list of control variables to be passed to \code{optim}
##' @return Returns a list with the following components \describe{
##'   \item{thetahat}{T-by-J matrix of the optimal shrinkage estimator}
##'   \item{mu_opt}{ T-by-J matrix of optimal shrinkage location. Unless
##'   \code{centering="cov"}, all columns are the same.}
##'   \item{Lambda_opt}{T-by-T matrix of the optimal \eqn{\Lambda} }
##'   \item{obj_val}{Minimized value of the objective function for the
##'   optimization problem solved to tune the hyperparameters}
##'   \item{all_vals}{Minimized value of the objective function for each initial
##'   value. Returns \code{NULL} if \code{all_init_vals==FALSE}}
##'   \item{all_pars}{Optimal hyperparameters for each initial value. Returns
##'   \code{NULL} if \code{all_init_vals==FALSE}} }
##'
##' @references{
##' \cite{Soonwoo Kwon. 2020.
##' "Optimal Shrinkage Estimation of Fixed Effects in Linear Panel Data Models." \emph{working paper}.}
##' }
###' @export
fe_shrink <- function(y, M, centering = c("0", "gen", "cov"), W = NULL,
                      Z = list(), type = c("URE", "EBMLE"), diag_lam = FALSE,
                      tau = .95, B = 1e03, n_init_vals = 1,
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

  # Save missing indices and replace with 0

  missing_cells <- any(is.na(y))
  
  if (missing_cells) {
    O <- vector("list", length = J)
    for (j in 1:J) {
      O[[j]] <- which(!is.na(y[, j]))
      y[-O[[j]], j] <- 0
    }
  } else {
    O <- NULL
  }
  
  if (is.null(W)) W <- diag(T) # W is always diagonal for now
  

  
  if (centering == "0") {
    grad <- gen_deriv(y, M, type = type, W, missing_cells, O, diag_lam)
    bounds <- NULL
    Z <- list()
  } else if (centering == "gen") {
    bounds <- sapply(1:T, function(t) quantile(y[t, ], tau))
    grad <- NULL
    Z <- list()
  } else if (centering == "cov") {
    bounds <-  B * abs(lm(as.vector(y) ~ do.call("rbind", Z))$coef)
  }
  
  obj <- gen_obj(y, M, type, W, missing_cells, O, diag_lam, centering,
                 bounds, Z)
  
  
  if (!diag_lam) {
    init_val_mat <- matrix(rnorm(n_init_vals * T * (T + 1) / 2),
                           nrow = n_init_vals)
  } else {
    init_val_mat <- matrix(rnorm(n_init_vals * T), nrow = n_init_vals)
  }
  
  opt <- NULL
  all_vals <- NULL
  all_pars <- NULL
  
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
      if (!diag_lam) {
        L <- opt_res$par
        Lambda_opt <- make_from_lowertri(L, T)
      } else {
        D <- opt_res$par
        Lambda_opt <- diag(D)^2
      }
    }
  }
  
  if (centering == "0") {
    mu_opt <- matrix(0, T, J)
  } else {  
    if (type == "URE") {
      mu_opt <- opt_mu_Lambda_URE(Lambda_opt, y, M, bounds, W,
                                  missing_cells, O, Z)
    } else if (type == "EBMLE") {
      mu_opt <- opt_mu_Lambda_nll(Lambda_opt, y, M, missing_cells, O, Z)
    }

    ## if (centering == "cov") {
    ##   mu_opt <- gam_to_mu(mu_opt, Z, T)
    ## }
  } 
  
  thetahat <- get_thetahat(mu_opt, Lambda_opt, y, M, missing_cells, O)
  
  
  
  return(list(thetahat = thetahat, mu_opt = mu_opt, Lambda_opt = Lambda_opt,
              obj_val = opt, all_vals = all_vals, all_pars = all_pars))
}
