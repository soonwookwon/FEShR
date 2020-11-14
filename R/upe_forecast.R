##' Calculate the UPE forecast
##'
##' This function calculates the UPE forecases as described in Kwon (2020)
##' 
##' @param y T-by-J data matrix, possibly correlated within each column. If a
##'   vector is provided, it is assumed that T=1 and will be coerced to a
##'   1-by-length(y) matrix.
##' @param M Length J list of the covaraince matrices.
##' @param centering Location to shrink toward. Currently, only shrinkage toward
##'   0 is supported
##' @param n_init_vals Number of initial values to try for optimization
##' @param all_init_vals If TRUE, print and return optimization results for all
##'   initial values tried
##' @param optim_control List of control variables to be passed to \code{optim}
##' @param K Constant used to define the parameter space \eqn{\mathcal{L}}
##' @return Returns a list with the following components \describe{
##'   \item{thetahat}{1-by-J matrix of the optimal forecasts}
##'   \item{LambTmT_opt}{T-1 vector of the optimal \eqn{\Lambda_{T,-T}} }
##'   \item{Lambda_opt}{(T-1)-by-T matrix of the optimal \eqn{\Lambda_{-T}} }
##'   \item{obj_val}{Minimized value of the objective function for the
##'   optimization problem solved to tune the hyperparameters}
##'   \item{all_vals}{Minimized value of the objective function for each initial
##'   value. Returns \code{NULL} if \code{all_init_vals==FALSE}}
##'   \item{all_pars}{Optimal hyperparameters for each initial value. Returns
##'   \code{NULL} if \code{all_init_vals==FALSE}}
##'   \item{all_thetahats}{Optimal
##'   forecasts for each initial value. Returns \code{NULL} if
##'   \code{all_init_vals==FALSE}} }
##'
##' @references{ \cite{Soonwoo Kwon. 2020.
##'   "Optimal Shrinkage Estimation of Fixed Effects in Linear Panel Data Models."
##'   \emph{working paper}.}  }
##' @export
upe_forecast <- function(y, M, centering = c("0"),
                         n_init_vals = 1, all_init_vals = FALSE,
                         optim_control = list(), K = 1e03) {

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
  
  print(paste0("centering=", centering))
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
  ## Fix this part later
  bounds <- K * 1
  if (centering == "0") {
    obj <-  function(L) { 
      Lambda <- make_from_lowertri(L, T - 1)
      return(UPE(Lambda, y, M, missing_cells, O, bounds))
    }
  }
  
  all_vals <- NULL
  all_pars <- NULL
  all_thetahats <- NULL
  opt <- NULL
  init_val_mat <- matrix(rnorm(n_init_vals * (T - 1) * T / 2),
                         nrow = n_init_vals)
  
  for (j in 1:n_init_vals) {
    init_val <- init_val_mat[j, ]
    opt_res <- optim(par = init_val, obj, method = "BFGS",
                     control = optim_control)
    
    if (all_init_vals) {
      L <- opt_res$par
      Lambda_opt <- make_from_lowertri(L, T - 1)
      LamTmT_opt <- opt_LamTmT_Lambda(Lambda_opt, y, M, missing_cells, O)
      thetahat <- get_thetahat_fc(LamTmT_opt, Lambda_opt, y, M, missing_cells, O)
      all_thetahats <- rbind(all_thetahats, thetahat)
      all_vals <- c(all_vals, opt_res$val)
      all_pars <- rbind(all_pars, opt_res$par)

      ## if (j == n_init_vals) {
        print(paste0(rep("=", 80), collapse = ""))
        print(all_vals)
        print(all_pars)
      ## }
    }
    
    opt <- min(opt_res$val, opt)

    if (opt == opt_res$val) {
      L <- opt_res$par
      Lambda_opt <- make_from_lowertri(L, T - 1)
    }
  }
  
  LamTmT_opt <- opt_LamTmT_Lambda(Lambda_opt, y, M, missing_cells, O)
  thetahat <- get_thetahat_fc(LamTmT_opt, Lambda_opt, y, M, missing_cells, O)
  
  return(list(thetahat = thetahat, LamTmT_opt = LamTmT_opt,
              Lambda_opt = Lambda_opt,
              obj_val = opt, all_vals = all_vals,
              all_pars = all_pars,
              all_thetahats = all_thetahats))
}
