##' ...
##'
##' ..
##' @param thetas 
##' @param y 
##' @param M 
##' @param diag_lam 
##' @param centering 
make_oracle_obj <- function(thetas, y, M, W, missing_cells, O, diag_lam,
                            centering, Z = list()) {
  
  T <- nrow(y)
  J <- ncol(y)

  if (centering == "0") {
    mu_opt <- matrix(0, nrow = T, ncol = J)
    if (!diag_lam) {
      obj <- function(L) {
        Lambda <- make_from_lowertri(L, T)
        return(oracle_loss(thetas, mu_opt, Lambda, y, M, W, missing_cells, O))
      }
    } else {
      obj <- function(D) {
        Lambda <- diag(D)^2
        return(oracle_loss(thetas, mu_opt, Lambda, y, M, W, missing_cells, O))
      }
    }
  } else if (centering == "gen" | centering == "cov") {

    if (!diag_lam) {
      obj <- function(L) {
        Lambda <- make_from_lowertri(L, T)
        loss <- oracle_loss(thetas,
                            opt_mu_Lambda_ol(thetas, Lambda, y, M, W,
                                             missing_cells, O, Z),
                            Lambda, y, M, W, missing_cells, O)
        return((1/J) * loss)
      }
    } else {
      obj <- function(D) {
        Lambda <- diag(D)^2
        loss <- oracle_loss(thetas,
                            opt_mu_Lambda_ol(thetas, Lambda, y, M, W,
                                             missing_cells, O, Z),
                            Lambda, y, M, W, missing_cells, O)
        return((1/J) * loss)
      }
    }
  }

  
  return(obj)
}


##' @title URE
##' @description Calculates URE(mu, Lambda) as defined in Kwon (2020)
##' 
##' @param thetas 
##' @param y 
##' @param M 
##' @param W 
##' @param missing_cells 
##' @param O 
##' @param diag_lam 
make_ol_deriv_lt <- function(thetas, y, M, W, missing_cells, O, diag_lam) {

  T <- nrow(y)
  J <- ncol(y)

  mu_opt <- matrix(0, nrow = T, ncol = J)
  
  if (!diag_lam) {
    ol_deriv_lt <- function(L) {
      
      Lambda <- make_from_lowertri(L, T)
      ol_deriv <- oracle_deriv(thetas, mu_opt, Lambda, y, M, W, missing_cells, O)
      
      L_mat <- matrix(0, nrow = T, ncol = T)
      L_mat[lower.tri(L_mat, diag = TRUE)] <- L
      ol_deriv_lt <- (ol_deriv + t(ol_deriv)) %*% L_mat 
      ol_deriv_lt <- ol_deriv_lt[lower.tri(ol_deriv_lt, diag = TRUE)]

      return(as.vector(ol_deriv_lt)) 
    }
  } else {
    ol_deriv_lt <- function(D) {

      Lambda <- diag(D)^2
      ol_deriv <- oracle_deriv(thetas, mu_opt, Lambda, y, M, W, missing_cells, O)
      
      return(2 * diag(ol_deriv) * D) 
    }
  }

  return(ol_deriv_lt)
}


##' ...
##'
##' ..
##' @param thetas 
##' @param y 
##' @param M 
##' @param centering 
##' @param diag_lam 
##' @param W 
##' @param n_init_vals 
##' @param optim_control 
get_theta_ol <- function(thetas, y, M, centering, diag_lam = FALSE, W = NULL,
                         Z = list(), n_init_vals, all_init_vals = FALSE,
                         optim_control) {

  print(paste0("Caculating the oracle for centering=", centering))
  if (centering != "cov") Z <- list()
  T <- nrow(y)
  J <- ncol(y)
  
  if (is.null(W)) W <- diag(T)

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

  obj <- make_oracle_obj(thetas, y, M, W, missing_cells, O, diag_lam, centering,
                         Z = Z)
 
  if (centering == "0") {
    grad <- make_ol_deriv_lt(thetas, y, M, W, missing_cells, O, diag_lam)
  } else {
    grad <- NULL
  }
  

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

    opt_res <- optim(par = init_val, obj, gr = grad, method = "BFGS",
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
      if (!diag_lam) {
        Lambda_ol <- make_from_lowertri(L, T)
      } else {
        Lambda_ol <- diag(L)^2
      }
    }
  }

  if (centering == "0") {
    mu_ol <- matrix(0, nrow = T, ncol = J)
  } else  {
    mu_ol <- opt_mu_Lambda_ol(thetas, Lambda_ol, y, M, W, missing_cells, O, Z)
  }

  theta_ol <- get_thetahat(mu_ol, Lambda_ol, y, M, missing_cells, O)
   
  return(list(theta_ol = theta_ol, Lambda_ol = Lambda_ol, obj_val = opt))
}
