##' ...
##'
##' ..
##' @param thetas 
##' @param mu 
##' @param y 
##' @param M 
##' @param diag 
make_oracle_cpp_obj <- function(thetas, mu, y, M, diag) {
  
  T <- nrow(y)
  J <- ncol(y)
  

  if (!diag) {
    oracle_obj <- function(L) {
      Lambda <- make_from_lowertri(L, T)
      return(oracle_loss_cpp(thetas, mu, Lambda, y, M))
    }
  } else {
    oracle_obj <- function(D) {
      Lambda <- diag(D)^2
      return(oracle_loss_cpp(thetas, mu, Lambda, y, M))
    }
  }
  
  return(oracle_obj)
}


##' @title URE
##' @description Calculates URE(mu, Lambda) as defined in Kwon (2020)
##' 
##' @param thetas 
##' @param mu 
##' @param y 
##' @param M 
make_ol_cpp_deriv_lt <- function(thetas, mu, y, M) {

  T <- nrow(y)
  J <- ncol(y)
  
  ol_deriv_lt <- function(L) {

    Lambda <- make_from_lowertri(L, T)
    ol_deriv <- oracle_cpp_deriv(thetas, mu, Lambda, y, M)
    
    L_mat <- matrix(0, nrow = T, ncol = T)
    L_mat[lower.tri(L_mat, diag = TRUE)] <- L
    ol_deriv_lt <- (ol_deriv + t(ol_deriv)) %*% L_mat 
    ol_deriv_lt <- ol_deriv_lt[lower.tri(ol_deriv_lt, diag = TRUE)]

    return(as.vector(ol_deriv_lt)) 
  }
}


##' ...
##'
##' ..
##' @param thetas 
##' @param y 
##' @param M 
##' @param centering 
##' @param diag 
##' @param n_init_vals 
get_theta_ol <- function(thetas, y, M, centering, diag = FALSE, n_init_vals,
                         optim_control) {

  T <- nrow(y)
  J <- ncol(y)
 
  if (centering == "0") {
    mu_ol = matrix(0, nrow = T, ncol = J)
    obj <- make_oracle_cpp_obj(thetas, mu = matrix(0, T, J), y, M, diag = diag)
    grad <- make_ol_cpp_deriv_lt(thetas, mu = matrix(0, T, J), y, M)
  } else if (centering == "gen") {
    
    if (!diag) {
      obj <- function(L) {
        Lambda <- make_from_lowertri(L, T)
        loss <- oracle_loss_cpp(thetas, opt_mu_Lambda_ol(thetas, Lambda, y, M),
                                Lambda, y, M)
        return((1/J) * loss)
      }
    } else {
      obj <- function(D) {
        Lambda <- diag(D)^2
        loss <- oracle_loss_cpp(thetas, opt_mu_Lambda_ol(thetas, Lambda, y, M),
                                Lambda, y, M)
        return((1/J) * loss)
      }
    }
    grad <- NULL
  }
  
  
  if (!diag) {
    init_val_mat <- matrix(rnorm(n_init_vals * T * (T + 1) / 2),
                           nrow = n_init_vals)
  } else {
    init_val_mat <- matrix(rnorm(n_init_vals * T), nrow = n_init_vals)
  }
  

  opt <- NULL
  for (j in 1:n_init_vals) {
    init_val <- init_val_mat[j, ]
    opt_res <- optim(par = init_val, obj, gr = grad, method = "BFGS",
                     control = optim_control)
 
    opt <- min(opt_res$val, opt)
    
    if (opt == opt_res$val) {
      L <- opt_res$par
      if (!diag) {
        Lambda_ol <- make_from_lowertri(L, T)
      } else {
        Lambda_ol <- diag(L)^2
      }
    }
  }

  if (centering == "gen") {
    mu_ol <- opt_mu_Lambda_ol(thetas, Lambda_ol, y, M)
  }

  theta_ol <- get_thetahat(mu_ol, Lambda_ol, y, M)
   
  return(list(theta_ol = theta_ol, Lambda_ol = Lambda_ol, obj_val = opt))
}
