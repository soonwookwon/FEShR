make_oracle_obj <- function(thetas, y, M, centering, diag) {
  
  T <- nrow(y)
  J <- ncol(y)
  
  if (centering == "0") {
    if (!diag) {
      oracle_obj <- function(L) {
        Lambda <- make_from_lowertri(L, T)
        oracle_loss <- 0
        for (j in 1:J) {
          thetahatj <- thetahat_j(0, Lambda, M[[j]], y[, j])
          oracle_loss_j <- sum((thetas[, j] - thetahatj)^2)
          oracle_loss <- oracle_loss + oracle_loss_j
        }
        return((1/J) * oracle_loss)
      }
    } else {
      oracle_obj <- function(D) {
        Lambda <- diag(D)^2
        oracle_loss <- 0
        for (j in 1:J) {
          thetahatj <- thetahat_j(0, Lambda, M[[j]], y[, j])
          oracle_loss_j <- sum((thetas[, j] - thetahatj)^2)
          oracle_loss <- oracle_loss + oracle_loss_j
        }
        return((1/J) * oracle_loss)
      }
    }
    
  } else if (centering == "gen") {
    ## insert oracle for general centering case
  } else {
    ##  insert oracle for covariate case
  }
  
  return(oracle_obj)
}


##' @title URE
##' @description Calculates URE(mu, Lambda) as defined in Kwon (2020)
##' 
##' @param y a T-by-J data matrix
##' @param M a length J list with the corresponding covariance matrices 
make_ol_deriv <- function(y, M, thetas) {

  T <- nrow(y)
  J <- ncol(y)
  
  ol_deriv <- function(Lambda_entries) {

    Lambda <- matrix(Lambda_entries, nrow = T)
    ol_deriv <- 0

    for (j in 1:J) {
      theta_j <- thetas[, j]
      y_j <- y[, j]
      M_j <- M[[j]]
      ol_deriv <- ol_deriv + ol_deriv_j(Lambda, y_j, M_j, theta_j)
    }
    
    return((1/J) * as.vector(t(ol_deriv))) 
  }
}

##' @title URE
##' @description Calculates URE(mu, Lambda) as defined in Kwon (2020)
##' 
##' @param y a T-by-J data matrix
##' @param M a length J list with the corresponding covariance matrices 
##' @param thetas 
make_ol_deriv_lt <- function(y, M, thetas) {

  T <- nrow(y)
  J <- ncol(y)
  
  ol_deriv_lt <- function(L) {

    Lambda <- make_from_lowertri(L, T)
    ol_deriv <- 0

    for (j in 1:J) {
      theta_j <- thetas[, j]
      y_j <- y[, j]
      M_j <- M[[j]]
      ol_deriv <- ol_deriv + ol_deriv_j(Lambda, y_j, M_j, theta_j)
    }

    L_mat <- matrix(0, nrow = T, ncol = T)
    L_mat[lower.tri(L_mat, diag = TRUE)] <- L
    ol_deriv_lt <- (ol_deriv + t(ol_deriv)) %*% L_mat 
    ol_deriv_lt <- ol_deriv_lt[lower.tri(ol_deriv_lt, diag = TRUE)]
    return((1/J) * as.vector(ol_deriv_lt)) 
  }
}


##' @title Derivative of URE
##' @description Makes the function that calculate the component of URE(mu,
##'   Lambda) that cooresponds to the jth observation ..
##' @param Lambda 
##' @param y_j 
##' @param M_j 
ol_deriv_j <- function(Lambda, y_j, M_j, theta_j) {
  
  inv_Lam_Mj <- chol2inv(chol(Lambda + M_j))
  Mj_inv_Lam_Mj <- M_j %*% inv_Lam_Mj
  ol_deriv_j <- 2 * t(Mj_inv_Lam_Mj) %*% (y_j - theta_j) %*% t(y_j) %*%
    inv_Lam_Mj - 2 * inv_Lam_Mj %*% tcrossprod(y_j) %*% crossprod(Mj_inv_Lam_Mj)

  return(ol_deriv_j)
}


#' @export
get_theta_ol <- function(thetas, y, M, centering, diag = FALSE, n_init_vals) {

  T <- nrow(y)
  J <- ncol(y)
 
  theta_ol <- matrix(0, nrow = T, ncol = J)
  obj <- make_oracle_obj(thetas, y, M, centering, diag = diag)
  grad <- make_ol_deriv_lt(y, M, thetas)

  if (!diag) {
    init_val_mat <- matrix(rnorm(n_init_vals * T * (T + 1) / 2),
                           nrow = n_init_vals)
  } else {
    init_val_mat <- matrix(rnorm(n_init_vals * T), nrow = n_init_vals)
  }
  

  opt <- NULL
  for (j in 1:nrow(init_val_mat)) {
    init_val <- init_val_mat[j, ]
    opt_res <- optim(par = init_val, obj, gr = grad, method = "BFGS")
 
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
  

  for (j in 1:J) {
    theta_ol[, j] <- thetahat_j(mu = 0, Lambda_ol, M[[j]], y[, j])
  }

  ## print(opt)
  ## print(Lambda)
  
  return(list(theta_ol = theta_ol, Lambda_ol = Lambda_ol, obj_val = opt))
}
