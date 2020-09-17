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
                      type = c("URE", "EBMLE"), init_vals = 1,
                      lam_range = c(.1, 1), all_init_val = FALSE,
                      use_EBMLE_opt = FALSE,
                      use_DEoptim = FALSE,
                      print_by_init_val = FALSE) {
  
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
        obj <- make_URE_obj(y, M, centering)
        grad <- make_URE_deriv_lt(y, M)
        opt_alg <- "BFGS"
      } else {
        obj <- make_nll_obj(y, M, centering)
        grad <- NULL
        opt_alg <- "BFGS"
      }

    all_vals <- NULL
    all_pars <- NULL
    
    if (!use_DEoptim) {
      
      opt <- NULL
      init_val_mat <- matrix(0, nrow = init_vals, ncol = T*(T+1)/2)

      ## for(j in 1:init_vals) {
      ##   if (j == 1) {
      ##     init_val_mat[j, ] <- .5 * rep(1, T*(T+1)/2)
      ##   } else {
      ##     init_val_mat[j, ] <-.1 * j + j * rnorm(T*(T+1)/2)
      ##   }
      ## }

      init_diags <- t(rep(lam_range[1], T) + (lam_range[2] - lam_range[1]) *
                        t(randtoolbox::sobol(init_vals, dim = T)))
      init_val_mat <- sapply(1:init_vals,
                             function(j) extract_lowertri(diag(init_diags[j, ])))
      init_val_mat <- t(init_val_mat)
      
      lam_seq <- seq(lam_range[1], lam_range[2], length.out = init_vals)
      init_vals_equl <- t(sapply(1:init_vals,
                                 function(j) {
                                   extract_lowertri(lam_seq[j] * diag(T))
                                 })
                          )
      init_val_mat <- rbind(init_vals_equl, init_val_mat)

      # EBMLE opt par
      if (use_EBMLE_opt) {
        EBMLE_obj <- make_nll_obj(y, M, centering)
        EBMLE_opt_par <-
          optim(par = init_val_mat[1, ],
                obj,
                method = "BFGS")$par
        init_val_mat <- rbind(EBMLE_opt_par, init_val_mat)
      }
      ## solver <- ifelse(use_optimx, optimx::optimx, optim)
      ## opt_res <- solver(par = init_val, obj)
      
      for (j in 1:nrow(init_val_mat)) {
        init_val <- init_val_mat[j, ]
        opt_res <- optim(par = init_val, gr = grad, obj, method = opt_alg)
        
        ## control=list(trace=TRUE))
        ## if (verbose) {
        ##   print(paste0(j, "th initial value: ",
        ##                paste0(init_val, collapse = ", ")))
        ##   print("Min")
        ##   print(opt_res$val)
        ##   print("Par")
        ##   print(opt_res$par)
        ## }

        if (all_init_val) {
          all_vals <- c(all_vals, opt_res$val)
          all_pars <- rbind(all_pars, opt_res$par)

          if (print_by_init_val) {
            if (j %% 5 == 0) {
              print(paste0(rep("=", 80), collapse = ""))
              print(all_vals)
            }
          }
          
          ## print(all_pars)
        }
        
        opt <- min(opt_res$val, opt)

        if (opt == opt_res$val) {
          L <- opt_res$par
          Lambda_opt <- make_from_lowertri(L, T)
        }
      }
    } else if (use_DEoptim) {
      opt_res <- DEoptim::DEoptim(obj,
                                  lower = rep(-2, T * (T+1) / 2 ),
                                  upper = rep(2, T * (T+1) / 2 ),
                                  control = list(itermax = 200))
      opt <- opt_res$optim$bestval
      ## print(opt_res)  
      L <- opt_res$optim$bestmem
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
  print(opt)
  print(Lambda_opt)
  return(list(thetahat = thetahat, Lambda_opt = Lambda_opt, obj_val = opt,
              all_vals = all_vals, all_pars = all_pars))
}
