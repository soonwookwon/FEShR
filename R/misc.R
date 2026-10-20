##' ...
##'
##' ..
##' @param L 
##' @param T 
make_from_lowertri <- function(L, T) {
  L_mat <- matrix(0, nrow = T, ncol = T)
  L_mat[lower.tri(L_mat, diag = TRUE)] <- L
  return(tcrossprod(L_mat))
}

##' ...
##'
##' ..
##' @param Lambda 
extract_lowertri <- function(Lambda) {
   return(as.vector(Lambda[lower.tri(Lambda, diag = TRUE)]))
}

.moment_init <- function(y, M, centering, Z, missing_cells, O, diag_lam,
                         periods = NULL) {
  T <- nrow(y)
  J <- ncol(y)
  obs <- matrix(TRUE, T, J)
  if (missing_cells) {
    obs[] <- FALSE
    for (j in 1:J) {
      obs[O[[j]], j] <- TRUE
    }
  }
  if (centering == "cov") {
    X <- do.call("rbind", Z)
    keep <- as.vector(obs)
    gam <- lm.fit(X[keep, , drop = FALSE], as.vector(y)[keep])$coefficients
    gam[is.na(gam)] <- 0
    mu <- matrix(X %*% gam, T, J)
  } else if (centering == "gen") {
    mu <- matrix(rowSums(y * obs) / pmax(rowSums(obs), 1), T, J)
  } else {
    mu <- 0
  }
  r <- y - mu
  r[!obs] <- 0
  Mbar <- matrix(0, T, T)
  for (j in 1:J) {
    o <- which(obs[, j])
    Mbar[o, o] <- Mbar[o, o] + as.matrix(M[[j]])[o, o]
  }
  N <- pmax(tcrossprod(obs), 1)
  Lambda0 <- (tcrossprod(r) - Mbar) / N
  Mavg <- Mbar / N
  if (!is.null(periods)) {
    Lambda0 <- Lambda0[periods, periods, drop = FALSE]
    Mavg <- Mavg[periods, periods, drop = FALSE]
  }
  floor_val <- 1e-6 * max(abs(Lambda0), diag(Mavg), .Machine$double.eps)
  if (diag_lam) return(sqrt(pmax(diag(Lambda0), floor_val)))
  ev <- eigen((Lambda0 + t(Lambda0)) / 2, symmetric = TRUE)
  vals <- pmax(ev$values, floor_val)
  Lambda0 <- ev$vectors %*% (vals * t(ev$vectors))
  return(extract_lowertri(t(chol(Lambda0))))
}
