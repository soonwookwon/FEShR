proj_spg <- function(x, T) {

  A <- matrix(x, nrow = T, ncol = T)
  A <- 0.5 * (t(A) + A) #symmetrize A

  r <- eigen(A, symmetric = TRUE)
  lambda <- r$values
  V <- r$vectors

  lambda <- pmax(0, lambda)
  
  A <- V %*% diag(lambda) %*% t(V)
  x[1:(T^2)] <- as.vector(A)

  return(x)
}
