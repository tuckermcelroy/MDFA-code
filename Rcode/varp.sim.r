varp.sim <- function(phi.array,innovar.matrix,T.sim)
{
  
  # varp.sim by Tucker McElroy
  #
  # Simulates VAR(p) process
  #
  # Inputs:
  #   phi.array: n x n x p dimensional array of coefficients
  #   innovar.matrix: n x n innovation covariance matrix
  
  p <- dim(phi.array)[3]
  gamma <- VARMAauto(phi.array,NULL,innovar.matrix,10)
  gamma.0 <- gamma[,,1]
  x.init <- t(chol(gamma.0)) %*% rnorm(2)
  x.next <- x.init
  x.sim <- NULL
  for(t in 1:T.sim)
  {
    x.next <- phi.matrix %*% x.next + rnorm(2)
    x.sim <- cbind(x.sim,x.next)
  }
  x.sim <- ts(t(x.sim))
}