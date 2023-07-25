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
  n <- dim(phi.array)[1]
  phi.matrix <- matrix(phi.array,nrow=n)
  comp.matrix <- diag(n*(p+1))[1:(n*p),(n+1):(n*p+n)]
  comp.matrix[1:n,] <- phi.matrix
  
  gamma <- VARMAauto(phi.array,NULL,innovar.matrix,p+1)
  gamma.init <- array(0,c(n,p,n,p))
  for(j in 1:p)
  {
    for(k in 1:p)
    {
      if(j >= k) { gamma.init[,j,,k] <- gamma[,,j-k+1] } else 
        { gamma.init[,j,,k] <- t(gamma[,,k-j+1]) }
    }
  }
  gamma.0 <- matrix(gamma.init,nrow=n*p)
  x.init <- t(chol(gamma.0)) %*% rnorm(n*p)
  x.next <- x.init
  x.sim <- NULL
  for(t in 1:T.sim)
  {
    x.next <- comp.matrix %*% x.next + matrix(c(1,rep(0,p-1)) %x% rnorm(n),ncol=1)
    x.sim <- cbind(x.sim,x.next[1:n])
  }
  x.sim <- ts(t(x.sim))
  return(x.sim)
}