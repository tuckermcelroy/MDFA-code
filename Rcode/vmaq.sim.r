vmaq.sim <- function(theta.array,innovar.matrix,T.sim)
{
  
  # vmaq.sim by Tucker McElroy
  #
  # Simulates VMA(q) process
  #
  # Inputs:
  #   theta.array: n x n x q dimensional array of coefficients
  #   innovar.matrix: n x n innovation covariance matrix
  
  q <- dim(theta.array)[3]
  n <- dim(theta.array)[1]
  theta.matrix <- matrix(theta.array,nrow=n)
  
  eps.init <- NULL
  for(j in 1:q)
  {
    eps.init <- rbind(eps.init,t(chol(innovar.matrix)) %*% rnorm(n))
  }
  eps.next <- eps.init
  x.sim <- NULL
  for(t in 1:T.sim)
  {
    x.next <- rnorm(n) + theta.matrix %*% eps.next
    eps.next <- rbind(t(chol(innovar.matrix)) %*% rnorm(n),eps.next[1:(n*q-n)])
    x.sim <- cbind(x.sim,x.next)
  }
  x.sim <- ts(t(x.sim))
  return(x.sim)
}