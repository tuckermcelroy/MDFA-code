vmaq.sim <- function(theta.array,innovar.matrix,T.sim)
{
  
  # vmaq.sim by Tucker McElroy
  #
  # Simulates VMA(q) process
  #
  # Inputs:
  #   theta.array: n x n x q dimensional array of coefficients
  #   innovar.matrix: n x n innovation covariance matrix
  
  
  
  eps.old <- rnorm(2)
x.sim <- NULL
for(t in 1:T)
{
  eps.next <- rnorm(2)
  x.next <- eps.next + theta.matrix %*% eps.old
  eps.old <- eps.next
  x.sim <- cbind(x.sim,x.next)
}
x.sim <- ts(t(x.sim))