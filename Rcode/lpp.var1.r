lpp.var1 <- function(theta,psi.array,acf.array)
{
  ###  lpp.var1
  #		theta: pre-parameters, N^2 vector of real numbers ;
  #			uses var.pre2par to generate phi matrix
  #		psi.array: the target filter at lags -1 through -L,
  #			as an 1 x N x L array object
  #		acf.array: given acf at lags 0 through L,
  #			as an N x N x L array object
  
  N <- dim(psi.array)[2]
  L <- dim(psi.array)[3]
  phi.matrix <- var.pre2par(theta,1,N)
  phi.matrix <- phi.matrix[,,1]
  phi.next <- phi.matrix
  A.psi.phi <- matrix(0,1,N)
  sum.single <- matrix(0,1,N)
  sum.double <- 0
  for(h in 1:L)
  {
    A.psi.phi <- A.psi.phi + psi.array[,,h] %*% phi.next
    phi.next <- phi.matrix %*% phi.next
    sum.single <- sum.single + psi.array[,,h] %*% acf.array[,,(h+1)]
    for(k in 1:L)
    {
      if(h >= k) { covar <- acf.array[,,(h-k+1)] } else { covar <- t(acf.array[,,(k-h+1)]) }
      sum.double <- sum.double + t(as.matrix(psi.array[,,h])) %*% covar %*% as.matrix(psi.array[,,k])
    }
  }
  lpp.criterion <- sum.double - A.psi.phi %*% t(sum.single) - 
    sum.single %*% t(A.psi.phi) + A.psi.phi %*% acf.array[,,1] %*% t(A.psi.phi)
  return(lpp.criterion)
}
