mdfa.nonstatfilter <- function(dpoly,rem.vec,filter.sharp)
{
  
  #######################################################
  #
  #	mdfa.nonstatfilter by Tucker McElroy
  #
  # computes R and Q constraint matrices for a non-stationary process
  #	inputs:
  #   dpoly is the differencing polynomial
  #   rem.vec: block column vector of d real matrices, 
  #     corresponding to the coefficients
  #     of Psi^{star} (z), starting with constant coefficient.
  #   filter.sharp is a N x N x q-d array of Psi^{sharp} filter coefficients
  #	outputs:
  #   filter.out is a N x N x q array of Psi filter coefficients
  #
  ##############################################################
  
  d <- length(dpoly)-1
  N <- dim(rem.vec)[2]
  q <- dim(filter.sharp)[3]+d
  
  R.mat <- toeplitz(c(dpoly/dpoly[d+1],rep(0,q-1-d)))
  R.mat[upper.tri(R.mat)] <- 0
  R.mat <- R.mat[,-seq(q-d+1,q)] %x% diag(N)
  rem.trans <- t(matrix(aperm(array(t(rem.vec),c(N,N,d)),c(2,1,3)),nrow=N))
  Q.mat <- rbind(rem.trans, matrix(0,nrow=N*(q-d),ncol=N))
  filter.out <- array(t(R.mat %*% t(matrix(filter.sharp,nrow=N)) + 
                              Q.mat),c(N,N,q))
  return(filter.out)
}  