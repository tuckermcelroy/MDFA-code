mdfa.getremainder <- function(frf,rootfreqs)
{
  
  #######################################################
  #
  #	mdfa.getremainder by Tucker McElroy
  #
  #	computes remainder polynomial after having divided by Delta(z).
  #
  # Background: consider Psi(z) = Delta(z)Psi^{sharp}(z) + Psi^{star}(z),
  #   where Delta(z) has degree d and Psi^{star}(z) has degree d-1.
  #   Suppose there are r roots zeta_j with multiplicity m_j.
  #   Then Psi^(k) (zeta_j) = {Psi^{star}}^k (zeta_j) for 0 <= k < m_j.
  #	inputs:
  #		frf is array N x N x Grid of complex entries, the target
  #			frequency response function Psi(e^{-i lambda})
  #			for lambda given by Grid number of Fourier frequencies 
  #		rootfreqs: lists frequencies x in [-1,1] such that lambda=pi*x, 
  #			with repeats for double roots, 
  #			such that exp(-i*lambda) is a root of Delta(z).
  #     Presumes both x and -x are included if root is non-real.
  #	outputs:
  #   rem.vec: block column vector of d real matrices, 
  #     corresponding to the coefficients
  #     of Psi^{star} (z), starting with constant coefficient.
  #
  #######################################################
  
  N <- dim(frf)[1]
  grid <- dim(frf)[3]
  m <- floor(grid/2)
  d <- length(rootfreqs)
  
  sig.lambdas <- unique(rootfreqs)
  sig.mults <- NULL
  if(length(rootfreqs)>0) 
  {
    for(j in 1:length(sig.lambdas))
    {
      sig.mults <- c(sig.mults,sum(rootfreqs == sig.lambdas[j]))
    } 
  }
  
  sig.vec <- NULL
  if(length(sig.lambdas) > 0) {
    for(k in 1:length(sig.lambdas))
    {
      j.star <- floor(sig.lambdas[k]*grid/2) + m+1
      if(j.star==(grid+1)) j.star <- 1
      sig.vec.new <- frf[,,j.star]  
      sig.vec <- rbind(sig.vec,sig.vec.new)
      if(sig.mults[k]==2)
      {
        if(j.star==grid) 
        { 
          sig.vec.new <- 1i*exp(1i*sig.lambdas[k]*pi)*
            (frf[,,1]-frf[,,j.star])*grid/(2*pi) 
        } else 
        { 
          sig.vec.new <- 1i*exp(1i*sig.lambdas[k]*pi)*
            (frf[,,j.star+1]-frf[,,j.star])*grid/(2*pi) 
        }
        sig.vec <- rbind(sig.vec,sig.vec.new)
      }
    } 
  }
  
  W.mat <- NULL
  for(j in 1:length(sig.lambdas))
  {
    zeta <- exp(1i*pi*sig.lambdas[j])
    for(k in 1:sig.mults[j])
    {
      weights <- c(rep(0,k-1),factorial(k-1)*choose(seq(k-1,d-1),k-1))
      W.mat <- rbind(W.mat,weights*zeta^(seq(1,d)-k))
    }
  }
  
  rem.vec <- solve(W.mat %x% diag(N),sig.vec)
  rem.vec <- Re(rem.vec)
  
  return(rem.vec)
  
}