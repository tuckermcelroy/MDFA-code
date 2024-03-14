mdfa.getquotient2 <- function(frf,rootfreqs,rem.vec)
{
  
  #######################################################
  #
  #	mdfa.getquotient by Tucker McElroy
  #
  #	computes quotient polynomial after having divided by Delta(z).
  #
  # Background: consider Psi(z) = Delta(z)Psi^{sharp}(z) + Psi^{star}(z),
  #   where Delta(z) has degree d and Psi^{star}(z) has degree d-1.
  #   Suppose there are r roots zeta_j with multiplicity m_j.
  #   Then Psi^(k) (zeta_j) = {Psi^{star}}^k (zeta_j) for 0 <= k < m_j.
  #   We compute Psi^{sharp}(z) at the roots of Delta(z)
  #	inputs:
  #		frf is array N x N x Grid of complex entries, the target
  #			frequency response function Psi(e^{-i lambda})
  #			for lambda given by Grid number of Fourier frequencies 
  #		rootfreqs: lists frequencies x in [-1,1] such that lambda=pi*x, 
  #			with repeats for double roots, 
  #			such that exp(-i*lambda) is a root of Delta(z).
  #     Presumes both x and -x are included if root is non-real.
  #   rem.vec: block column vector of d real matrices, 
  #     corresponding to the coefficients
  #     of Psi^{star} (z), starting with constant coefficient.
  #	outputs:
  #		frf.psi.sharp is array N x N x Grid of complex entries, the  
  #			frequency response function Psi^{sharp}(e^{-i lambda})
  #			for lambda given by Grid number of Fourier frequencies 
  #
  #######################################################
  
  ceps2wold <- function(ceps,q)
  {
    m <- length(ceps)
    if(q > m) { ceps <- c(ceps,rep(0,q-m)) }
    wold <- 1
    wolds <- wold
    for(j in 1:q)
    {
      wold <- sum(seq(1,j)*ceps[1:j]*wolds[j:1])/j
      wolds <- c(wolds,wold)
    }
    return(wolds)
  }
  
  roots2ceps <- function(roots,m)
  {
    p <- length(roots)
    ceps <- rep(0,m)
    for(k in 1:m)
    {	
      ceps[k] <- -1*sum(roots^(-k))/k
    }
    return(ceps)
  }
  
  N <- dim(frf)[1]
  grid <- dim(frf)[3]
  m <- floor(grid/2)
  d <- length(rootfreqs)

  root.ceps <- Re(roots2ceps(exp(-1i*pi*rootfreqs),length(rootfreqs)))
  delta <- ceps2wold(root.ceps,length(rootfreqs))
  delta <- delta/delta[d+1]
  frf.star <- mdfa.frf(aperm(array(rem.vec,c(N,d,N)),c(1,3,2)),0,grid)
  frf.delta <- mdfa.frf(array(delta,c(1,1,d+1)),0,grid)
  
  sig.lambdas <- unique(rootfreqs)
  sig.mults <- NULL
  if(length(rootfreqs)>0) 
  {
    for(j in 1:length(sig.lambdas))
    {
      sig.mults <- c(sig.mults,sum(rootfreqs == sig.lambdas[j]))
    } 
  }
  
  j.stars <- NULL
  sig.vec <- NULL
  sigstar.vec <- NULL
  sigdelta.vec <- NULL
  if(length(sig.lambdas) > 0) {
    for(k in 1:length(sig.lambdas))
    {
      j.star <- floor(sig.lambdas[k]*grid/2) + m+1
      if(j.star==(grid+1)) j.star <- 1
      j.stars <- c(j.stars,j.star) 
      if(sig.mults[k]==1)
      {
        if(j.star==grid) 
        { 
          sig.vec.new <- (frf[,,1]-frf[,,j.star])*grid/(2*pi) 
          sigstar.vec.new <- (frf.star[,,1]-frf.star[,,j.star])*grid/(2*pi) 
          sigdelta.vec.new <- (frf.delta[,,1]-frf.delta[,,j.star])*grid/(2*pi)
        } else 
        { 
          sig.vec.new <- (frf[,,j.star+1]-frf[,,j.star])*grid/(2*pi) 
          sigstar.vec.new <- (frf.star[,,j.star+1]-frf.star[,,j.star])*grid/(2*pi) 
          sigdelta.vec.new <- (frf.delta[,,j.star+1]-frf.delta[,,j.star])*grid/(2*pi) 
        }
        sig.vec <- rbind(sig.vec,sig.vec.new)        
        sigstar.vec <- rbind(sigstar.vec,sigstar.vec.new)
        sigdelta.vec <- rbind(sigdelta.vec,sigdelta.vec.new)
      }
      if(sig.mults[k]==2)
      {
        if(j.star==grid) 
        { 
          sig.vec.new <- (frf[,,1]-2*frf[,,j.star]+frf[,,j.star-1])*(grid/(2*pi))^2
          sigstar.vec.new <- (frf.star[,,1]-2*frf.star[,,j.star]+frf.star[,,j.star-1])*(grid/(2*pi))^2 
          sigdelta.vec.new <- (frf.delta[,,1]-2*frf.delta[,,j.star]+frf.delta[,,j.star-1])*(grid/(2*pi))^2 
        } else 
        if(j.star==1)
        {
          sig.vec.new <- (frf[,,j.star+1]-2*frf[,,j.star]+frf[,,grid])*(grid/(2*pi))^2 
          sigstar.vec.new <- (frf.star[,,j.star+1]-2*frf.star[,,j.star]+frf.star[,,grid])*(grid/(2*pi))^2 
          sigdelta.vec.new <- (frf.delta[,,j.star+1]-2*frf.delta[,,j.star]+frf.delta[,,grid])*(grid/(2*pi))^2 
        } else
        { 
          sig.vec.new <- (frf[,,j.star+1]-2*frf[,,j.star]+frf[,,j.star-1])*(grid/(2*pi))^2 
          sigstar.vec.new <- (frf.star[,,j.star+1]-2*frf.star[,,j.star]+frf.star[,,j.star-1])*(grid/(2*pi))^2 
          sigdelta.vec.new <- (frf.delta[,,j.star+1]-2*frf.delta[,,j.star]+frf.delta[,,j.star-1])*(grid/(2*pi))^2 
        }
        sig.vec <- rbind(sig.vec,sig.vec.new)        
        sigstar.vec <- rbind(sigstar.vec,sigstar.vec.new)
        sigdelta.vec <- rbind(sigdelta.vec,sigdelta.vec.new)
      }
    } 
  }
  
  frf.psi.sharp <- frf.psi - frf.star
  for(j in 1:grid)
  {
    if(!is.element(j,j.stars)) 
    {
      frf.psi.sharp[,,j] <- frf.psi.sharp[,,j]/frf.delta[1,1,j]
    } else
    {
      index <- which(j.stars==j)
      frf.psi.sharp[,,j] <- (sig.vec[index] - sigstar.vec[index])/sigdelta.vec[index]
    }
  }
  
  return(frf.psi.sharp)
  
}  