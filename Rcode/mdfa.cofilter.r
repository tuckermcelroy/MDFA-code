mdfa.cofilter <- function(frf,spec,noisefreqs,beta,R,Q)
{
  
  #######################################################
  #
  #	mdfa.cofilter by Tucker McElroy
  #
  #	computes optimal concurrent moving average filter
  #		as approx of given target filter, using MDFA method,
  #		based upon the moving average filter class of length q
  #	inputs:
  #		frf is array N x N x Grid of complex entries, the target
  #			frequency response function Psi(e^{-i lambda})
  #			for lambda given by Grid number of Fourier frequencies 
  #		spec is array N+1 x N x Grid of complex entries, f(lambda):
  #     the process/data cross spectrum of coint process with data process,
  #			followed by process/data spectral density matrix 
  #		noisefreqs lists frequencies lambda=pi*x, with x in [-1,1], 
  #			with repeats for double roots, 
  #			such that exp(-i*lambda) is a root
  #			of the noise differencing polynomial
  #   beta is N vector of co-integrating relations
  #		R is array N x q x N x M of q-M constraints on the moving 
  #			average filter solution.
  #		Q is array N x q x N of constraints, theta = R phi + Q.
  #	outputs:
  #		opt.array is array N x N x q of filter coefficients
  #
  ##############################################################
  
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
  
  N <- dim(spec)[1]
  Grid <- dim(frf)[3]
  m <- floor(Grid/2)
  q <- dim(R)[2]
  M <- dim(R)[4]
  R.mat <- matrix(R,nrow=N*q,ncol=N*M)
  Q.mat <- matrix(Q,ncol=N)
  
  noise.ceps <- Re(roots2ceps(exp(-1i*pi*noisefreqs),length(noisefreqs)))
  delta.noise <- ceps2wold(noise.ceps,length(noisefreqs))
  
  spec.zx <- array(0,c(1,N,Grid))
  spec.xx <- array(0,c(N,N,Grid))
  for(i in 1:Grid)
  {
    spec.xz[,,i] <- spec[1,,i,drop=FALSE]
    spec.xx[,,i] <- spec[2:(N+1),,i,drop=FALSE]
  }  
  
  fpsi <- NULL
  fmat <- NULL
  lambda.ft <- exp(-1i*2*pi*Grid^{-1}*(seq(1,Grid) - (m+1)))	## this is e^{-i lambda}
  noise.ft <- delta.noise[1]*lambda.ft^0
  for(j in 1:length(noisefreqs))
  {
    noise.ft <- noise.ft + delta.noise[j+1]*lambda.ft^j
  }  
  
  for(k in 0:(q-1))
  {
    fpsi.new <- do.call(cbind,lapply(seq(1,Grid),function(i) frf[,,i] %*% spec.xx[,,i]))
    fpsi.new <- Grid^{-1}*fpsi.new %*% (lambda.ft^{-k} %x% diag(N))
    fpsi.coint <- do.call(cbind,lapply(seq(1,Grid),function(i) spec.zx[,,i] - t(beta) %*% spec.xx[,,i] ))
    fpsi.coint <- Grid^{-1}*fpsi.coint %*% (noise.ft*lambda.ft^{-k} %x% diag(N))
    fpsi <- cbind(fpsi,fpsi.new + fpsi.coint)
    fmat.new <- Grid^{-1}*matrix(spec,nrow=N) %*% (lambda.ft^{-k} %x% diag(N))
    if(k==0) { 
      fmat <- fmat.new 
      fzero <- fmat.new
    } else {
      if(k==1) {
        fmat <- cbind(fmat,fmat.new)
        fmat <- rbind(fmat,cbind(t(fmat.new),fzero))
      } else {
        side.mat <- fmat[1:(dim(fmat)[2]-N),(dim(fmat)[2]+1-N):dim(fmat)[2],drop=FALSE]
        fmat <- cbind(fmat,rbind(fmat.new,side.mat))
        fmat <- rbind(fmat,cbind(t(fmat.new),t(side.mat),fzero))
      }
    }
  }
  fpsi <- Re(fpsi)
  fmat <- Re(fmat)
  
  opt <- solve(t(R.mat) %*% fmat %*% R.mat) %*% t(R.mat) %*% (t(fpsi) - fmat %*% Q.mat)
  opt <- R.mat %*% opt + Q.mat
  opt.array <- array(t(opt),c(N,N,q))
  
  return(opt.array) 
}



