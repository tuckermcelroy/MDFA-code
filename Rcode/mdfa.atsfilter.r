mdfa.atsfilter <- function(frf,spec,lambda,eta,mu,q)
{
  
  #######################################################
  #
  #	mdfa.atsfilter by Tucker McElroy
  #
  #	computes optimal concurrent moving average filter
  #		as approx of given target filter, using I-ATS MDFA method,
  #		based upon the moving average filter class of length q
  #	inputs:
  #		frf is array N x N x Grid of complex entries, the target
  #			frequency response function Psi(e^{-i lambda})
  #			for lambda given by Grid number of Fourier frequencies 
  #		spec is array N x N x Grid of complex entries, the
  #			process/data spectral density matrix f(lambda)
  #   lambda is a non-negative vector that tunes the Timeliness component
  #   eta is a non-negative vector that tunes the Smoothness component
  #   mu is a scalar in (0,pi) that gives cutoff for stop-band  
  #   q integer order of MDFA moving average filter
  #	outputs:
  #		opt.array is array N x N x q of filter coefficients
  #		opt.val is N x N matrix corresponding to minimal MSE
  #
  ##############################################################
  
  ## first code up N=1 case...
  
  N <- dim(spec)[1]
  grid <- dim(frf)[3]
  m <- floor(grid/2)
  lambda.ft <- exp(-1i*2*pi*grid^{-1}*(seq(1,grid) - (m+1)))	## this is e^{-i lambda}
  
  W.fcn <- list()
  scaling <- list()
  for(i in 1:N)
  {
    W.fcn[[i]] <- (1 + pmax(0,abs(2*pi*grid^{-1}*(seq(1,grid) - (m+1))) - mu))^eta[i]
    Amp.frf <- abs(frf[i,i,,drop=FALSE])
    Phase.frf <- -Arg(frf[i,i,,drop=FALSE])
    scaling[[i]] <- sqrt(1 + 4*lambda[i]*Amp.frf^2)
  }

  
  fpsi <- NULL
  fmat <- NULL
  g.fcn <- list()
#  opt.val <- do.call(cbind,lapply(seq(1,grid),function(i) frf[,,i] %*% spec[,,i] %*% Conj(t(frf[,,i]))))
#  opt.val <- grid^{-1}*opt.val %*% (rep(1,grid) %x% diag(N))
  for(k in 0:(q-1))
  {
    g.piece <- exp(-1i*Phase.frf)*lambda.ft^k
    g.fcn[[k+1]] <- Re(g.piece) + 1i*scaling*Im(g.piece)
    
    
  }
  fpsi <- Re(fpsi)
  fmat <- Re(fmat)

  opt <- solve(fmat,fpsi)
  
  
}


