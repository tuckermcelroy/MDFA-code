mdfa.atsfilter <- function(frf,spec,lambda,eta,mu,q)
{
  
  #######################################################
  #
  #	mdfa.atsfilter by Tucker McElroy
  #
  #	computes optimal concurrent moving average filter
  #		as approx of given target filter, using I-ATS MDFA method,
  #		based upon the moving average filter class of length q.
  #   This code is applied univariately, to each of N series.
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
  
  N <- dim(spec)[1]
  grid <- dim(frf)[3]
  m <- floor(grid/2)
  lambda.ft <- exp(-1i*2*pi*grid^{-1}*(seq(1,grid) - (m+1)))	## this is e^{-i lambda}
  
  opts <- list()
  opt.vals <- list()
  for(i in 1:N)
  {
    W.fcn <- (1 + pmax(0,abs(2*pi*grid^{-1}*(seq(1,grid) - (m+1))) - mu))^eta[i]
    Amp.frf <- abs(frf[i,i,])
    Phase.frf <- -Arg(frf[i,i,])
    scaling <- sqrt(1 + 4*lambda[i]*Amp.frf^2)
    
    g.fcn <- NULL
    opt.val <- mean(W.fcn*Amp.frf^2*spec[i,i,])
    for(k in 0:(q-1))
    {
      g.piece <- exp(-1i*Phase.frf)*lambda.ft^k
      g.fcn <- rbind(g.fcn,Re(g.piece) + 1i*scaling*Im(g.piece))
    }

    fpsi <- matrix(0,nrow=q,ncol=1)
    fmat <- matrix(0,nrow=q,ncol=q)
    for(j in 0:(q-1))
    {
      fpsi[(j+1),1] <- mean(g.fcn[j+1,]*W.fcn*Amp.frf*spec[i,i,])
      for(k in 0:(q-1))
      {
        fmat[(j+1),(k+1)] <- mean(g.fcn[j+1,]*Conj(g.fcn[k+1,])*W.fcn)
      }
    }
    fpsi <- Re(fpsi)
    fmat <- Re(fmat)
    opt <- solve(fmat,fpsi)
    opt.val <- opt.val - t(fpsi) %*% opt
  
    opts[[i]] <- opt
    opt.vals[[i]] <- opt.val
  }
  
  return(list(opts,opt.vals))
}


