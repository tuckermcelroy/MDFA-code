mdfa.atsfilter2 <- function(frf,spec,lambda1,lambda2,eta,mu,q)
{
  
  #######################################################
  #
  #	mdfa.atsfilter2 by Tucker McElroy
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
  #   lambda1 is a non-negative vector that tunes the Timeliness component
  #   lambda2 is a non-negative vector that tunes the Smoothness component
  #   eta is a non-negative vector that further tunes the Smoothness component
  #   mu is a scalar in (0,pi) that gives cutoff for stop-band  
  #   q integer order of MDFA moving average filter
  #	outputs:
  #		opts is list of N univariate length q  filter coefficients
  #		opt.vals is list of N scalars corresponding to minimal MSE
  #
  ##############################################################
  
  N <- dim(spec)[1]
  grid <- dim(frf)[3]
  m <- floor(grid/2)
  freq.ft <- 2*pi*grid^{-1}*(seq(1,grid) - (m+1))
  lambda.ft <- exp(-1i*freq.ft)	## this is e^{-i lambda}

  pass.frf <- rep(1,grid)
  pass.frf[abs(freq.ft) > mu] <- 0
  stop.frf <- 1 - pass.frf

  opts <- NULL
  opt.vals <- NULL
  for(n in 1:N)
  {
    W.fcn <- (1 + pmax(0,abs(2*pi*grid^{-1}*(seq(1,grid) - (m+1))) - mu))^eta[n]
    opt.val <- mean(abs(frf[n,n,])^2 * spec[n,n,])
    fpsi <- NULL
    fmat <- NULL
    for(k in 1:q)
    {
      fpsi.new <- mean(frf[n,n,] * spec[n,n,] * lambda.ft^{1-k})
      fpsi <- cbind(fpsi,fpsi.new)
      fmat.col <- NULL
      for(j in 1:q)
      {
        fmat.new <- mean((1 + W.fcn * lambda2[n] * stop.frf) * lambda.ft^{j-k} * spec[n,n,])
        spec.temp <- mean((W.fcn * lambda1[n] * pass.frf) * lambda.ft^{-j-k} * frf[n,n,]^2 * spec[n,n,])
        fmat.new <- fmat.new - spec.temp
        spec.temp <- mean((W.fcn * lambda1[n] * pass.frf) * lambda.ft^{j+k} * Conj(frf[n,n,])^2 * spec[n,n,])
        fmat.new <- fmat.new - spec.temp
        spec.temp <- mean((W.fcn * lambda1[n] * pass.frf) * lambda.ft^{j-k} * Mod(frf[n,n,])^2 * spec[n,n,])
        fmat.new <- fmat.new + 2*spec.temp
        fmat.col <- rbind(fmat.col,fmat.new)
      }
      fmat <- cbind(fmat,fmat.col)
    }  
    fpsi <- Re(fpsi)
    fmat <- Re(fmat)
    opt <- solve(fmat,t(fpsi))
    opt.val <- Re(opt.val) - fpsi %*% opt
    
    opts[[n]] <- opt
    opt.vals[[n]] <- opt.val
  }    
  
  return(list(opts,opt.vals))
}


