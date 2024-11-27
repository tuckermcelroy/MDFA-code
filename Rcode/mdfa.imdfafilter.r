mdfa.imdfafilter <- function(frf,spec,lambda,eta,mu,q)
{
  
  #######################################################
  #
  #	mdfa.imdfafilter by Tucker McElroy
  #
  #	computes optimal concurrent moving average filter
  #		as approx of given target filter, using I-ATS MDFA method,
  #		based upon the moving average filter class of length q.
  #   This code is applied multivariately, one component at a time
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
  #		opts is list of N univariate length q  filter coefficients
  #		opt.vals is list of N scalars corresponding to minimal MSE
  #
  ##############################################################
  
  N <- dim(spec)[1]
  grid <- dim(frf)[3]
  m <- floor(grid/2)
  lambda.ft <- exp(-1i*2*pi*grid^{-1}*(seq(1,grid) - (m+1)))	## this is e^{-i lambda}

  opts <- NULL
  opt.vals <- NULL
  for(n in 1:N)
  {
    W.fcn <- (1 + pmax(0,abs(2*pi*grid^{-1}*(seq(1,grid) - (m+1))) - mu))^eta[n]
    opt.val <- do.call(cbind,lapply(seq(1,grid),
                            function(i) matrix(frf[n,,i],nrow=1) %*% spec[,,i] %*% 
                                        Conj(t(matrix(frf[n,,i],nrow=1)))))
	  opt.val <- grid^{-1}*opt.val %*% W.fcn 
		fpsi <- NULL
	  fmat <- NULL
	  for(k in 1:q)
	  {
		  fpsi.new <- do.call(cbind,lapply(seq(1,grid),
		                        function(i) matrix(frf[n,,i],nrow=1) %*% spec[,,i]))
 		  fpsi.new <- grid^{-1}*fpsi.new %*% ((W.fcn * lambda.ft^{1-k}) %x% diag(N))
		  fpsi <- cbind(fpsi,fpsi.new)
		  fmat.col <- NULL
		  for(j in 1:q)
		  {
		    fmat.new <- grid^{-1}*matrix(spec,nrow=N) %*% ((W.fcn * lambda.ft^{j-k}) %x% diag(N))
		    spec.temp <- do.call(cbind,lapply(seq(1,grid),
		                          function(i) t(spec[,,i]) %*% t(matrix(frf[n,,i],nrow=1)) %*% 
		                                      matrix(frf[n,,i],nrow=1) %*% spec[,,i] ))
		    fmat.temp <- grid^{-1}*spec.temp %*% ((W.fcn * lambda.ft^{-j-k}) %x% diag(N))
		    fmat.new <- fmat.new - lambda[n]*fmat.temp
		    spec.temp <- do.call(cbind,lapply(seq(1,grid),
		                          function(i) spec[,,i] %*% Conj(t(matrix(frf[n,,i],nrow=1))) %*% 
		                                      Conj(matrix(frf[n,,i],nrow=1)) %*% t(spec[,,i]) ))
		    fmat.temp <- grid^{-1}*spec.temp %*% ((W.fcn * lambda.ft^{j+k}) %x% diag(N))
		    fmat.new <- fmat.new - lambda[n]*fmat.temp
		    spec.temp <- do.call(cbind,lapply(seq(1,grid),
		                          function(i) spec[,,i] %*% Conj(t(matrix(frf[n,,i],nrow=1))) %*% 
		                                      matrix(frf[n,,i],nrow=1) %*% spec[,,i] ))
		    fmat.temp <- grid^{-1}*spec.temp %*% ((W.fcn * lambda.ft^{j-k}) %x% diag(N))
		    fmat.new <- fmat.new + 2*lambda[n]*fmat.temp
		    fmat.col <- rbind(fmat.col,fmat.new)
		  }
		  fmat <- cbind(fmat,fmat.col)
	  }  
		fpsi <- Re(fpsi)
    fmat <- Re(fmat)
    opt <- solve(fmat,t(fpsi))
    opt.val <- opt.val - fpsi %*% opt
  
    opts <- rbind(opts,opt)
    opt.vals <- c(opt.vals,opt.val)  
  }    
    
  return(list(opts,opt.vals))
}


