mdfa.imdfafilter <- function(frf,spec,lambda,eta,mu,q)
{
  
  #######################################################
  #
  #	mdfa.imdfafilter by Tucker McElroy
  #
  #	computes optimal concurrent moving average filter
  #		as approx of given target filter, using I-ATS MDFA method,
  #		based upon the moving average filter class of length q.
  #   This code is applied multivariately
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

  W.fcn <- (1 + pmax(0,abs(2*pi*grid^{-1}*(seq(1,grid) - (m+1))) - mu))^eta
     
  fpsi <- NULL
	fmat <- NULL
	lambda.ft <- exp(-1i*2*pi*grid^{-1}*(seq(1,grid) - (m+1)))	## this is e^{-i lambda}

	opt.val <- do.call(cbind,lapply(seq(1,grid),function(i) W.fcn[i] * 
	                                  frf[,,i] %*% spec[,,i] %*% Conj(t(frf[,,i]))))
	opt.val <- grid^{-1}*opt.val %*% (rep(1,grid) %x% diag(N))
	for(k in 0:(q-1))
	{
		fpsi.new <- do.call(cbind,lapply(seq(1,grid),function(i) W.fcn[i] * frf[,,i] %*% spec[,,i]))
 		fpsi.new <- grid^{-1}*fpsi.new %*% (lambda.ft^{-k} %x% diag(N))
		fpsi <- cbind(fpsi,fpsi.new)
		fmat.new <- grid^{-1}*matrix(spec,nrow=N) %*% (lambda.ft^{-k} %x% diag(N))
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
  opt <- solve(fmat,fpsi)
  opt.val <- opt.val - t(fpsi) %*% opt
    
  return(list(opts,opt.vals))
}


