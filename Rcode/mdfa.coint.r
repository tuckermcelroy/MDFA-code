mdfa.coint <- function(frf,spec.main,spec.coint,alpha,q)
{

	#######################################################
	#
	#	mdfa.coint by Tucker McElroy
	#
	#	computes optimal concurrent moving average filter
	#		as approx of given target filter, using MDFA method,
	#		based upon the moving average filter class of length q
	#	inputs:
	#		frf is array N x N x Grid of complex entries, the target
	#			frequency response function Psi(e^{-i lambda})
	#			for lambda given by Grid number of Fourier frequencies 
	#		spec.main is array N x N x Grid of complex entries, the
	#			process/data spectral density matrix f(lambda)
  #		spec.coint is array r x N x Grid of complex entries, the
  #			cross-spectral density matrix of Delta^N (L)Z_t with Delta (L)X_t
  #   alpha: a N x r dimensional matrix, where r is number of
  #     co-integrating relations
  #   q: desired order of MA filter
	#	outputs:
	#		opt.array is array N x N x q of filter coefficients
	#
	##############################################################

	N <- dim(spec.main)[1]
	grid <- dim(frf)[3]
	m <- floor(grid/2)
#	q <- dim(R)[2]
#	M <- dim(R)[4]
#	R.mat <- matrix(R,nrow=N*q,ncol=N*M)
#	Q.mat <- matrix(Q,ncol=N)

	fpsi <- NULL
	fmat <- NULL
	lambda.ft <- exp(-1i*2*pi*grid^{-1}*(seq(1,grid) - (m+1)))	## this is e^{-i lambda}

#	opt.val <- do.call(cbind,lapply(seq(1,grid),function(i) frf[,,i] %*% spec.main[,,i] %*% Conj(t(frf[,,i]))))
#	opt.val <- grid^{-1}*opt.val %*% (rep(1,grid) %x% diag(N))
	for(k in 0:(q-1))
	{
		fpsi.new <- do.call(cbind,lapply(seq(1,grid),
                       function(i) frf[,,i] %*% spec.main[,,i] + alpha %*% spec.coint[,,i]))
 		fpsi.new <- grid^{-1}*fpsi.new %*% (lambda.ft^{-k} %x% diag(N))
		fpsi <- cbind(fpsi,fpsi.new)
		fmat.new <- grid^{-1}*matrix(spec.main,nrow=N) %*% (lambda.ft^{-k} %x% diag(N))
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
 
	opt <- solve(fmat,t(fpsi))
	opt.array <- array(t(opt),c(N,N,q))

	return(opt.array) 
}



