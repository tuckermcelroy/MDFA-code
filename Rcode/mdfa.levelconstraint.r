mdfa.levelconstraint <- function(frf,spec,q)
{

	#######################################################
	#
	#	mdfa.levelconstraint by Tucker McElroy
	#
	#	computes optimal concurrent moving average filter
	#		as approx of given target filter, using MDFA method.
	#		uses mdfa.filter with level constraint
	#	inputs:
	#		frf is array N x N x Grid of complex entries, the target
	#			frequency response function Psi(e^{-i lambda})
	#			for lambda given by Grid number of Fourier frequencies 
	#		spec is array N x N x Grid of complex entries, the
	#			process/data spectral density matrix f(lambda)
      #           q: desired order of MA filter
	#	outputs:
      #           opt is list consisting of:
      #           opt.array is array N x N x q of filter coefficients
      #           opt.val is N x N matrix corresponding to minimal MSE
 	#
	##############################################################

	N <- dim(frf)[1]
	Grid <- dim(frf)[3]
	m <- floor(Grid/2)
	R.mat <- diag(q-1) %x% diag(N)
	R.mat <- rbind(-t(rep(1,q-1)) %x% diag(N),R.mat)
	Q.mat <- matrix(0,nrow=N*(q-1),ncol=N)
	lc <- Re(frf[,,m+1])	# level constraint
	Q.mat <- rbind(t(lc),Q.mat)

	R <- array(R.mat,c(N,q,N,q-1))
	Q <- array(Q.mat,c(N,q,N))

	out <- mdfa.filter(frf,spec,R,Q)

	return(out)
}


	