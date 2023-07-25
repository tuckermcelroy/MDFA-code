mdfa.tsconstraint <- function(frf,spec,q)
{

	#######################################################
	#
	#	mdfa.tsconstraint by Tucker McElroy
	#
	#	computes optimal concurrent moving average filter
	#		as approx of given target filter, using MDFA method.
	#		uses mdfa.filter with time shift constraint
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
	R.mat[1:N,] <- R.mat[1:N,] - t(seq(1,q-1)) %x% diag(N)
	R.mat <- rbind( t(c(1,rep(0,q-2))) %x% diag(N), R.mat)
	delta <- 2*pi/Grid
	tsc <- Re(1i*(frf[,,m+2]-frf[,,m+1])/delta)	# approximate derivative
	Q.mat <- matrix(0,nrow=N*(q-2),ncol=N)
	Q.mat <- rbind(t(tsc),Q.mat)
	Q.mat <- rbind(0*diag(N),Q.mat)
	
	R <- array(R.mat,c(N,q,N,q-1))
	Q <- array(Q.mat,c(N,q,N))

	out <- mdfa.filter(frf,spec,R,Q)

	return(out)
}


	