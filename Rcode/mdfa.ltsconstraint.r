mdfa.ltsconstraint <- function(frf,spec,q)
{

	#######################################################
	#
	#	mdfa.ltsconstraint by Tucker McElroy
	#
	#	computes optimal concurrent moving average filter
	#		as approx of given target filter, using MDFA method.
	#		uses mdfa.filter with level and time shift constraint
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
	R.mat <- diag(q-2) %x% diag(N)
	R.mat <- rbind( -t(1+seq(1,q-2)) %x% diag(N), R.mat)
	R.mat <- rbind( t(seq(1,q-2)) %x% diag(N), R.mat)
	lc <- Re(frf[,,m+1])	# level constraint
	delta <- 2*pi/Grid
	tsc <- Re(1i*(frf[,,m+2]-frf[,,m+1])/delta)	# approximate derivative
	ltsc <- -tsc + lc
	Q.mat <- matrix(0,nrow=N*(q-2),ncol=N)
	Q.mat <- rbind(t(tsc),Q.mat)
	Q.mat <- rbind(t(ltsc),Q.mat)

	R <- array(R.mat,c(N,q,N,q-2))
	Q <- array(Q.mat,c(N,q,N))

	out <- mdfa.filter(frf,spec,R,Q)

	return(out)
}


	