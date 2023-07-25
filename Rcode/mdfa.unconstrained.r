mdfa.unconstrained <- function(frf,spec,q)
{

	#######################################################
	#
	#	mdfa.unconstrained by Tucker McElroy
	#
	#	computes optimal concurrent moving average filter
	#		as approx of given target filter, using MDFA method.
	#		uses mdfa.filter with no constraints
	#	inputs:
	#		frf is array N x N x Grid of complex entries, the target
	#			frequency response function Psi(e^{-i lambda})
	#			for lambda given by Grid number of Fourier frequencies 
	#		spec is array N x N x Grid of complex entries, the
	#			process/data spectral density matrix f(lambda)
  #   q: desired order of MA filter
  # outputs:
  #   opt is list consisting of:
  #   opt.array is array N x N x q of filter coefficients
  #   opt.val is N x N matrix corresponding to minimal MSE
 	#
	##############################################################

	N <- dim(frf)[1]
	R.mat <- diag(q) %x% diag(N)
	Q.mat <- matrix(0,nrow=N*q,ncol=N)

	R <- array(R.mat,c(N,q,N,q))
	Q <- array(Q.mat,c(N,q,N))

	out <- mdfa.filter(frf,spec,R,Q)

	return(out)
}


	