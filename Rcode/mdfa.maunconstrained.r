mdfa.maunconstrained <- function(psi.array,acf.array,q)
{

	#######################################################
	#
	#	mdfa.maunconstrained by Tucker McElroy
	#
	#	computes optimal concurrent moving average filter
	#		as approx of given target filter, using MDFA method.
	#		uses mdfa.mafilter with no constraints
	#	inputs:
	#		psi.array: the target filter at lags r+q-2 through 1-r,
	#			as an N x N x (2r+q-2) array object
	#		acf.array: given acf at lags 0 through r-1,
	#			as an N x N x r array object
	#		q: desired order of MA filter
	#	outputs:
	#		opt is list consisting of:
	#		opt.array is array N x N x q of filter coefficients
	#		opt.val is N x N matrix corresponding to minimal MSE
	#
	##############################################################

	N <- dim(psi.array)[2]
	
	R.mat <- diag(q) %x% diag(N)
	Q.mat <- matrix(0,nrow=N*q,ncol=N)

	R <- array(R.mat,c(N,q,N,q))
	Q <- array(Q.mat,c(N,q,N))

	opt <- mdfa.mafilter(psi.array,acf.array,R,Q)

	return(opt)
}


	