mdfa.pergram <- function(x.data,delta)
{

	#######################################################
	#
	#	mdfa.pergram by Tucker McElroy
	#
	#	computes uncentered periodogram from x.data,
	#		at Fourier frequencies
	#	inputs:
	#		x.data is T x N matrix of reals
	#			There are T Fourier frequencies, space 2pi/T apart
	#		delta is differencing polynomial for pseudo-spectrum:
	#			delta=1 for stationary, delta=c(1,-1) for I(1), 
	#			delta=c(1,-2,1) for I(2), etc
	#	output:
	#		periodogram is N x N x T array of complex entries
	#
	###################################################

	T <- dim(x.data)[1]
	N <- dim(x.data)[2]
 	Grid <- T 
	m <- floor(Grid/2)
	lambda.ft <- exp(-1i*2*pi*Grid^(-1)*(seq(1,Grid) - (m+1)))	## this is e^{-i lambda}
	d <- length(delta)-1
	diff.op <- delta[1]*lambda.ft^0
	if(d>0) {
	for(i in 2:(d+1))
	{
		diff.op <- diff.op + delta[i]*lambda.ft^(i-1)
	} }
	diff.op <- Mod(diff.op)^2
	
	x.dft <- mdfa.dft(x.data)
	pergram <- do.call(cbind,lapply(seq(1,T),function(i) x.dft[i,] %*% Conj(t(x.dft[i,])) * diff.op[i]^(-1)))
	pergram <- array(pergram,c(N,N,T))
#	pergram[,,m+1] <- 0*diag(N)
	return(pergram)
}

