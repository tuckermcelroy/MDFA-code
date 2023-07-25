mdfa.idft <- function(x.dft)
{

	#############################################
	#
	#	mdfa.idft by Tucker McElroy
	#
	#	computes inverse multivariate dft from x.dft
	#		at Fourier frequencies
	#	inputs: 
	#		x.dft is T x N matrix of complex dft
	#			There are T Fourier frequencies, space 2pi/T apart
	#	outputs: 
	#		x.data is T x N matrix of real numbers
	#
	#################################################

	T <- dim(x.dft)[1]
	m <- floor(T/2)

	Q.mat <- exp(-1i*2*pi*T^{-1}*matrix(seq(1,T) - (m+1),ncol=1) %x% matrix(seq(1,T),nrow=1))*T^{-1/2}
	x.data <- t(Conj(Q.mat)) %*% x.dft

	return(Re(x.data))
}

