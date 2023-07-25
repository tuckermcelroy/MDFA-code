mdfa.dft <- function(x.data)
{

	#############################################
	#
	#	mdfa.dft by Tucker McElroy
	#
	#	computes multivariate dft from x.data,
	#		at Fourier frequencies
	#	inputs: 
	#		x.data is T x N matrix of reals
	#			There are T Fourier frequencies, space 2pi/T apart
	#	outputs: 
	#		x.dft is T x N matrix of complex numbers
	#
	#################################################

	T <- dim(x.data)[1]
	m <- floor(T/2)

	Q.mat <- exp(-1i*2*pi*T^{-1}*matrix(seq(1,T) - (m+1),ncol=1) %x% matrix(seq(1,T),nrow=1))*T^{-1/2}
	x.dft <- Q.mat %*% x.data

	return(x.dft)
}
