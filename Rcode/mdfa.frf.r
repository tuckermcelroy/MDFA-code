mdfa.frf <- function(filter,shift,grid)
{

	#######################################################
	#
	#	mdfa.frf by Tucker McElroy
	#
	#	computes  frequency response function (frf) of a given
	#		filter with given time shift, at Grid Fourier frequencies
	#	inputs:
	#		filter is array N x N x L of real entries, 
	#			corresponding to filter coefficients of index
	#			-shift,...,L-1-shift.
	#		shift gives the integer offset; shift=0 for a causal filter
	#		grid is integer number of Fourier frequencies	
	#	outputs:
	#		frf is array N x N x grid of complex entries
	#
	################################################# 

	N <- dim(filter)[1]
	L <- dim(filter)[3]
	m <- floor(grid/2)
	lambda.ft <- exp(-1i*2*pi*grid^(-1)*(seq(1,grid) - (m+1)))	## this is e^{-i lambda}

	frf.mat <- matrix(0,nrow=N,ncol=(N*grid))
	for(l in 1:L)
	{
		frf.mat <- frf.mat + (matrix(lambda.ft^(l-1-shift),nrow=1) %x% filter[,,l]) 
	}
	frf <- array(frf.mat,c(N,N,grid))

	return(frf)
}


