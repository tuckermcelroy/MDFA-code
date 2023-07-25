mdfa.coeff <- function(frf,aft,fore)
{

	#######################################################
	#
	#	mdfa.coeff by Tucker McElroy
	#
	#	computes filter coefficients from a given
	#		frequency response function (frf), for indices
	#		aft through fore
	#	inputs:
	#		frf is array N x N x grid of complex entries
  #		aft <= fore are integer indices
	#	outputs:
	#		filter is array N x N x (fore-aft+1) of real entries
	#
	################################################# 

	N <- dim(frf)[1]
	grid <- dim(frf)[3]
	m <- floor(grid/2)
	lambda.ft <- exp(-1i*2*pi*grid^{-1}*(seq(1,grid) - (m+1)))	## this is e^{-i lambda}
	filter <- array(0,c(N,N,(fore-aft+1)))

	for(l in aft:fore)
	{
		val <- do.call(cbind,lapply(seq(1,grid),function(i) frf[,,i] * lambda.ft[i]^{-l}))
		val <- grid^{-1}*val %*% (rep(1,grid) %x% diag(N))
		filter[,,l-aft+1] <- Re(val)
	}

	return(filter)
}
