mdfa.target2i1 <- function(filter,shift)
{

	#######################################################
	#
	#	mdfa.target2i1 by Tucker McElroy
	#
	#	takes given target filter and computes alternate form
	#		for I(1) input series
	#	inputs:
	#		filter is array N x N x L of real entries, 
	#			corresponding to filter coefficients of index
	#			-shift,...,L-1-shift.
	#		shift gives the integer offset; shift=0 for a causal filter,
	#			shift = L-1 for anti-causal filter,
	#			shift = (L-1)/2 for symmetric filter
      #     outputs:
  	#		out.filter is N x N x L-1 of real entries,
	#			corresponding to filter coefficients of index
	#			-shift,...,L-2-shift.
	##############################################################

	N <- dim(filter)[1]
	L <- dim(filter)[3]

	r <- shift
	s <- L-1-shift

	forward.coeff <- NULL
	if(s > 0) {
		forward.coeff <- array(0,c(N,N,s))
		forward.coeff[,,s] <- filter[,,L]
		if(s > 1) {
			for(l in (s-1):1) 
			{
				forward.coeff[,,l] <- forward.coeff[,,l+1] + filter[,,r+1+l]
			}
		}
	}

	backward.coeff <- NULL
	if(r > 0) {
		backward.coeff <- array(0,c(N,N,r))
		backward.coeff[,,r] <- filter[,,1]
		if(r > 1) {
			for(l in (r-1):1) 
			{
				backward.coeff[,,l] <- backward.coeff[,,l+1] + filter[,,r-l+1]
			}
		}
	}
	
	out.filter <- array(0,c(N,N,L-1))
	if(r > 0) {
		for(l in 1:r)
		{
			out.filter[,,l] <- backward.coeff[,,(r+1-l)]
		}
	} 
	if(s > 0) {
		for(l in 1:s)
		{
			out.filter[,,(r+l)] <- forward.coeff[,,l]
		}
	}

	return(out.filter)	
}
	