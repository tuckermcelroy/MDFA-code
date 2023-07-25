mdfa.i22target <- function(filter,level,ts)
{

	#######################################################
	#
	#	mdfa.i22target by Tucker McElroy
	#
	#	takes given concurrent filter in its I(2) form, and
	#		computes alternate form, as a inverse of target2i2
	#	inputs:
	#		filter is array N x N x L-2 of real entries, 
	#			corresponding to filter coefficients of index
	#			0,...,L-3.
	#		level equals sum of original target filter coefficients,
	#			and is added to the first coefficient
	#		ts equals time shift from original target filter coefficients
      #     outputs:
  	#		out.filter is N x N x L of real entries,
	#			corresponding to filter coefficients of index
	#			0,...,L-1.
	##############################################################

	N <- dim(filter)[1]
	L <- dim(filter)[3]+2

	out.filter <- array(0,c(N,N,L))
	out.filter[,,1] <- filter[,,1] + level - ts
	out.filter[,,2] <- filter[,,2] - 2*filter[,,1] + ts
	for(l in 3:(L-2))
	{
		out.filter[,,l] <- filter[,,l] - 2*filter[,,l-1] + filter[,,l-2]
	}
	out.filter[,,L-1] <- -2*filter[,,L-2] + filter[,,L-3]
	out.filter[,,L] <- filter[,,L-2] 

	return(out.filter)	
}
	