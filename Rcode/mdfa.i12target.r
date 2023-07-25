mdfa.i12target <- function(filter,level)
{

	#######################################################
	#
	#	mdfa.i12target by Tucker McElroy
	#
	#	takes given concurrent filter in its I(1) form, and
	#		computes alternate form, as a inverse of target2i1
	#	inputs:
	#		filter is array N x N x L-1 of real entries, 
	#			corresponding to filter coefficients of index
	#			0,...,L-2.
	#		level equals sum of original target filter coefficients,
	#			and is added to the first coefficient
      #     outputs:
  	#		out.filter is N x N x L of real entries,
	#			corresponding to filter coefficients of index
	#			0,...,L-1.
	##############################################################

	N <- dim(filter)[1]
	L <- dim(filter)[3]+1

	out.filter <- array(0,c(N,N,L))
	out.filter[,,1] <- -1*filter[,,1] + level
	for(l in 2:(L-1))
	{
		out.filter[,,l] <- filter[,,l-1] - filter[,,l]
	}
	out.filter[,,L] <- filter[,,L-1] 

	return(out.filter)	
}
	