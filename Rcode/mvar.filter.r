mvar.filter <- function(data.ts,psi.filter)
{
  #######################################################
  #
  #	mvar.filter by Tucker McElroy
  #
  #	computes multivariate filter output
  #	inputs:
  #   data.ts is T x N time series
  #		psi.filter is array N x N x L of real entries,
  #			corresponding to filter coefficients of index
  #			-shift,...,L-1-shift, where shift gives the integer offset;
  #      shift=0 for a causal filter
  #	outputs:
  #		out.ts is array (T-L+1) x N output series
  #
  #################################################

  N <- dim(data.ts)[2]
  T <- dim(data.ts)[1]
  L <- dim(psi.filter)[3]

  out.ts <- NULL
  for(j in 1:N)
  {
    output.j <- rep(0,T-L+1)
    for(k in 1:N)
    {
      output.k <- stats::filter(data.ts[,k],psi.filter[j,k,],
                       method="convolution",sides=1)
      output.k <- output.k[L:T]
      output.j <- output.j + output.k
    }
    out.ts <- cbind(out.ts,output.j)
  }

  return(out.ts)
}
