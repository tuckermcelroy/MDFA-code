mdfa.embed <- function(x.series,embed.f)
{
  
  #######################################################
  #
  #	mdfa.embed by Tucker McElroy
  #
  #  embeds a given time series at another sampling frequency
  #	inputs:
  #   x.series is a T x N multivariate time series of frequency 
  #      given by samp.f
  #   embed.f is sampling frequency of embedding
  #	outputs:
  #   x.mat is T/s x Ns, where s = samp.f/embed.f 
  #
  #  NOTE: s must be integer
  ##############################################################
 
  N <- dim(x.series)[2]
  T <- dim(x.series)[1]
  samp.f <- frequency(x.series)
  s.embed <- samp.f/embed.f
  T <- s.embed * (T %/% s.embed)
  x.mat <- NULL
  for(j in 1:N)
  {
    x.svf <- t(matrix(x.series[1:T,j,drop=FALSE],nrow=s.embed))
    x.mat <- cbind(x.mat,x.svf)
  }  
  return(x.mat) 
}