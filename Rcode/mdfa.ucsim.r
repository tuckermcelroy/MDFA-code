mdfa.ucsim <- function(delta,innovations)
{
  
  ###########################
  #	mdfa.ucsim
  #		by Tucker McElroy
  #
  #	generates (T-d) x N simulations from input innovations,
  #   which are T x N dimensional white noise.
  #  delta is the degree d differencing polynomial,
  #   as array N x N x (d+1)
  #		
  ##############################
  
  N <- dim(innovations)[2]
  T <- dim(innovations)[1]
  d <- dim(delta)[3] - 1
  delta.0 <- delta[,,1]
  delta.mat <- matrix(-1*delta[,,-1],nrow=N)
  sims <- solve(delta.0) %*% t(innovations[d:1,,drop=FALSE])
  for(t in (d+1):T)
  {
    new.sim <- solve(delta.0) %*% delta.mat %*% matrix(sims[,1:d],ncol=1) +
        solve(delta.0) %*% t(innovations[t,,drop=FALSE])
    sims <- cbind(new.sim,sims)
  }
  sims <- t(sims[,seq(T-d,1)])
                  
  return(sims)
}


