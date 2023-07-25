mdfa.wkfrf <- function(delta.noise,delta.signal,spec.noise,spec.signal)
{
  ########################################
  ## mdfa.wkfrf
  #   computes model-based  frequency response function
  #     for optimal symmetric (WK) signal extraction filter
  #   delta.noise: N x N x (d.noise+1) array of noise differencing
  #   delta.signal: N x N x (d.signal+1) array of signal differencing
  #		spec.noise: array of dimension N x N x grid, consisting of 
  #     spectrum for differenced noise process at Fourier frequencies
  #		spec.signal: array of dimension N x N x grid, consisting of 
  #     spectrum for differenced signal process at Fourier frequencies
  #	Outputs:
  #   frf.wk: array of dimension N x N x grid of frf, where
  #		  grid is the desired number of frequencies 
  #  Note: if no delta.noise or delta.signal portions, pass in identity matrices!
  #  Conventions: normal polynomial conventions for differencing polynomials
  #
  # Requires: getGCD  [in MDFA code, this GCD handles complex matrices]
  ###########################################################################
  
  d.noise <- dim(delta.noise)[3] - 1
  d.signal <- dim(delta.signal)[3] - 1
  N <- dim(spec.noise)[1]
  grid <- dim(spec.noise)[3]
  m <- floor(grid/2)
  lambda <- 2*pi*(seq(1,grid) - (m+1))/grid
  
  frf.noise <- t(rep(1,grid)) %x% delta.noise[,,1]  
  if(d.noise > 0) {
    for(i in 1:d.noise)
    {
      frf.noise <- frf.noise + t(exp(-1i*lambda*i)) %x% delta.noise[,,i+1]  
    } }
  frf.noise <- array(frf.noise,c(N,N,grid))
  
  frf.signal <- t(rep(1,grid)) %x% delta.signal[,,1]  
  if(d.signal > 0) {
    for(i in 1:d.signal)
    {
      frf.signal <- frf.signal + t(exp(-1i*lambda*i)) %x% delta.signal[,,i+1]  
    } }
  frf.signal <- array(frf.signal,c(N,N,grid))
  
  spec.signal.del <- array(t(rep(1,grid) %x% diag(N)),c(N,N,grid))
  spec.noise.del <- array(t(rep(1,grid) %x% diag(N)),c(N,N,grid))
  spec.data.del <- array(t(rep(1,grid) %x% diag(N)),c(N,N,grid))
  frf.wk <- array(t(rep(1,grid) %x% diag(N)),c(N,N,grid))
  for(k in 1:grid)
  {
    spec.noise.del[,,k] <- frf.signal[,,k] %*% spec.noise[,,k] %*% Conj(t(frf.signal[,,k]))
    spec.signal.del[,,k] <- frf.noise[,,k] %*% spec.signal[,,k] %*% Conj(t(frf.noise[,,k]))
    spec.data.del[,,k] <- spec.noise.del[,,k] + spec.signal.del[,,k]
    if( sum(Mod(frf.noise[,,k])) == 0 )
    {
      spec.rank <- qr(spec.noise[,,k])$rank
      gcd.out <- getGCD(spec.noise[,,k],spec.rank)
      spec.chol <- gcd.out[[1]] %*% diag(sqrt(gcd.out[[2]]),ncol=spec.rank)
      spec.sig.inv <- solve(spec.signal[,,k])
      frf.wk[,,k] <- diag(N) - spec.chol %*% solve( as.matrix(t(Conj(spec.chol)) %*% 
                      spec.sig.inv %*% spec.chol) ) %*% t(Conj(spec.chol)) %*% spec.sig.inv 
    }  else 
    {
      if( sum(Mod(frf.signal[,,k])) > 0 )
      {
        frf.wk[,,k] <- spec.signal.del[,,k] %*% solve(spec.data.del[,,k])
      } else
      {
        spec.rank <- qr(spec.signal[,,k])$rank
        gcd.out <- getGCD(spec.signal[,,k],spec.rank)
        spec.chol <- gcd.out[[1]] %*% diag(sqrt(gcd.out[[2]]),ncol=spec.rank)
        spec.noise.inv <- solve(spec.noise[,,k])
        frf.wk[,,k] <- spec.chol %*% solve( as.matrix(t(Conj(spec.chol)) %*% 
                        spec.noise.inv %*% spec.chol) ) %*% t(Conj(spec.chol)) %*% 
                        spec.noise.inv 
      }
    }
  }
  
  return(frf.wk)
}