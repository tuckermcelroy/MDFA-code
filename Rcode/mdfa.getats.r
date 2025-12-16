mdfa.getats <- function(frf.psi,frf.hat,spec,mu)
{
  
  #######################################################
  #
  #	mdfa.getats by Tucker McElroy
  #
  #	computes Accuracy, Timeliness, Smoothness, and Residual components
  #		for two filters, the target (frf.psi) and the approximation (frf.hat).
  #   Works univariately for each of the series.
  #	inputs:
  #		frf.psi is array N x N x Grid of complex entries, the target
  #			frequency response function Psi(e^{-i lambda})
  #			for lambda given by Grid number of Fourier frequencies 
  #		frf.hat is array N x N x Grid of complex entries, the approximation
  #			frequency response function \hat{Psi} (e^{-i lambda})
  #			for lambda given by Grid number of Fourier frequencies 
  #		spec is array N x N x Grid of complex entries, the
  #			process/data spectral density matrix f(lambda)
  #   mu is a scalar in (0,pi) that gives cutoff for stop-band  
  #	outputs:
  #
  ##############################################################
  
  N <- dim(frf.psi)[1]
  grid <- dim(frf.psi)[3]
  m <- floor(grid/2)
  freq.ft <- 2*pi*grid^{-1}*(seq(1,grid) - (m+1))
  
  pass.frf <- rep(1,grid)
  pass.frf[abs(freq.ft) > mu] <- 0
  stop.frf <- 1 - pass.frf
  
  psi.amps <- NULL
  psi.phases <- NULL
  hat.amps <- NULL
  hat.phases <- NULL
  ats <- NULL
  for(i in 1:N)
  {
    psi.amp <- abs(frf.psi[i,i,])
    psi.phase <- -Arg(frf.psi[i,i,])
    hat.amp <- abs(frf.hat[i,i,])
    hat.phase <- -Arg(frf.hat[i,i,])
    psi.amps <- cbind(psi.amps,psi.amp)
    psi.phases <- cbind(psi.phases,psi.phase)
    hat.amps <- cbind(hat.amps,hat.amp)
    hat.phases <- cbind(hat.phases,hat.phase)
    
    accur <-  Re(mean((psi.amp - hat.amp)^2*pass.frf*spec[i,i,]))
    timely <- Re(4*mean(psi.amp*hat.amp*(sin((psi.phase - hat.phase)/2))^2*
                         pass.frf*spec[i,i,]))
    smooth <-  Re(mean((psi.amp - hat.amp)^2*stop.frf*spec[i,i,]))
    resi <- Re(4*mean(psi.amp*hat.amp*(sin((psi.phase - hat.phase)/2))^2*
                       stop.frf*spec[i,i,]))
    ats <- cbind(ats,c(accur,timely,smooth,resi))
  }
  
  return(list(psi.amps,psi.phases,hat.amps,hat.phases,ats))
  
}  
  