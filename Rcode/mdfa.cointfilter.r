mdfa.cointfilter <- function(frf,spec,q,deltad)
{
  
  #######################################################
  #
  #	mdfa.cointfilter by Tucker McElroy
  #
  #	computes optimal concurrent moving average filter
  #		as approx of given target filter, using MDFA method,
  #		based upon the moving average filter class of length q,
  #   for the co-integrated estimation problem
  #	inputs:
  #		frf is array N x N x Grid of complex entries, the target
  #			frequency response function Psi^{sharp} (e^{-i lambda})
  #			for lambda given by Grid number of Fourier frequencies 
  #		spec is array (N+r) x (N+r) x Grid of complex entries, the
  #			process/data spectral density matrix f(lambda), in block form
  #     [ f_zz, f_zx ]
  #     [ f_xz, f_xx ]
  #     where f_zz is rxr for the noise-differenced co-integrated process Z_t,
  #     and f_xx is NxN for the fully-differenced process X_t
  #   q is length of moving average filter
  #   deltad is final coefficient of Delta(L)
  #	outputs:
  #		opt.array is array N x N x q of filter coefficients
  #		opt.val is N x N matrix corresponding to minimal MSE
  #
  ##############################################################
  
  N <- dim(frf)[1]
  grid <- dim(frf)[3]
  r <- dim(spec)[1] - N
  m <- floor(grid/2)
  spec.zz <- spec[1:r,1:r,,drop=FALSE]
  spec.xz <- spec[(r+1):(r+N),1:r,,drop=FALSE]
  spec.xx <- spec[(r+1):(r+N),(r+1):(r+N),,drop=FALSE]

  lambda.ft <- exp(-1i*2*pi*grid^{-1}*(seq(1,grid) - (m+1)))	## this is e^{-i lambda}
  fpsi <- do.call(cbind,lapply(seq(1,grid),function(i) frf[,,i] %*% spec.xz[,,i]))
  fpsi <- -1*grid^{-1}*deltad^{-1}*fpsi %*% (rep(1,grid) %x% diag(r))
  opt.val <- do.call(cbind,lapply(seq(1,grid),function(i) frf[,,i] %*% 
                      spec.xx[,,i] %*% Conj(t(frf[,,i]))))
  opt.val <- grid^{-1}*deltad^{-2}*opt.val %*% (rep(1,grid) %x% diag(N))
  fzzmat <- grid^{-1}*matrix(spec.zz,nrow=r) %*% (lambda.ft^{0} %x% diag(r))
  for(k in 0:(q-1))
  {
    fpsi.new <- do.call(cbind,lapply(seq(1,grid),function(i) frf[,,i] %*% spec.xx[,,i]))
    fpsi.new <- grid^{-1}*deltad^{-2}*fpsi.new %*% (lambda.ft^{-k} %x% diag(N))
    fpsi <- cbind(fpsi,fpsi.new)
    fxxmat.new <- grid^{-1}*deltad^{-2}*matrix(spec.xx,nrow=N) %*% (lambda.ft^{-k} %x% diag(N))
    if(k==0) { 
      fxxmat <- fxxmat.new 
      fxxzero <- fxxmat.new
    } else {
      if(k==1) {
        fxxmat <- cbind(fxxmat,fxxmat.new)
        fxxmat <- rbind(fxxmat,cbind(t(fxxmat.new),fxxzero))
      } else {
        sidexx.mat <- fxxmat[1:(dim(fxxmat)[2]-N),(dim(fxxmat)[2]+1-N):dim(fxxmat)[2],drop=FALSE]
        fxxmat <- cbind(fxxmat,rbind(fxxmat.new,sidexx.mat))
        fxxmat <- rbind(fxxmat,cbind(t(fxxmat.new),t(sidexx.mat),fxxzero))
      }
    }
    fxzmat.new <- -1*grid^{-1}*deltad^{-1}*matrix(spec.xz,nrow=N) %*% (lambda.ft^{k} %x% diag(r))
    if(k==0) { 
      fxzmat <- fxzmat.new 
    } else { fxzmat <- rbind(fxzmat,fxzmat.new) }
  }
  fmat <- cbind(fzzmat,t(fxzmat))
  fmat <- rbind(fmat,cbind(fxzmat,fxxmat))
  fpsi <- Re(fpsi)
  fmat <- Re(fmat)
  
  opt <- solve(fmat) %*% t(fpsi)
  opt.val <- Re(opt.val) - fpsi %*% solve(fmat) %*% t(fpsi)
  alpha <- opt[seq(1,r),,drop=FALSE]
  opt <- opt[-seq(1,r),]
  opt.array <- array(t(opt),c(N,N,q))
  
  return(list(alpha,opt.array,opt.val)) 
}



