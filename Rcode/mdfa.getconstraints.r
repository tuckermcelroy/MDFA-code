mdfa.getconstraints <- function(frf,sigfreqs,noisefreqs,coint,q)
{

	#######################################################
	#
	#	mdfa.getconstraints by Tucker McElroy
	#
	#	computes filter constraints needed for nonstationary series,
	#		based on sigfreqs and noisefreqs.
	#	inputs:
	#		sigfreqs lists frequencies x in [-1,1] such that lambda=pi*x, 
	#			with repeats for double roots, 
	#			such that exp(-i*lambda) is a root
	#			of the signal differencing polynomial
	#		noisefreqs lists frequencies x in [-1,1] such that lambda=pi*x,  
	#			with repeats for double roots, 
	#			such that exp(-i*lambda) is a root
	#			of the noise differencing polynomial
  #   coint is N x N matrix of co-integrating row vectors, 
  #     assumed to be the same for all signal frequencies.
  #     set to zero if no co-integration constraints are assumed
	#		frf is array N x N x Grid of complex entries, the target
	#			frequency response function Psi(e^{-i lambda})
	#			for lambda given by Grid number of Fourier frequencies 
  #   q: desired order of MA filter
	#	outputs:
	#		R is array N x q x N x M of q-M constraints on the moving 
	#			average filter solution.
	#		Q is array N x q x N of constraints, theta = R phi + Q.
	#
	##############################################################

  ceps2wold <- function(ceps,q)
  {
    m <- length(ceps)
    if(q > m) { ceps <- c(ceps,rep(0,q-m)) }
    wold <- 1
    wolds <- wold
    for(j in 1:q)
    {
      wold <- sum(seq(1,j)*ceps[1:j]*wolds[j:1])/j
      wolds <- c(wolds,wold)
    }
    return(wolds)
  }
  
  roots2ceps <- function(roots,m)
  {
    p <- length(roots)
    ceps <- rep(0,m)
    for(k in 1:m)
    {	
      ceps[k] <- -1*sum(roots^(-k))/k
    }
    return(ceps)
  }
  
  N <- dim(frf)[1]
	grid <- dim(frf)[3]
	m <- floor(grid/2)

	sig.lambdas <- unique(sigfreqs)
	sig.mults <- NULL
	if(length(sigfreqs)>0) 
	{
	  for(j in 1:length(sig.lambdas))
	  {
		  sig.mults <- c(sig.mults,sum(sigfreqs == sig.lambdas[j]))
	  } 
	}

	delta.noise <- 1
	noise.lambdas <- unique(noisefreqs)
	noise.mults <- NULL
	if(length(noisefreqs)>0) 
	{
	  for(j in 1:length(noise.lambdas))
	  {
		  noise.mults <- c(noise.mults,sum(noisefreqs == noise.lambdas[j]))
	  } 
		noise.ceps <- Re(roots2ceps(exp(-1i*pi*noisefreqs),length(noisefreqs)))
	  delta.noise <- ceps2wold(noise.ceps,length(noisefreqs))
	}
	
	sig.mat <- NULL
	sig.vec <- NULL
	if(length(sig.lambdas) > 0) {
	for(k in 1:length(sig.lambdas))
	{
		j.star <- floor(sig.lambdas[k]*grid/2) + m+1
		if(j.star==(grid+1)) j.star <- 1
		zeta <- exp(-1i*sig.lambdas[k]*pi)
		sig.mat.new <- matrix(zeta^(seq(0,q-1)),ncol=q)
		sig.vec.new <- frf[,,j.star] - coint * sum(delta.noise * zeta^seq(0,length(noisefreqs)))
		sig.mat <- rbind(sig.mat,Re(sig.mat.new))
		sig.vec <- rbind(sig.vec,t(Re(sig.vec.new)))
		if( (sig.lambdas[k] != 0) && (sig.lambdas[k] != 1) )
		{ 
			sig.mat <- rbind(sig.mat,Im(sig.mat.new))
			sig.vec <- rbind(sig.vec,t(Im(sig.vec.new)))
		}
		if(sig.mults[k]==2)
		{
			sig.mat.new <- seq(0,q-1)*zeta^(seq(-1,q-2))
			if(j.star==grid) 
			{ 
				sig.vec.new <- 1i*exp(1i*sig.lambdas[k]*pi)*
					(frf[,,1]-frf[,,j.star])*grid/(2*pi) 
			} else 
			{ 
				sig.vec.new <- 1i*exp(1i*sig.lambdas[k]*pi)*
					(frf[,,j.star+1]-frf[,,j.star])*grid/(2*pi) 
			}
			sig.mat <- rbind(sig.mat,Re(sig.mat.new))
			sig.vec <- rbind(sig.vec,t(Re(sig.vec.new)))
			if( (sig.lambdas[k] != 0) && (sig.lambdas[k] != 1) )
			{ 
				sig.mat <- rbind(sig.mat,Im(sig.mat.new))
				sig.vec <- rbind(sig.vec,t(Im(sig.vec.new)))
			}
		}
 	} }
	
	noise.mat <- NULL
	noise.vec <- NULL
	if(length(noise.lambdas) > 0) {
	for(k in 1:length(noise.lambdas))  
	{
		j.star <- floor(noise.lambdas[k]*grid/2) + m+1
		if(j.star==(grid+1)) j.star <- 1
		zeta <- exp(-1i*noise.lambdas[k]*pi)
		noise.mat.new <- matrix(zeta^(seq(0,q-1)),ncol=q)
		noise.vec.new <- frf[,,j.star]
		noise.mat <- rbind(noise.mat,Re(noise.mat.new))
		noise.vec <- rbind(noise.vec,t(Re(noise.vec.new)))
		if( (noise.lambdas[k] != 0) && (noise.lambdas[k] != 1) )
		{ 
			noise.mat <- rbind(noise.mat,Im(noise.mat.new))
			noise.vec <- rbind(noise.vec,t(Im(noise.vec.new)))
		}
		if(noise.mults[k]==2)
		{
			noise.mat.new <- seq(0,q-1)*zeta^(seq(-1,q-2))
			if(j.star==grid) 
			{ 
				noise.vec.new <- 1i*exp(1i*noise.lambdas[k]*pi)*
					(frf[,,1]-frf[,,j.star])*grid/(2*pi) 
			} else 
			{ 
				noise.vec.new <- 1i*exp(1i*noise.lambdas[k]*pi)*
					(frf[,,j.star+1]-frf[,,j.star])*grid/(2*pi) 
			}	
			noise.mat <- rbind(noise.mat,Re(noise.mat.new))
			noise.vec <- rbind(noise.vec,t(Re(noise.vec.new)))
			if( (noise.lambdas[k] != 0) && (noise.lambdas[k] != 1) )
			{ 
				noise.mat <- rbind(noise.mat,Im(noise.mat.new))
				noise.vec <- rbind(noise.vec,t(Im(noise.vec.new)))
			}
		}
 	} }
 
	constraint.mat <- rbind(sig.mat,noise.mat)
	constraint.vec <- rbind(sig.vec,noise.vec)

	## compute decomposition
	constraint.qr <- qr(constraint.mat)
	constraint.q <- qr.Q(constraint.qr)
	constraint.r <- qr.R(constraint.qr)
	constraint.pivot <- constraint.qr$pivot
	constraint.ipivot <- sort.list(constraint.pivot)

	M <- q - dim(constraint.r)[2] + dim(constraint.q)[2]
	R.mat <- rbind(-solve(constraint.r[,1:(dim(constraint.q)[2]),drop=FALSE],
		constraint.r[,(dim(constraint.q)[2]+1):q,drop=FALSE]),diag(q-M))
	R.mat <- R.mat[constraint.ipivot,] %x% diag(N)
	Q.mat <- rbind(solve(constraint.r[,1:(dim(constraint.q)[2]),drop=FALSE]) %*% 
		solve(constraint.q),matrix(0,q-M,M)) 
	Q.mat <- (Q.mat[constraint.ipivot,] %x% diag(N)) %*% constraint.vec
	R <- array(R.mat,c(N,q,N,q-M))
	Q <- array(Q.mat,c(N,q,N))
 
	return(list(R,Q))
}



