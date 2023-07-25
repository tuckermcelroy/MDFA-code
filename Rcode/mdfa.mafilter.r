mdfa.mafilter <- function(psi.array,acf.array,R,Q)
{

	#######################################################
	#
	#	mdfa.mafilter by Tucker McElroy
	#
	#	computes optimal concurrent moving average filter
	#		as approx of given target filter, using MDFA method,
	#		based upon the moving average filter class of length q
	#	inputs:
	#		psi.array: the target filter at lags r+q-2 through 1-r,
	#			as an N x N x (2r+q-2) array object
	#		acf.array: given acf at lags 0 through r-1,
	#			as an N x N x r array object
	#		R is array N x q x N x M of M constraints on the moving 
	#			average filter solution.
	#		Q is array N x q x N of constraints, theta = R phi + Q.
	#	outputs:
	#		opt.array is array N x N x q of filter coefficients
	#		opt.val is N x N matrix corresponding to minimal MSE
	#
	##############################################################

	N <- dim(psi.array)[2]
	r <- dim(acf.array)[3]
	q <- dim(R)[2]
	M <- dim(R)[4]
	R.mat <- matrix(R,nrow=N*q,ncol=N*M)
	Q.mat <- matrix(Q,ncol=N)

	acf.col <- rbind(t(matrix(aperm(acf.array[,,r:1],c(2,1,3)),nrow=N)),
		t(matrix(acf.array[,,2:r],nrow=N)))
	acf.newcol <- rbind(acf.col,matrix(0,N*(q-1),N))
	acf.toep <- acf.newcol
	for(j in 2:q)
	{	
		acf.newcol <- rbind(0*diag(N),acf.newcol[1:(dim(acf.newcol)[1]-N),,drop=FALSE])
		acf.toep <- cbind(acf.toep,acf.newcol)
	}
  	psiacf <- matrix(psi.array[,,(2*r+q-2):1],nrow=N) %*% acf.toep	 
	acf.toep <- acf.toep[((r-1)*N+1):((r+q-1)*N),]

	opt <- solve(t(R.mat) %*% acf.toep %*% R.mat) %*% t(R.mat) %*% (t(psiacf) - acf.toep %*% Q.mat)
	opt <- R.mat %*% opt + Q.mat

	opt.val <- 0*diag(N)
	for(j in 1:(2*r+q-2))
	{
		for(k in 1:(2*r+q-2))
		{
			ind <- j-k
			if(abs(ind) >= r) { cov.mat <- 0*diag(N) } else {
				if(ind >= 0) { cov.mat <- acf.array[,,ind+1] } else {
					cov.mat <- t(acf.array[,,-ind+1]) } }
			opt.val <- opt.val + psi.array[,,j] %*% cov.mat %*% t(psi.array[,,k])
		}
	}
	opt.val <- opt.val - t(opt) %*% acf.toep %*% opt
	opt.array <- array(t(opt),c(N,N,q))

	return(list(opt.array,opt.val)) 
}



