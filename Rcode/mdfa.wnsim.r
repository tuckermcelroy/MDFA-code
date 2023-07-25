mdfa.wnsim <- function(psi,ranks,T,dof)
{
  #	generates T x N white noise simulation with student t disturbances
  #		of dof degrees of freedom, with covariance matrix determined
  #		by parameters psi.  
  #  ranks is a N-vector of ones and zeroes, gives rank configuration
  
  N <- length(ranks)
  A.mat <- matrix(0,N,N)
  A.mat[lower.tri(A.mat)] <- 1
  D.dim <- length(ranks[ranks==1])
  L.dim <- sum(A.mat[,seq(1,N)[ranks==1],drop=FALSE])
  L.psi <- NULL
  if(L.dim > 0) L.psi <- psi[1:L.dim]
  D.psi <- psi[(L.dim+1):(L.dim+D.dim)]
  L.mat <- diag(N)
  L.sub <-  L.mat[,seq(1,N)[ranks==1],drop=FALSE]
  L.sub[lower.tri(L.mat)[,seq(1,N)[ranks==1]]] <- L.psi
  Sigma <- L.sub %*% diag(exp(D.psi),nrow=D.dim) %*% t(L.sub)
  if(dof == Inf) { eps <- matrix(rnorm(D.dim*T),nrow=D.dim) } else {
    eps <- matrix(rt(D.dim*T,df=dof),nrow=D.dim) }
  eps <- L.sub %*% diag(exp(D.psi/2),nrow=D.dim) %*% eps
  
  return(list(t(eps),Sigma)) 
}