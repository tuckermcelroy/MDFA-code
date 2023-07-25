var1.par2psi <- function(param,trunc)
{
	# takes VAR(1) NxN phi matrix (param) and computes 
	#	N^2 vector psi, with delta = \pm 1,
	#	given the stable parametrization

N <- dim(param)[1]
phi.mat <- param
v.mat <- diag(0,N)
for(k in 1:trunc) {
v.mat <- v.mat + (phi.mat %^% k) %*% t(phi.mat %^% k)
}
l.mat <- t(chol(v.mat))
d.mat <- diag(l.mat)
l.mat <- l.mat %*% solve(diag(d.mat))
Q.mat <- solve(sqrtm(v.mat)) %*% phi.mat %*% sqrtm(diag(N) + v.mat)
delta <- det(Q.mat)
E.mat <- diag(c(delta, rep(1, (N - 1))))
s.mat <- solve(diag(N) + E.mat %*% Q.mat) %*% (diag(N) - E.mat %*% Q.mat)
psi <- c(l.mat[2,1],2*log(d.mat),s.mat[2,1])

return(list(psi,delta))
}


