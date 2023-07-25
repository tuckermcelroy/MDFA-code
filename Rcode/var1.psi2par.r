var1.psi2par <- function(psi,delta)
{
	# computes VAR(1) NxN phi matrix from N^2 vector psi, with delta = \pm 1,
	#	yielding stable parametrization

N <- sqrt(length(psi))
l.mat <- diag(N)
l.mat[lower.tri(l.mat)] <- psi[1:choose(N,2)]
d.mat <- diag(exp(psi[(choose(N, 2) + 1):choose(N + 1, 2)]))
v.mat <- l.mat %*% d.mat %*% t(l.mat)
s.mat <- diag(0,N)
s.mat[lower.tri(s.mat)] <- psi[(choose(N + 1, 2) + 1):(N^2)]
s.mat <- s.mat - t(s.mat)
E.mat <- diag(c(delta, rep(1, (N - 1))))
Q.mat <- E.mat %*% (diag(N) - s.mat) %*% solve(diag(N) + s.mat)

phi.mat <- sqrtm(v.mat) %*% Q.mat %*% solve(sqrtm(diag(N) + v.mat))
#print(Mod(eigen(phi.mat)$values)) 	# check stability

return(phi.mat)
}

