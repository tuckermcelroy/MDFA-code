sqrtm <- function(A){
  # Function returns the square root of the matrix
  # Input:
  #   A: m x m non-negative definite matrix
  # Output:
  #   m x m matrix

  return(eigen(A)$vectors %*% diag(sqrt(eigen(A)$values)) %*% t(eigen(A)$vectors))
}
