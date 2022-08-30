# Takes a matrix and adjusts it to be positive definite
# by setting any non-positive eigen values to a small positive number

CalcAdjustMadsenInverseEigen <- function(A, tol = 1e-5){
  # A is the matrix to adjust
  # tol is the small positive to set negative eigen values to be
  
  eigens <- eigen(A)
  eigVal <- eigens$values
  eigVec <- eigens$vectors
  
  # A = Q %*% Lambda %*% Q^-1
  Q <- eigVec
  Qinv <- solve(Q)
  
  # Adjust the lamda by setting any negatives to some small
  # positive number
  eigValAdj <- ifelse(eigVal <= 0, tol, eigVal)
  lamadj <- diag(eigValAdj, 
              ncol = length(eigVal), nrow = length(eigVal))
  
  # Construct adjusted matrix
  Aadj <- Q %*% lamadj %*% Qinv
  
  # Return matrix as correlation matrix
  return(cov2cor(Aadj))
}

# matVec <- c(1.0000000, 0.8254860, 0.8693058, 
#          0.8254860, 1.0000000, 0.9983935, 
#          0.8693058, 0.9983935, 1.0000000)
# A <- matrix(matVec, 3,3,T)  
# A
# tol <- 1e-5
# CalcMadsenInverseEigenAdjust(A, tol)


