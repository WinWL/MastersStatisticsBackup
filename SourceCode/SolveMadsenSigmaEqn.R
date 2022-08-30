# Numerically solves the what the entry in the correlation matrix
# of the multivariate normal should be to get the desired correlation
# after transformation to correlated binomial. 

# Uses optimize to solve the equation. 
# Equation is from Madsen and Birkes (2013)
solveSigmaEqnForNormToBin <- function(corDesire, sigmaI2, sigmaJ2, 
                                      uI, uJ, nI, nJ, pI, pJ){
  # corDesire is the correlation desired between two binomial r.v.
  # sigmaI2, uI, nI, pI are specified parameters for the first binomial r.v.
  # sigmaJ2, uJ, nJ, pJ are specified parameters for the second binomial r.v.
  require(mvtnorm)
  
  # Initialize the values we will be using
  sigmaI <- sqrt(sigmaI2)
  sigmaJ <- sqrt(sigmaJ2)
  
  # Support of distributions. 2 Binomials
  rVals <- c(0,1,2,0,1,2,0,1,2)
  sVals <- c(0,0,0,1,1,1,2,2,2)
  
  # Cdf of R and S values. 
  FiRVec <- pbinom(rVals, nI, pI)
  FjSVec <- pbinom(sVals, nJ, pJ)
  
  # Inverse of Standard Normal of FiR and FiS
  NInvFiRVec <- qnorm(FiRVec)
  NInvFjSVec <- qnorm(FjSVec)
  
  # Function we will be optimizing to find delta
  # delta is what the value should be in the correlation matrix for 
  # the multivariate normal
  fx.opt <- function(deltaInit){
    delta <- deltaInit[1]
    
    # Function to create bivariate normal functions with the delta for optimization
    fx.BivarNorm <- function(nFiR, nFiS) {
      sigDelta <- matrix(c(1,delta,delta,1), nrow = 2, ncol = 2)
      pmvnorm(upper = c(nFiR, nFiS), mean = c(0,0), sigma = sigDelta)
    }
    
    # Constructs vector so we can sum the bivariate normal values up
    BiNormVec <- vector(mode = "double", length = length(rVals))
    for (i in c(1:9)) {
      BiNormVec[i] <- fx.BivarNorm(NInvFiRVec[i], NInvFjSVec[i])
    }
    
    # Assembling the equation
    eqSumTermVec <- 1 - FiRVec - FjSVec + BiNormVec
    
    rightEqTerm <- 1/(sigmaI*sigmaJ) * (sum(eqSumTermVec) - uI*uJ)
    leftTerm <- corDesire
    
    # What is being optimizd
    (rightEqTerm - leftTerm)^2
  }
  
  # Delta is between -1 to 1
  opt <- optimize(fx.opt, interval = c(-1,1))
  estimatedDelta <- opt$minimum
  return(estimatedDelta)
}
