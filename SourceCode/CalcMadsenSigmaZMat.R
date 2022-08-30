# Calculates the correlation matrix used for the multivariate normal 
# given the target correlation matrix for the correlated binomial random vector

CalculateSigmaZ <- function(nBin, pAlleleMinor, corDesire) {
  # nBin = Should be a vector of just 2's. The number of "trials"
  # pAlleleMinor = vector of minor allele frequencies
  # corDesire = correlation matrix for the target r.v.
  
  # Calculate correlation matrix for normals
  corLen <- length(corDesire[,1])
  sigmaZ <- diag(length(corDesire[, 1]))
  
  # Note: Can speed up by running loop in parallel
  # Calculates the correlation matrix
  for (i in c(1:corLen)) {
    for (j in c(1:corLen)) {
      
      # Initialize values needed for calculating delta
      corD <- corDesire[i, j]
      pI <- pAlleleMinor[i]
      pJ <- pAlleleMinor[j]
      nI <- nBin[i]
      nJ <- nBin[j]
      sigmaI2 <- pI * nI * (1 - pI)
      sigmaJ2 <- pJ * nJ * (1 - pJ)
      uI <- nI * pI
      uJ <- nJ * pJ

      # Calculate the matrix entry by entry
      # The matrix is symmetric and main diagonal is all 1's
      if (i == j) {
        # diagonal elements of correlation matrix
        sigmaZ[i,j] <- 1
      } else if (i < j) {
        # Compute the entry in the upper triangle part 
        sigmaZ[i,j] <- solveSigmaEqnForNormToBin(corD, sigmaI2, sigmaJ2, 
                                                 uI, uJ, nI, nJ, pI, pJ)
      } else if (i > j) {
        # Copy the upper triangle part for the lower triangle 
        sigmaZ[i,j] <- sigmaZ[j,i]
      }
         
    }
  }
  
  return(sigmaZ)
}

# nObs <- 1
# nBin <- c(2,2)
# pAlleleMinor
# corDesire
# CalculateSigmaZ(nObs, nBin, pAlleleMinor, corDesire)