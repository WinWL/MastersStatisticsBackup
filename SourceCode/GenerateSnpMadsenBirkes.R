# Simulate observations from correlated binomials using 
# the Madsen Birkes method

# Returns the generated binomial values and values used in the
# calculated correlation matrix in a matrix
GenerateSnpMadsenBirkes <- function(pAlleleMinor, targetBinCorMat, 
                                    nBinObs, useEigAdjust = F){
  # pAlleleMinor vector of minor allele frequencies
  # targetBinCorMat Matrix of desired correlations
  # nBinObs Number of observations to generate
  
  require(gdata)
  require(matrixcalc)
  require(MASS)
  nBin <- rep(2, length(pAlleleMinor))

  # Calculate correlation matrix for multivariate normal
  sigmaZ <- CalculateSigmaZ(nBin, pAlleleMinor, targetBinCorMat)
  
  # Ensures the correlation matrix is positive definite if desired
  if (useEigAdjust) {
    sigmaZ <- CalcAdjustMadsenInverseEigen(sigmaZ)
  }
  
  if (is.symmetric.matrix(sigmaZ)){
    if (!is.positive.definite(sigmaZ)){
      llen <- length(pAlleleMinor) + ncol(sigmaZ)
      retVal <- matrix(NA, ncol = llen)
      return(retVal)
    }
  }

  normMeans <- nBin * 0 # mean is just 0
  
  # Simulate from multivariate normal
  generatedNormalObs <- mvrnorm(nBinObs, normMeans, sigmaZ)
  
  # Convert values from multivariate normal to binomial 
  transformedNormToBinMat <- as.matrix(
    TransformNormToBin(
      pAlleleMinor = pAlleleMinor,
      normalGeneratedMat = generatedNormalObs)
  )
  
  # Return value
  retVal <- cbind(transformedNormToBinMat,
                  matrix(upperTriangle(sigmaZ), byrow = T, 
                         ncol = length(upperTriangle(sigmaZ)), 
                         nrow = nBinObs))
  return(as.matrix(retVal))
}
# GenerateSnpMadsenBirkes( pAlleleMinor, targetBinCorMat, nBinObs)



# cbind(transformedNormToBinMat,
#       matrix(upperTriangle(sigmaZ), byrow = T, ncol = ncol(sigmaZ), nrow = nBinObs))
# 
# matrix(upperTriangle(sigmaZ), byrow = T, ncol = ncol(sigmaZ), nrow = nBinObs)

##############
# 2 Snps
# p1Allele <- 0.1
# p2Allele <- 0.15
# targetBinCor <- 0.3
# nBinObs <- 100
# # 
# pAlleleMinor <- c(p1Allele, p2Allele)
# targetBinCorMat <- matrix(c(1, targetBinCor, targetBinCor, 1),
#                           nrow = 2,
#                           ncol = 2)
# # 
# GenerateSnpMadsenBirkes( pAlleleMinor, targetBinCorMat, nBinObs)

##############
# 3 Snps

# p1Allele <- 0.1
# p2Allele <- 0.1
# p3Allele <- 0.1
# tc12 <- 0.1
# tc13 <- 0.2
# tc23 <- 0.3
# nBinObs <- 1000000
# 
# pAlleleMinor <- c(p1Allele, p2Allele, p3Allele)
# targetBinCorMat <- matrix(c(1, tc12, tc13,
#                             tc12, 1, tc23,
#                             tc13, tc23, 1),
#                           nrow = 3, ncol = 3, byrow = T)
# 
# 
# # Simulate Using Madsen and Birkes method
# GenerateSnpMadsenBirkes(pAlleleMinor, targetBinCorMat, nBinObs)

# p <- 0.1
# t <- 0.15
# nBinObs <- 1000
# nSnp <- 4
# 
# pAlleleMinor <- rep(p, nSnp)
# targetBinCorMat <- matrix(c(t),
#                           nrow = nSnp, ncol = nSnp, byrow = T)
# diag(targetBinCorMat) <- 1
# GenerateSnpMadsenBirkes(pAlleleMinor, targetBinCorMat, nBinObs)
