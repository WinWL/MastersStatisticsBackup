# Simulate Binomial values using
# the pmf originating from optimizing a system of linear equations

# Returns the binomial values in a matrix and with the values used in the
# calculated correlation matrix for the multivariate normal
GenerateSnpPmf3SnpConstrained <- function(p1, p2, p3, t12, t13, t23, nBinObs)
{
  # Calculate the pmf
  # p1 <- 0.1; p2 <- 0.12; p3 <- 0.15
  # t12 <- 0.1; t13 <- 0.13; t23 <- 0.15
  # nBinObs <- 100
  pmfMat <- calcPmfParam3Snp(p1, p2, p3, t12, t13, t23)
  matSnpVals <- NA
  
  if (all(pmfMat >= 0)){
    # Assemble matrix that describess what each value of the
    # multinomial corresponds to a single snp genotype value
    
    # Generate multinomial values using sample.int
    genMultinomVals <- genMultinomial(nBinObs, pmfMat)
    
    
    # Transform the multinomial values into SNP genotype values
    # by using the mapping above
    pABmarginalMap <- cbind(
      c(0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2),
      c(0,0,0,1,1,1,2,2,2,0,0,0,1,1,1,2,2,2,0,0,0,1,1,1,2,2,2),
      c(0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2,0,1,2)
    )
    
    SNP1 <- pABmarginalMap[,1][genMultinomVals]
    SNP2 <- pABmarginalMap[,2][genMultinomVals]
    SNP3 <- pABmarginalMap[,3][genMultinomVals]
    matSnpVals <- cbind(SNP1, SNP2, SNP3)
  }
  return(matSnpVals)
}

# p1 <- 0.1; p2 <- 0.12; p3 <- 0.15
# t12 <- 0.1; t13 <- 0.13; t23 <- 0.15
# nBinObs <- 1000*1000
# GenerateSnpPmf3SnpConstrained(p1, p2, p3, t12, t13, t23, nBinObs)
