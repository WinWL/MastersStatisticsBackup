# Simulate Binomial values using
# the pmf originating from optimizing a system of linear equations

# Returns the binomial values in a matrix and with the values used in the
# calculated correlation matrix for the multivariate normal
GenerateSnpPmf2SnpConstrained <- function(pCol, pRow, targetCor, nBinObs)
  {
  # Calculate the pmf
  pmfRnd <- calcOptPmf2SnpBinomial(pCol, pRow, targetCor)
  
  # pmfRnd
  # Adjusts the pmf so no small negatives (if they exist)
  # negtol <- -0.05
  # pmfRndAdjpos <- ifelse(pmfRnd < 0 & pmfRnd > negtol,
  #                        0, pmfRnd)
  # pmfRndAdjposOne <- pmfRndAdjpos/sum(pmfRndAdjpos)
  
  # Convert pmf to a matrix for 2 SNPs
  pmfMat <- matrix(pmfRnd, ncol = 3, nrow = 3, byrow = T)
  matSnpVals <- NA
  
  if (all(pmfMat >= 0)){
    # Assemble matrix that describess what each value of the
    # multinomial corresponds to a single snp genotype value
    
    # Generate multinomial values using sample.int
    genMultinomVals <- genMultinomial(nBinObs, pmfMat)
    
    
    # Transform the multinomial values into SNP genotype values
    # by using the mapping above
    pABmarginalMap <- cbind(
      c(0,1,2,0,1,2,0,1,2), 
      c(0,0,0,1,1,1,2,2,2) 
    )
    SNP1 <- pABmarginalMap[,1][genMultinomVals]
    SNP2 <- pABmarginalMap[,2][genMultinomVals]
    matSnpVals <- cbind(SNP1, SNP2)
  }
  return(matSnpVals)
}

# pCol <- 0.1
# pRow <- 0.15
# nBinObs <- 1000
# targetCor <- 0.1

# GenerateSnpPmf2SnpConstrained(pCol, pRow, targetCor, nBinObs)
