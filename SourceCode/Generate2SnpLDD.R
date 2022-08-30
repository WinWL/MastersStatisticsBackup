# Generates Binomial SNP values using the
# Method from LD measure D


Generate2SnpLDD <- function(p1, p2, D, nObs){
  # pa is the minor allele frequency of snp 1
  # pb is the minor allele frequency of snp 2
  # D is the target D
  
  # Major Frequency from minor allele frequency (the input)
  pA <- 1 - p1
  pB <- 1 - p2

  # 2x2 table of frequencies (Haplotype)
  pAB <- pA*pB + D
  paB <- (1-pA)*pB - D
  pAb <- pA*(1-pB) - D
  pab <- (1-pA)*(1-pB) + D
  
  # 3x3 table of frequencies (Genotype)
  pAABB <- pAB^2
  pAaBB <- 2*pAB*paB
  paaBB <- paB^2
  pAABb <- 2*pAB*pAb
  pAaBb <- 2*pAB*pab + 2*pAb*paB
  paaBb <- 2*paB*pab
  pAAbb <- pAb^2
  pAabb <- 2*pAb*pab
  paabb <- pab^2
  
  # Assemble into single vector
  # Column is AA, Aa, aa
  # Row is BB, Bb, bb
  # Numbered 1-9 in multinomial in order
  pGeno <- c(
    pAABB, pAaBB, paaBB,
    pAABb, pAaBb, paaBb,
    pAAbb, pAabb, paabb
  )
  

  # Check if any of the probabilities are negative
  if (Reduce('|',pGeno < 0)) {
    return(cbind(NA, NA))
  }
  
  pGenoMat <- matrix(pGeno, ncol = 3, nrow = 3, byrow = T)
  
  
  # Assemble matrix that describess what each value of the
  # multinomial corresponds to a single snp genotype value
  pABmarginalMap <- cbind(
    c(0,1,2,0,1,2,0,1,2), 
    c(0,0,0,1,1,1,2,2,2) 
  )
  
  # Generate multinomial values from weighted sampling
  genMultinomVals <- genMultinomial(nObs, pGeno)
  
  # Transform the multinomial values into SNP genotype values
  # by using the mapping above
  matSnpVals <- apply(as.matrix(genMultinomVals), 1, function(multinomVal) {
    
    # For each SNP, we apply the mapping to genotype
    retval <- apply(pABmarginalMap, 2, function(snpMap) snpMap[multinomVal] )
    
    return(retval)
  }) %>% t()
  
  return(matSnpVals)
}

# p1 <- 0.1
# p2 <- 0.1
# D  <- 0.5
# nObs <- 1000
# 
# Generate2SnpLDD(p1, p2, D, nObs)

