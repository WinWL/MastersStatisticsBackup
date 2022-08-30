# Generates integer given the probability for each integer
# nObs = number of observations to generate
# pGeno = vector of probability weights
genBag <- function(nObs, pGeno){
  if (!all.equal(sum(pGeno), 1)) print("Given probabilities don't add to 1 ")
  
  # Generate and return values
  genBagVals <- sample.int(length(pGeno), nObs , replace=TRUE, prob=pGeno)
  return(genBagVals)
}

# genUnifVal, pmfMat

# cc <- sample.int(length(pGeno), nObs , replace=TRUE, prob=pGeno)
# pGeno
# table(cc)/length(cc)
# 
# # sample.int with replacement might be far better.
# 
# # Major Frequency from minor allele frequency (the input)
# pA <- 1 - 0.1
# pB <- 1 - 0.15
# D  <- 0
# 
# # 2x2 table of frequencies (Haplotype)
# pAB <- pA*pB + D
# paB <- (1-pA)*pB - D
# pAb <- pA*(1-pB) - D
# pab <- (1-pA)*(1-pB) + D
# 
# # 3x3 table of frequencies (Genotype)
# pAABB <- pAB^2
# pAaBB <- 2*pAB*paB
# paaBB <- paB^2
# pAABb <- 2*pAB*pAb
# pAaBb <- 2*pAB*pab + 2*pAb*paB
# paaBb <- 2*paB*pab
# pAAbb <- pAb^2
# pAabb <- 2*pAb*pab
# paabb <- pab^2
# 
# # Assemble into single vector
# pGeno <- c(
#   pAABB, pAaBB, paaBB,
#   pAABb, pAaBb, paaBb,
#   pAAbb, pAabb, paabb
# )
# 
# 
# genMultinomVals <- genMultinomial(1e7, pGeno)
# zz <- genMultinomVals %>% table
# ((zz/1e7 - pGeno)/pGeno*100) %>% abs %>% signif(2)
# 
# zz/1e7 - pGeno
# pGeno
# zz
