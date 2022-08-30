# Simulates a number of observations from 
# two correlated binomials

# Outputs a matrix of observations
Generate2SnpBernoulli <- function(p1, p2, rho, nObs){
  # p1 minor allele frequency of snp 1
  # p2 minor allele frequency of snp 2
  # rho correlation between SNP1 and SNP2
  # nObs number of observations to generate
  
  
  # Compute the four probabilities for the joint distribution.
  a0 <-  rho * sqrt(p1*p2*(1-p1)*(1-p2)) + (1-p1)*(1-p2)
  prob <-
    c(
      `(0,0)` = a0,
      `(1,0)` = 1 - p2 - a0,
      `(0,1)` = 1 - p1 - a0,
      `(1,1)` = a0 + p1 + p2 - 1
    )
  
  # Checks for valid probabilities, skips over if invalid
  # Negative probabilities due to choice of p1, p2 and 
  # correlation coefficient can cause this
  if (min(prob) < 0) {
    return(as.matrix(cbind(NA, NA)))
  }
  
  nbin <- 2
  
  # Sample integers with specific weights
  u <- sample.int(4, nBinObs * nbin, replace=TRUE, prob=prob)
  
  # Split into bernoulli representation
  y <- floor((u-1)/2)
  x <- 1 - u %% 2
  
  # Sum up "columns" to get binomial
  xx <- colSums(matrix(x, nrow=nbin)) # Sum in groups of `n`
  yy <- colSums(matrix(y, nrow=nbin)) # Sum in groups of `n`

  return(as.matrix(cbind(xx,yy)))  
}

# Generate2SnpBernoulli(p1, p2, 1, nBinObs)
