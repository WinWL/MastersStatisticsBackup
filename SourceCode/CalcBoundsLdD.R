# Calculates the Bounds on D from LD
# Bounds are given by finding the max/min values
# such that the probabilities of all SNP allele combinations
# are still between 0 and 1

getLdDRangeEqs <- function(p1, p2){
  
  # Testing values
  # p1 <- 0.01
  # p2 <- 0.015
  
  LowerD <- max(
    -(1-p1)*(1-p2),
    p1*(1-p2) - 1,
    p2*(1-p2) - 1,
    -p1*p2
  )
  
  UpperD <- min(
    1 - (1-p1)*(1-p2),
    p1*(1-p2),
    p2*(1-p2),
    1 - p1*p2
  )
  
  lrUr <- c(LowerD, UpperD)
  names(lrUr) <- c("Lower", "Upper")
  return(lrUr)
}


# getLdDRangeEqs(0.1, 0.9)
# getLdDRangeEqs(0.2, 0.9)
# getLdDRangeEqs(0.3, 0.9)
