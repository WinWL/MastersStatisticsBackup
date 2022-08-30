# Gets the name for each entry in the PMF of correlated binomials
namingForBinomialPmf <- function(numSnp){
  # numSNP = number of SNPs
  # Naming convention follows a specific pattern
  # they are the ordered (lowest to highest) of what 
  # value (1,2,3) each snp (indicated by index) takes on 
  # For example, if snp1 is 1, snp2 is 3 and snp3 is 2, 
  # then the corresponding name is 132.

  # Construct the names from right to left.
  # The pattern of 1,2,3 is a repetitive one.
  nmedoub <- rep(0, 3^numSnp)
  for (snpInd in c(1:numSnp)) {
    # snpInd <- 1
    inner <- 3 ^ (snpInd - 1)
    outer <- 3 ^ (numSnp - snpInd)
    eexp <- snpInd - 1
    
    toAdd <- (10^eexp)*rep(
      c(rep(1, inner),
        rep(2, inner),
        rep(3, inner)),
      outer
    )
    toAdd
    nmedoub <- nmedoub + toAdd
  }
  
  nmestri <- paste("x", nmedoub, sep = "")
  nmestri
}

# Gets the name for each entry in the PMF of correlated Bernoullis
namingForBernoullipmf <- function(numSnp){
  # Naming convention follows a specific pattern
  # they are the ordered (lowest to highest) of what value (1,2) 
  # each snp (indicated by index) takes on 
  # For example, if snp1 is 1, snp2 is 1 and snp3 is 2, 
  # then the corresponding name is 112.

  # Construct the names from right to left.
  # The pattern of 1,2 is a repetitive one. 
  nmedoub <- rep(0, 2^numSnp)
  for (snpInd in c(1:numSnp)) {
    # snpInd <- 1
    inner <- 2 ^ (snpInd - 1)
    outer <- 2 ^ (numSnp - snpInd)
    eexp <- snpInd - 1
    
    toAdd <- (10^eexp)*rep(
      c(rep(1, inner),
        rep(2, inner)),
      outer
    )
    toAdd
    nmedoub <- nmedoub + toAdd
  }
  
  nmestri <- paste("x", nmedoub, sep = "")
  nmestri
}
