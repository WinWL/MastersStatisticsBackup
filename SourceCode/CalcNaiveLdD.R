# Calculates D and D' based on the total number of individuals
# to get p00 (person has both major alles on chromosome). 
# Note that this approach does not give the same values as using 
# mle to find p00 as seen when using the R genetics package

# Input: Vector of p1 and p2 minor allele frequencies
# Input: Vector of 2Snp frequencies

# Output: matrix of D and D'
calcNaiveLdD <- function(p1MinorAllele, p2MinorAllele, 
                         Freq00Samp, Freq01Samp, Freq02Samp,
                         Freq10Samp, Freq11Samp, Freq12Samp,
                         Freq20Samp, Freq21Samp, Freq22Samp){
  
  
  
  # Counting how many of each different combination of Aa and Bb (Freq00,..etc)
  # contribute to seeing AB, Ab, aB, bb
  p00 <- (1*Freq00Samp + 0.5*Freq10Samp + 0.5*Freq01Samp + 0.25*Freq11Samp)
  p02 <- (1*Freq02Samp + 0.5*Freq12Samp + 0.5*Freq01Samp + 0.25*Freq11Samp)
  p20 <- (1*Freq20Samp + 0.5*Freq10Samp + 0.5*Freq21Samp + 0.25*Freq11Samp)
  p22 <- (1*Freq22Samp + 0.5*Freq12Samp + 0.5*Freq21Samp + 0.25*Freq11Samp)
  
  # For ease of notation
  pA <- (1-p1MinorAllele)
  pB <- (1-p2MinorAllele)
  
  # Calculate D
  Dsam <- p00 - pA*pB
  
  # Calculate D' 
  DmaxBel0 <- pmax(-pA*pB, -(1-pA)*(1-pB))
  DmaxAbv0 <- pmin(pA*(1-pB), (1-pA)*pB)
  DPrimesam <- ifelse(Dsam > 0, Dsam/DmaxAbv0, Dsam/DmaxBel0)
  
  # return value
  retval <- cbind(Dsam, DPrimesam)
  colnames(retval) <- c("Dsam", "DPrimeSam")
  return(retval)
}
