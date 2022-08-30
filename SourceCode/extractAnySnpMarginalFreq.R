# Winfield Lai
# Extracts the Binomial (2, p) or SNP frequencies from a sample.

# Input: binMat: n column matrix, entries are 0,1,2
# Input: pAllele: Vector of minor allele frequencies

# Output: single row df whose columns are the sample
# and expected frequncies of binomial 0,1,2 values

# rm(list=ls())
extractAnySnpMarginalFreq <- function(binMat, pAllele) {
  require(reshape2)
  nObs <- nrow(binMat)
  
  # Snp Frequencies of generated Sample
  freqMat.sam <- apply(binMat, 2, function(binCol) {
    freq0sam <- sum(binCol == 0) / nObs
    freq1sam <- sum(binCol == 1) / nObs
    freq2sam <- sum(binCol == 2) / nObs
    
    return(t(c(freq0sam, freq1sam, freq2sam)))
  })

  # Expected SNP frequencies
  freq0exp <- sapply(pAllele, function(p)
    (1 - p) ^ 2)
  freq1exp <- sapply(pAllele, function(p)
    2 * p * (1 - p))
  freq2exp <- sapply(pAllele, function(p)
    p ^ 2)
  
  # Combine into single matrix
  freqMat <- rbind(freqMat.sam, freq0exp, freq1exp, freq2exp)
  
  # Set up the row and column names so we can get a data frame with a column of frequencies and
  colnames(freqMat) <- paste(rep("SNP", ncol(binMat)),
                             c(1:ncol(binMat)),
                             sep = "")
  rownames(freqMat) <- c("Freq0Samp",
                         "Freq1Samp",
                         "Freq2Samp",
                         "Freq0Expe",
                         "Freq1Expe",
                         "Freq2Expe")
  meltFreq <- melt(freqMat)
  
  # Convert to single row data frame
  singleRowFreq <- data.frame(t(meltFreq$value))
  colnames(singleRowFreq) <- paste(meltFreq[, 2],
                                   meltFreq[, 1],
                                   sep = "")
  
  return(singleRowFreq)
}

# p1A <- 0.1; p2A <- 0.15; p3A <- 0.2
# pAllele <- c(p1A, p2A, p3A)
# binMat <- cbind(rbinom(n = 1e7, size = 2, prob = p1A),
#                 rbinom(n = 1e7, size = 2, prob = p2A),
#                 rbinom(n = 1e7, size = 2, prob = p3A))
# df <- extractAnySnpMarginalFreq(binMat, pAllele)
# df %>% t
# df
# library(dplyr)
