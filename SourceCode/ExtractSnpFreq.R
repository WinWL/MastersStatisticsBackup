# Winfield Lai
# Extracts the Binomial (2, p) or SNP frequencies from a sample. 

# Input: binMat: 2 column matrix, entries are 0,1,2
# Input: p1A: probability of minor allele
# Input: p2A: probability of minor allele

# Output: single row df whose columns are the sample 
# and expected frequncies of binomial 0,1,2 values 

# rm(list=ls())
extractSnpFreqDfVec <- function(binMat, p1A, p2A) {
    require(reshape2)
  
  # Snp Frequencies of generated Sample
  freqMat.sam <- apply(binMat, 2, function(binCol) {
    freq0sam <- sum(binCol == 0)/length(binCol)
    freq1sam <- sum(binCol == 1)/length(binCol)
    freq2sam <- sum(binCol == 2)/length(binCol)
    
    return(t(c(freq0sam, freq1sam, freq2sam)))
  })
  
  # Expected SNP frequencies
  freq0exp <- cbind((1-p1A)^2, (1-p2A)^2)
  freq1exp <- cbind(2*p1A*(1-p1A), 2*p2A*(1-p2A))
  freq2exp <- cbind(p1A^2, p2A^2)
  
  freqMat <- rbind(freqMat.sam, freq0exp, freq1exp, freq2exp)
  
  # Set up the row and column names so we can get a data frame with a column of frequencies and
  colnames(freqMat) <- paste(
    rep("SNP", ncol(binMat)), 
    c(1:ncol(binMat)), 
    sep = "")
  rownames(freqMat) <- c("Freq0Samp", "Freq1Samp", "Freq2Samp",
                         "Freq0Expe", "Freq1Expe", "Freq2Expe")
  meltFreq <- melt(freqMat)
  
  # Convert to single row data frame
  singleRowFreq <- data.frame(t(meltFreq$value))
  colnames(singleRowFreq) <- paste(
    meltFreq[,2],
    meltFreq[,1],
    sep = ""
  )
  
  return(singleRowFreq)
}

# p1A <- 0.1
# p2A <- 0.15
# bin1 <- rbinom(n = 1e7, size = 2, prob = p1A)
# bin2 <- rbinom(n = 1e7, size = 2, prob = p2A)
# binMat <- cbind(bin1, bin2)
# df <- extractSnpFreqDfVec(binMat, p1A, p2A)
# df %>% t
# df
# library(dplyr)

