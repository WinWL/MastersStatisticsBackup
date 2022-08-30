# Input: binMat: n column matrix, entries are 0,1,2

# Output: single row df whose columns are the sample 
# and expected frequencies of multivariate binomial 0,1,2 values

# e.g. in the 2 SNP case
# (0,0), (0,1), (0,2), 
# (1,0), (1,1), (1,2)
# (2,0), (2,1), (2,2)

extractSnpMultivarDfVec <- function(binMat) {
  require(plyr)
  
  # Multivariate Frequencies
  # Add dummy frequencies so the count function will incorporate all frequencies even
  # if they are 0 observations
  b <- na.omit(as.data.frame(addDummyFreq(binMat)))
  nSnp <- ncol(binMat)
  
  # Subtraction to get rid of the dummy values
  multiSamFreq <- (plyr::count(b)[, nSnp + 1] - 1) / nrow(binMat)
  
  
  names(multiSamFreq)  <- paste("Freq",
                                substring(namingForBinomialPmf(nSnp), 2),
                                "Samp",
                                sep = "")
  
  return(multiSamFreq)
}

# binMat <- matBinVal
# 
# pAllele <- c(0.1, 0.2, 0.3)
# binMat <- cbind(
#   rbinom(n = 1e5, size = 2, prob = pAllele[1]),
#   rbinom(n = 1e5, size = 2, prob = pAllele[2]),
#   rbinom(n = 1e5, size = 2, prob = pAllele[3])
# )
# 
# df <- extractSnpMultivarDfVec(binMat)
# df
# library(dplyr)
# df %>% t
# df
# sum(df)
# z1 <- binMat[,1] %>% head(1e8)
# z2 <- binMat[,2] %>% head(1e8)
# 
# Paste all columns together with no separation
# zp <- paste(z1,z2,sep="")
# 
# Convert back to integer and use table to count
# zp %>% as.integer() %>% table()/1e8
# Rename table headers using the naming system for the pmf
# Easy way to get the multivariate pmf for more than 1 snp

