# Input: binMat: n column matrix, entries are 0,1,2

# Output: single row df whose columns are the sample 
# and expected frequencies of marginal binomial 0,1,2 values
extractSnpMarginalSamDfVec <- function(binMat) {
  require(plyr)
  
  # Marginal Frequencies
  # Add dummy frequencies so the count function will incorporate all frequencies even
  # if they are 0 observations
  b <- na.omit(as.data.frame(addDummyFreq(binMat)))
  nSnp <- ncol(binMat)
  nObs <- nrow(binMat)
  
  marginalFreqsMat <- apply(b, 2, function(x) {
    # Subtraction to get rid of the dummy values
    (plyr::count(x)[, 2] - 3^(nSnp-1)) / nObs
  })
  
  marginalFreqsVec <- as.vector(marginalFreqsMat)
  
  freqName <- rep( c("Freq0Samp", "Freq1Samp", "Freq2Samp"),nSnp)
  snpName <- as.vector(sapply(c(1:nSnp), function(x) rep(x,3)))
  snpFreqName <- paste("SNP", snpName, freqName, sep="")
  names(marginalFreqsVec) <- snpFreqName
  return(marginalFreqsVec)
}

# Returns a matrix of values of the possible correlated binomial values
# This function is used to ensure there are no 0 observations for the count function. 
# This fixes the issue of the count function being unable to be used for
# frequency as it cannot count values with 0 observations
addDummyFreq <- function(binMat){
  nSnp <- ncol(binMat)
  dummyMat <- vector("double", 0)
  
  for (snpInd in c(1:nSnp)) {
    # snpInd <- 3
    inner <- 3 ^ (snpInd - 1)
    outer <- 3 ^ (nSnp - snpInd)
    
    toAdd <- rep(
      c(rep(0, inner),
        rep(1, inner),
        rep(2, inner)),
      outer
    )
    
    dummyMat <- cbind(toAdd, dummyMat)
  }
  colnames(dummyMat) <- paste("Snp", c(1:nSnp), "")
  colnames(binMat) <- paste("Snp", c(1:nSnp), "")
  
  matWithDummyFreq <- rbind(binMat, dummyMat)
  return(matWithDummyFreq)
}

# Extracts the sample minor allele frequencies
extractSamPAllele <- function(binMat){
  samPAllele <- apply(binMat, 2, function(x) sum(x))
  samPAllele <- samPAllele/(2*nrow(binMat))
  
  names(samPAllele) <- paste("p", c(1:ncol(binMat)), "s", sep="")
  
  return(samPAllele)
}


# Same as above but calculates/extracts the expected margins instead of the sample margins
extractSnpMarginalExpDfVec <- function(pAllele) {

  # Multivariate Frequencies
  nSnp <- length(pAllele)

  marginalFreqsMat <- sapply(pAllele, function(x) {
    c((1 - x) ^ 2, 2 * x * (1 - x), x ^ 2)
  })
  
  marginalFreqsVec <- as.vector(marginalFreqsMat)
  
  freqName <- rep( c("Freq0Expe", "Freq1Expe", "Freq2Expe"),nSnp)
  snpName <- as.vector(sapply(c(1:nSnp), function(x) rep(x,3)))
  snpFreqName <- paste("SNP", snpName, freqName, sep="")
  names(marginalFreqsVec) <- snpFreqName
  return(marginalFreqsVec)
}


# 
# pAllele <- c(0.1, 0.2, 0.3)
# binMat <- cbind(
#   rbinom(n = 1e5, size = 2, prob = pAllele[1]),
#   rbinom(n = 1e5, size = 2, prob = pAllele[2]),
#   rbinom(n = 1e5, size = 2, prob = pAllele[3])
# )
# # 
# df1 <- extractSnpMarginalSamDfVec(binMat)
# df1
# df2 <- extractSnpMarginalExpDfVec(pAllele)
# df2
# 
# cbind(df2,df1)

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

