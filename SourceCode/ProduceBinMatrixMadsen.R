# Generates a sample of binomial values using the method from Madsen and Birkes 2013

# rm(list=ls())

# Winfield Lai
# Transforms normal values to binomial values based on minor allele frequencies
TransformNormToBin <- function(pAlleleMinor, normalGeneratedMat) {
  # pAlleleMinor = vector of minor allele frequencies for each binomial variable
  # normalGeneratedMat = matrix of normal values to transform to binomial
  
  # Transform normal marginals to unif(0,1) then to binomial
  # indexBin <- 1
  
  # transformedBinomialMat <- apply(normalGeneratedMat, 2, function(x){
  transformedBinomialMat <- sapply(c(1:2), function(indexBin){
    
    x <- normalGeneratedMat[,indexBin]
    # Transform normal marginal to standard uniform
    unifVal <- pnorm(x, mean = 0, sd = 1)
    
    # Transform to specified marginal binomial
    # Cutoffs for mapping uniform to binomial
    minAllele <- pAlleleMinor[indexBin]
    binCutHomoMinor <- minAllele ^ 2
    binCutHomoMajor <- (1 - minAllele) ^ 2
    binCutHetero <- 2 * minAllele * (1 - minAllele)
    
    # Map uniform to binoimal using the cutoffs
    binValue <-
      ifelse(
        unifVal < binCutHomoMinor,
        2,
        ifelse(
          unifVal < binCutHetero + binCutHomoMinor,
          1,
          0
        )
      )
    
    # indexBin <- indexBin + 1
    return(binValue)
  })
  
  
  
  return(data.frame(transformedBinomialMat))
}

# p1A <- 0.1
# p2A <- 0.1
# corre <- 0
# corMat <-  matrix(c(1, corre, corre, 1), nrow = 2, ncol = 2)
# nMat <- mvrnorm(100000, c(0,0), corMat)
# resNtBM <- TransformNormToBin(c(p1A, p2A), nMat)
# tt <- table(resNtBM)
# tt

# Calculate SigmaZ
# The correlation matrix used for simulation (multivariate normal)
CalculateSigmaZ <- function(nObs, nBin, pAlleleMinor, corDesire) {
  # nObs = number of observations to generate
  # nBin = Should be a vector of just 2's. The number of "trials" for binomial
  # pAlleleMinor = vector of minor allele frequencies for each binomial variable
  # corDesire = correlation matrix for the multivariate binomial distribution we want
  
  # Calculate correlation matrix for normals
  corLen <- length(corDesire[,1])
  sigmaZ <- diag(length(corDesire[, 1]))
  
  # Note: Can speed up by running loop in parallel
  # Calculates the correlation matrix
  for (i in c(1:corLen)) {
    for (j in c(1:corLen)) {
      
      # Initialize values needed for calculating delta
      corD <- corDesire[i, j]
      pI <- pAlleleMinor[i]
      pJ <- pAlleleMinor[j]
      nI <- nBin[i]
      nJ <- nBin[j]
      sigmaI2 <- pI * nI * (1 - pI)
      sigmaJ2 <- pJ * nJ * (1 - pJ)
      uI <- nI * pI
      uJ <- nJ * pJ
      
      # For Later: can speed up by only calculating only the lower triangular elements and 
      # retriving them to make the upper triangular elements, delta(i,j) = delta (j,i)
      if (i != j) {
        # Off diagonal elements of correlation matrix
        solvedResult <- solveSigmaEqnForNormToBin(corD, sigmaI2, sigmaJ2, uI, uJ, nI, nJ, pI, pJ)
        sigmaZ[i,j] <- solvedResult
        
      } else {
        # Diagonal elements of correlation matrix
        sigmaZ[i,j] <- 1
      }
      
    }
  }
  
  return(sigmaZ)
}

# Winfield Lai
# Function that numerically estimates delta

# Numerically solve variance for normal marginals
# Returns the correlation supposed to be used
# Assumes R and S are binomial
solveSigmaEqnForNormToBin <- function(corDesire, sigmaI2, sigmaJ2, uI, uJ, nI, nJ, pI, pJ){
  # corDesire is the correlation desired between two binomials
  # sigmaI2, uI, nI, pI are specified parameters for one of the binomials
  # sigmaJ2, uJ, nJ, pJ are specified parameters for one of the binomials
  
  
  # Initialize the values we will be using
  sigmaI <- sqrt(sigmaI2)
  sigmaJ <- sqrt(sigmaJ2)
  
  # Support of distributions. 2 Binomials
  rVals <- c(0,1,2,0,1,2,0,1,2)
  sVals <- c(0,0,0,1,1,1,2,2,2)
  
  # Cdf of R and S values. 
  FiRVec <- pbinom(rVals, nI, pI)
  FjSVec <- pbinom(sVals, nJ, pJ)
  
  # Inverse of Standard Normal of FiR and FiS
  NInvFiRVec <- qnorm(FiRVec)
  NInvFjSVec <- qnorm(FjSVec)
  
  # Function we will be optimizing to find delta
  fx.opt <- function(deltaInit){
    delta <- deltaInit[1]
    
    # Function to create bivariate normal functions with the delta for optimization
    fx.BivarNorm <- function(nFiR, nFiS) {
      sigDelta <- matrix(c(1,delta,delta,1), nrow = 2, ncol = 2)
      pmvnorm(upper = c(nFiR, nFiS), mean = c(0,0), sigma = sigDelta)
    }
    
    # Constructs vector so we can sum the bivariate normal values up
    BiNormVec <- vector(mode = "double", length = length(rVals))
    for (i in c(1:9)) {
      BiNormVec[i] <- fx.BivarNorm(NInvFiRVec[i], NInvFjSVec[i])
    }
    
    # Assembling the equation
    eqSumTermVec <- 1 - FiRVec - FjSVec + BiNormVec
    
    rightEqTerm <- 1/(sigmaI*sigmaJ) * (sum(eqSumTermVec) - uI*uJ)
    leftTerm <- corDesire
    
    # What is being optimizd
    (rightEqTerm - leftTerm)^2
  }
  
  # Delta can be -1 to 1
  opt <- optimize(fx.opt, interval = c(-1,1))
  estimatedDelta <- opt$minimum
  return(estimatedDelta)
}

# Winfield Lai
# Extracts the Binomial (2, p) or SNP frequencies from a sample. 

# Input: binMat: 2 column matrix, entries are 0,1,2
# Input: p1A: probability of minor allele
# Input: p2A: probability of minor allele

# Output: single row df whose columns are the sample 
# and expected frequncies of binomial 0,1,2 values 


# Speed Up: Remove dependency on reshape2. Use table function instead of apply. 
extractSnpFreqDfVec <- function(binMat, p1A, p2A) {
  # library(reshape2)
  
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
# bin1 <- rbinom(n = 30000, size = 2, prob = p1A)
# bin2 <- rbinom(n = 30000, size = 2, prob = p2A)
# binMat <- cbind(bin1, bin2)
# df <- extractSnpFreqDfVec(binMat, p1A, p2A)
# t(df)
# df
# ddf <- cbind(t(df[,c(1:3)]), t(df[,c(4:6)]),t(df[,c(7:9)]), t(df[,c(10:12)]))
# colnames(ddf) <- c("SNP1 Sample", "SNP1 Expected", "SNP2 Sample", "SNP2 Expected")
# rownames(ddf) <- c("Freq0", "Freq1", "Freq2")
# ddf


# library(reshape2)
# library(MASS)
# library(mvtnorm)


# Madsen and birkes method to generate correlated binomial values
# Uses the multivariate normal with a particular correlation matrix
# Correlation matrix for multivariate normal is calculated based of Target Correlation, p1 and p2
# Multivariate normal values are transformed to binomial values based on mapping areas of the marginal
# normal density to marginal binomial density.

# nBinObs <- 1000000
# targetBinCor <- 0
# nBin <- c(2, 2)
# p1Allele <- 0.1
# p2Allele <- 0.15
# pAlleleMinor <- c(p1Allele, p2Allele)
# targetBinCorMat <- matrix(c(1, targetBinCor, targetBinCor, 1),
#                           nrow = 2,
#                           ncol = 2)
# 
# # Calculate correlation matrix for normals
# sigmaZ <- CalculateSigmaZ(nBinObs, nBin, pAlleleMinor, targetBinCorMat)
# 
# # Generate from standard multivaraiate normal with the calculated multivariate normal correlation matrix
# normMeans <- nBin * 0 # mean is just 0
# generatedNormalObs <- mvrnorm(nBinObs, normMeans, sigmaZ)
# 
# # Transform the normal values to binoimal values
# dfBinValues <- TransformNormToBin(
#   pAlleleMinor = pAlleleMinor,
#   normalGeneratedMat = generatedNormalObs)
# 
# tt <- table(dfBinValues)
# tt

# Marginal for rows and their expected frequency
# (sum(tt[1,]))/(sum(tt))
# (1-p1Allele)^2
# 
# (sum(tt[2,]))/(sum(tt))
# 2*p1Allele*(1-p1Allele)
# 
# (sum(tt[3,]))/(sum(tt))
# p1Allele^2


# Marginal for columns and their expected frequency
# (sum(tt[,1]))/(sum(tt))
# (1-p2Allele)^2
# 
# (sum(tt[,2]))/(sum(tt))
# 2*p2Allele*(1-p2Allele)
# 
# (sum(tt[,3]))/(sum(tt))
# p2Allele^2


# Same as above frequencies but using the function I normally use. Some transformations done to make a usable table
# df <- extractSnpFreqDfVec(dfBinValues, p1Allele, p2Allele)
# ddf <- cbind(t(df[,c(1:3)]), t(df[,c(4:6)]),t(df[,c(7:9)]), t(df[,c(10:12)]))
# colnames(ddf) <- c("SNP1 Sample", "SNP1 Expected", "SNP2 Sample", "SNP2 Expected")
# rownames(ddf) <- c("Freq0", "Freq1", "Freq2")
# ddf

# getRhoRangeFreq(0.099285, 0.14971)

