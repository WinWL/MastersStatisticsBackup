
# Input: a matrix whose rows are observations and columns are binomial values or extra identification information
# numerical vectors indicating which columns are binomial values, extra information, etc 
# dfAllVal: matrix of values
# binValCol: Vector indicating columns of binomial values
# extValCol: vector indicating all columns of non-binomial values
# pAlValCol: vector indicating the columns listing the minor allele frequencies


# Output: A row vector of summary statistics
# extractSummaryRowVec: Row vector containing sample correlation, target correlations,
# marginal sample and expected frequencies, multivariate sample frequencies, sample and expected
# minor allele frequencies

# rm(list=ls())
source("../SourceCode/ExtractSnpMarginalFreq.R")
source("../SourceCode/ExtractSnpMultivarFreq.R")
source("../SourceCode/NamingForBinomialAndBernoulliPmfVec.R")

# dfAllVal <- dfGeneratedBinVals
# head(dfAllVal)
# binValCol <- c(1:3)
# extValCol <- c(4:16)
# pAlValCol <- c(5:7)

extractSummaryRowVec <- function(dfAllVal, binValCol, extValCol, pAlValCol){
  
  # A row of extra information indicated by input values
  # print(dfAllVal[1,])dfGeneratedBinVals
  # browser()
  
  extRowVec <-  as.double(dfAllVal[1, extValCol])
  names(extRowVec) <- colnames(dfAllVal[,extValCol])
  
  # The matrix of Binomial values
  matBinVal <- dfAllVal[, binValCol]
  
  
  # Expected Marginal Frequencies
  pAllele <- as.vector(dfAllVal[1, pAlValCol])
  marginalExpFreq <- extractSnpMarginalExpDfVec(pAllele)

  # Sample Marginal and Multivariate Frequencies
  marginalSamFreq <- extractSnpMarginalSamDfVec(matBinVal)
  multivarSamFreq <- extractSnpMultivarDfVec(matBinVal)
  
  # sample Minor allele frequencies
  pSamAllele <- extractSamPAllele(matBinVal)
  
  # Sample Correlations between all unique pairs 
  corPairs <- expand.grid(binValCol, binValCol)
  corPairs <- corPairs[corPairs$Var1 < corPairs$Var2, ]
  # corPairs
  
  corrValsVec <- apply(corPairs, 1, function(rowVal){
    # print(rowVal)
    i <- rowVal[1]
    j <- rowVal[2]
    
    cor(matBinVal[,i], matBinVal[,j])
  })
  
  names(corrValsVec) <- paste("SampleCor", corPairs[,1], corPairs[,2], sep="")

  # Return a named double vector 
  retVec <- c(corrValsVec,
              marginalSamFreq,
              marginalExpFreq,
              multivarSamFreq,
              pSamAllele,
              extRowVec)
  return(retVec)
}

# dfAllVal <- dfGeneratedBinVals
# head(dfAllVal)
# binValCol <- c(1:3)
# extValCol <- c(4:16)
# pAlValCol <- c(5:7)
# z1 <- extractSummaryRowVec(dfAllVal, binValCol, extValCol, pAlValCol)
# z2 <- extractSummaryRowVec(dfAllVal, binValCol, extValCol, pAlValCol)
# 
# 
# rbind(z1,z2)


extractSummaryDataFrameInParallel <-
  function(pAllTarCorSim,
           dfGeneratedBinVals,
           pAlleleIndex,
           nBinObs,
           coresForRegularUse = 2) {
    
    require(doParallel)
    require(foreach)
    
    # cores = detectCores()
    # cl <-
    #   makeCluster(cores[1] - coresForRegularUse) #not to overload your computer. Leaves for regular use
    # registerDoParallel(cl)
    
    DfCorreAlleleSummary <-
      foreach(
        index = 1:nrow(pAllTarCorSim),
        # Loops over the choice of alleles
        .combine = rbind,
        .multicombine = TRUE,
        .packages = c("MASS", "mvtnorm", "foreach", "reshape2", "dplyr"),
        .export = c("extractSummaryRowVec", "extractSnpMarginalExpDfVec", "extractSnpMarginalSamDfVec", "extractSnpMultivarDfVec",
                    "extractSamPAllele", "namingForBinomialPmf", "addDummyFreq")
      ) %do% {
        #
        # index <- 495
        
        # browser()
        # Calculates summary statistics on each sample
        # Sample size is nObsInSample, so we can simply use batches of row indicies to
        # get the next sample
        
        # nBinObs = number of observations in each sample
        sampleIndexStart <- 1 + nBinObs * (index - 1)
        sampleIndexEnd <- (nBinObs * index)
        
        # Filter out
        dfFilter <-
          dfGeneratedBinVals[c(sampleIndexStart:sampleIndexEnd), ]
        dfAllVal <- dfFilter
        
        
        # Indicies used for to get extract summary statistics
        binValCol <- c(1:length(pAlleleIndex))
        
        extValCol <-
          c((length(pAlleleIndex) + 1):ncol(dfGeneratedBinVals))
        
        pAlValCol <-
          c((ncol(dfGeneratedBinVals) - (length(pAlleleIndex)) + 1):ncol(dfGeneratedBinVals))
        
        rowVec <-
          extractSummaryRowVec(dfAllVal, binValCol, extValCol, pAlValCol)
        return(rowVec)
      }
    
    
  }

# pa <- p123TarCorSim
# ddf <- dfGeneratedBinVals
# pAlleI <- pAlleleIndex
# nbi <- 50
# 
# cores = detectCores()
# cl <-
#   makeCluster(cores[1] - 2) #not to overload your computer. Leaves for regular use
# registerDoParallel(cl)
# system.time(extractSummaryDataFrameInParallel(pa, ddf, pAlleI, nbi))
# stopCluster(cl)
# 
# p1Choice <- c(0.1)
# p2Choice <- c(0.05)
# p3Choice <- c(0.05, 0.1, 0.2)
# 
# nBinObs <- 50 # number of observations to generate in a single sample/generation
# nSamCorr <- 50 # Number of samples/correlations to generate for each target correlation
# tc12Vec <- 0.05
# tc13Vec <- 0.60
# tc23Vec <- c(0, 0.125)
# 
# # All Combinations of p1,p2 and target correlations we're interested in
# p123TarCor.full <- expand.grid(p1Choice, p2Choice, p3Choice,
#                                tc12Vec, tc13Vec, tc23Vec)
# 
# # Filter out redundant choices
# p123TarCor <- p123TarCor.full[p123TarCor.full$Var2 <= p123TarCor.full$Var1, ]
# p123TarCor
# 
# pAlleleIndex <- c(1:3)
# tarCorIndex <- c(4:6)

