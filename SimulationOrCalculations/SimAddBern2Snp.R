# https://stats.stackexchange.com/questions/284996/generating-correlated-binomial-random-variables
# Generating Correlated Binomials by summing up correlated bernoulli

library(dplyr)
library(MASS)
library(mvtnorm)
library(ggplot2)
library(ggthemes)
library(reshape)
library(doParallel)
library(foreach)
library(gtools)
library(gridExtra)

rm(list=ls())
# setwd("D:/Dropbox/Dropbox/School/SS_Thes/Exp/SimulationOrCalculations")
# setwd("E:/Dropbox/Dropbox/School/SS_Thes/Exp/SimulationOrCalculations")

source("../SourceCode/ExtractSnpFreq.R")
source("../SourceCode/ExtractSnpBivarFreq.R")
source("../SourceCode/ExtractLD.R")
source("../SourceCode/Generate2SnpBernoulli.R")


# pAlleleChoices <- c(0.05, 0.1, 0.15, 0.3, 0.45)
p1Choice <- c(0.15, 0.1)
p2Choice <- c(0.1)
nBinObs <- 1000 # number of observations to generate in a single sample/generation
nSamCorr <- 1000 # Number of samples/correlations to generate for each target correlation
targetBinCorrVec <- seq(0, 1, length.out = 11)# Target Binomial correlations we're interested in
targetBinCorrVec <- c(0, 0.1, 0.3, 0.5, 0.75, 0.793, 0.805, 0.85, 0.9)

# All Combinations of p1,p2 and target correlations we're interested in
p12TarCor.full <- expand.grid(p1Choice, p2Choice, targetBinCorrVec)

# Filter out redundant choices
p12TarCor <- p12TarCor.full[p12TarCor.full$Var2 <= p12TarCor.full$Var1, ]
p12TarCor

cores = detectCores()
cl <- makeCluster(cores[1] - 1) #not to overload your computer. Leaves for regular use
registerDoParallel(cl)
set.seed(102030)


# Loops over the different combination of p1, p2 alleles
# For each p1, p2 allele combination, and each target correlation interested in, generate a distribution of correlation coefficients from simulated data. 
startParT <- Sys.time()
dfSim.raw <- 
  foreach(
    index = 1:nrow(p12TarCor), # Loops over the choice of alleles
    .combine = rbind,
    .multicombine = TRUE,
    .packages = c("foreach", "dplyr", "reshape2")
  ) %dopar% {

    # Input values for simulation
    # index <- 2
    p1 <- p12TarCor[index, 1]
    p2 <- p12TarCor[index, 2]
    rho <- p12TarCor[index, 3]
    
    
    # nSamCorr Number repeated simulations/samples for each combination of input values
    # Sectioned this way so some summary statistics are calculated at the same time
    corGenFreqSum <- sapply(c(1:nSamCorr), function(sampleNumIndex) {
      
      # Simulate binomial values
      genBinVals <- Generate2SnpBernoulli(p1, p2, rho, nBinObs)
      x <- genBinVals[,1]
      y <- genBinVals[,2]
      
      # Calculate sample correlation
      GeneratedCorrelation <- cor(x, y)
      
      # CalculateSNP Freq
      snpFreq <- extractSnpFreqDfVec(cbind(x,y), p1, p2) %>% as.matrix()
      snpBivarFreq <- extractSnpBivarDfVec(cbind(x,y), p1, p2) %>% as.matrix()
      
      # Calculate LD
      LdDDp <- getLD2Snp(as.matrix(cbind(x,y)), p1, p2)
      
      # Renaming and binding of values for return
      retDf <- data.frame(GeneratedCorrelation, rho, snpFreq, snpBivarFreq, LdDDp, 
                          p1, p2, nBinObs, sampleNumIndex)
      colnames(retDf) <- c("GeneratedCorrelation", "TargetCorrelation",
                           colnames(snpFreq), colnames(snpBivarFreq), colnames(LdDDp), 
                           "p1Allele", "p2Allele", "NumSampleObs", "SampleNum")
      
      return(retDf)
    }) %>% t()
    
    
    # Return generated dataframe
    return(corGenFreqSum)
  }

endParT <- Sys.time()
startParT - endParT
stopCluster(cl)

dfSim.raw %>% head
dfSim.raw %>% dim

# zz <- apply(dfSim.raw, 2, function(x) unlist(x))
# zz[!complete.cases(zz),] 


# setwd("D:/Dropbox/Dropbox/School/SS_Thes/Exp/Results")
# setwd("E:/Dropbox/Dropbox/School/SS_Thes/Exp/Results")

# saveRDS(dfSim.raw, "AddBernT0T1Size1000by1000p015p01q01sumvals.rds")

# Load in file. 
# dfSim.raw.Load <- readRDS("dfSimP01P15FocusCorr7850To7935.rds")
# dfSim.raw.Load$TargetCorrelation %>% unlist %>% unique


