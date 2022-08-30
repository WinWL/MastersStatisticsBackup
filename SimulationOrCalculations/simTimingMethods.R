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
library(gdata)
library(matrixcalc)

# rm(list=ls())
# setwd("D:/Dropbox/Dropbox/School/SS_Thes/Exp/SimulationOrCalculations")
# setwd("F:/Dropbox/Dropbox/School/SS_Thes/Exp/SimulationOrCalculations")


source("../SourceCode/CalcAdjustMadsenInverseEigen.R")
source("../SourceCode/CalcMadsenSigmaZMat.R")
source("../SourceCode/SolveMadsenSigmaEqn.R")
source("../SourceCode/TransNormMatToBin.R")
source("../SourceCode/GenerateSnpMadsenBirkes.R")

source("../SourceCode/NamingForBinomialAndBernoulliPmfVec.R")
source("../SourceCode/ExtractSnpMarginalFreq.R")
source("../SourceCode/ExtractSnpMultivarFreq.R")
source("../SourceCode/NamingForBinomialAndBernoulliPmfVec.R")
source("../SourceCode/ExtractSummaryDataframe.R")

source("../SourceCode/GetColNamesSimData.R")

source("../SourceCode/GenerateBag.R")
source("../SourceCode/GenerateAnySnpOptParamPmfBinBern.R")
source("../SourceCode/calcOptParamPmfBinBernAnySnp.R")

source("../SourceCode/NamingForBinomialAndBernoulliPmfVec.R")
source("../SourceCode/NamingForBinomialAndBernoulliPmfVec.R")



###########################################
# Calculating Simulating Time for methods generating 2 SNPs
p1Choice <- c(0.1)
p2Choice <- c(0.3)

nBinObsChoice <- c(1e3, 1e4, 1e6) # number of observations to generate in a single sample/generation
targetBinCorrVec <- c(0.2)
# All Combinations of p1,p2 and target correlations we're interested in
p12TarCor.full <- expand.grid(p1Choice, p2Choice, targetBinCorrVec, nBinObsChoice)

# Filter out redundant choices
p12TarCor <- p12TarCor.full[p12TarCor.full$Var2 >= p12TarCor.full$Var1, ]
p12TarCor

dfTime2Snp <- data.frame()
for (rowind in c(1:nrow(p12TarCor))) {
  # rowind <- 1
  rowval <- p12TarCor[rowind,]
  p1 <- rowval[1] %>% unlist
  p2 <- rowval[2] %>% unlist
  rho <- rowval[3] %>% unlist
  # p1 <- 0.1; p2 <- 0.1; rho <- 1; nBinObs <- 1000;
  # zz <- Generate2SnpBernoulli(p1, p2, rho, nBinObs)
  # zz <- GenerateSnpMadsenBirkes(pAlleleMinor, targetBinCorMat, nBinObs)
  # zz <- GenerateSnpPmf2SnpConstrained(p1, p2, rho, nBinObs)
  # sum(abs(zz[,1] - zz[,2]))
  
  nBinObs <- rowval[4] %>% unlist
  pAlleleMinor <- c(p1, p2)
  targetBinCorMat <- matrix(c(1, rho, rho, 1),
                            nrow = 2,
                            ncol = 2)
  
  # Timing OptBernoulli Method
  timeAll <- sapply(c(1:100), function(x){
    BernoulliTime <- 
      system.time(GenerateAnySnpOptParamBernoulli(pAlleleMinor, targetBinCorMat, nBinObs))
    
    # Timing Madsen and Birkes Method
    MadsenBirkesTime <- 
      system.time(GenerateSnpMadsenBirkes(pAlleleMinor, targetBinCorMat, nBinObs))
    
    # Timing OptBinomial Method
    pmfCalcTime <-
      system.time(GenerateAnySnpOptParamBinomial(pAlleleMinor, targetBinCorMat, nBinObs))
    
    timeAll <- c(MadsenBirkesTime[3], BernoulliTime[3], pmfCalcTime[3])
    return(timeAll)
  })

  # colMeans(t(timeAll))
  # BernoulliTime; MadsenBirkesTime; pmfCalcTime
  
  time <- colMeans(t(timeAll))
  # time
  # methods <- c("MadsenBirkes", "OptBernoulli", "OptBinomial")
  # 
  # dfRow <- data.frame(time, methods, p1, p2, rho, nBinObs)
  # dfTime2Snp <- rbind(dfTime2Snp, dfRow)
  
  dfRow <- c(nBinObs, time)
  dfTime2Snp <- rbind(dfTime2Snp, dfRow)
  
}

# colnames(dfTime2Snp) <- c("Time(s)", "Method", 
#                           "SNP1MAF", "SNP2MAF", "TargetCor12", "NumObs")
colnames(dfTime2Snp) <- c("SampleSize",
                          "MadsenBirkes(s)", "OptBernoulli","OptBinomial")
dfTime2Snp

# setwd("D:/Dropbox/Dropbox/School/SS_Thes/Exp/Results")
# setwd("F:/Dropbox/Dropbox/School/SS_Thes/Exp/Results")
# saveRDS(dfTime2Snp, "simTimingMethods2SnpAug15.rds")

###########################################
# Calculating Simulating Time for methods generating multiple SNPs

nSnps <- c(2, 5, 8)
p1Choice <- c(0.1)
nBinObsChoice <- c(1000) # number of observations to generate in a single sample/generation
targetBinCorrVec <- c(0.1)

# All Combinations of p1,p2 and target correlations we're interested in
p123TarCor.full <- expand.grid(p1Choice, nSnps, targetBinCorrVec,
                               nBinObsChoice)

# Filter out redundant choices
p123TarCor <- p123TarCor.full
p123TarCor

dfTimeManySnp <- data.frame()
# rowind <- 2
for (rowind in c(1:nrow(p123TarCor))) {
  # rowind <- 1
  rowval <- p123TarCor[rowind,]
  p1 <- rowval[1] %>% unlist
  nSnp <- rowval[2] %>% unlist
  tCor <- rowval[3] %>% unlist
  nBinObs <- rowval[4] %>% unlist

  
  # Timing Madsen and Birkes Method
  pAlleleMinor <- rep(p1, nSnp)
  targetBinCorMat <- matrix(tCor,nSnp,nSnp)
  diag(targetBinCorMat) <- 1
  
  # Timing OptBernoulli Method
  timeAll <- sapply(c(1:100), function(x){
    BernoulliTime <- 
      system.time(GenerateAnySnpOptParamBernoulli(pAlleleMinor, targetBinCorMat, nBinObs))
    
    # Timing Madsen and Birkes Method
    MadsenBirkesTime <- 
      system.time(GenerateSnpMadsenBirkes(pAlleleMinor, targetBinCorMat, nBinObs))
    
    # Timing OptBinomial Method
    pmfCalcTime <-
      system.time(GenerateAnySnpOptParamBinomial(pAlleleMinor, targetBinCorMat, nBinObs))
    
    timeAll <- c(MadsenBirkesTime[3], BernoulliTime[3], pmfCalcTime[3])
    return(timeAll)
  })
  
  time <- colMeans(t(timeAll))
  
  dfRow <- c(nSnp, time)
  dfTimeManySnp <- rbind(dfTimeManySnp, dfRow)
#   
#   time <- c(BernoulliTime[3], MadsenBirkesTime[3],pmfCalcTime[3])
#   methods <- c("Bernoulli", "Madsen", "pmfCalc")
#   
#   dfRow <- data.frame(time, methods, nSnp,
#                       p1, 
#                       tCor, nBinObs)
#   
# dfTimeManySnp <- rbind(dfTimeManySnp, dfRow)
}

# colnames(dfTimeManySnp) <- c("Time(s)", "Method", 
#                           "NumberofSnps", "SNPMAF","TargetCor", "NumObs")
colnames(dfTimeManySnp) <- c("TotalSNPs",
                          "MadsenBirkes(s)", "OptBernoulli","OptBinomial")
dfTimeManySnp


# setwd("D:/Dropbox/Dropbox/School/SS_Thes/Exp/Results")
# setwd("F:/Dropbox/Dropbox/School/SS_Thes/Exp/Results")
# saveRDS(dfTimeManySnp, "simTimingMethodsManySnpAug15.rds")

dfTime2SnpLoad <- readRDS("simTimingMethodsBernMadsenPmf2SNP.rds") %>% as.data.frame()
dfTime3SnpLoad <- readRDS("simTimingMethodsBernMadsenPmf3SNP.rds") %>% as.data.frame()

dfTime2SnpLoad %>% colnames()
dfTime3SnpLoad %>% colnames()

dfTime <- rbind(dfTime2SnpLoad, dfTime3SnpLoad)
dfTime %>% View

dfTime
dfTime %>% filter(SnpMAFs == unique(dfTime$SnpMAFs)[2] | SnpMAFs == unique(dfTime$SnpMAFs)[5]) %>% View

 
