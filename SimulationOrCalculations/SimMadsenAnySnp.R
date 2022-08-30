library(plyr)
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
# library(tidyr)

# rm(list=ls())
# setwd("D:/Dropbox/Dropbox/School/SS_Thes/Exp/SimulationOrCalculations")
# setwd("F:/Dropbox/Dropbox/School/SS_Thes/Exp/SimulationOrCalculations")

source("../SourceCode/SolveMadsenSigmaEqn.R")
source("../SourceCode/CalcMadsenSigmaZMat.R")
source("../SourceCode/CalcAdjustMadsenInverseEigen.R")
source("../SourceCode/TransNormMatToBin.R")
source("../SourceCode/GenerateSnpMadsenBirkes.R")

source("../SourceCode/NamingForBinomialAndBernoulliPmfVec.R")
source("../SourceCode/ExtractSnpMarginalFreq.R")
source("../SourceCode/ExtractSnpMultivarFreq.R")
source("../SourceCode/NamingForBinomialAndBernoulliPmfVec.R")
source("../SourceCode/ExtractSummaryDataframe.R")

source("../SourceCode/GetColNamesSimData.R")


###################
# Parallellized Loop for simulation

cores = detectCores()
cl <- makeCluster(cores[1] - 1) #not to overload your computer. Leaves for regular use
registerDoParallel(cl)
# set.seed(102030)
set.seed(102035)

startParT <- Sys.time()
dfGeneratedBinVals <-
  foreach(
    index = 1:nrow(p123TarCor),
    .combine = rbind,
    .multicombine = TRUE,
    .packages = c("MASS", "mvtnorm", "foreach", "reshape2", "dplyr", "gdata")
  ) %dopar% {

    # Initialize values
    # index <- 6
    pAlleleMinor <- as.double(p123TarCor[index, pAlleleIndex])
    tcVec <- as.double(p123TarCor[index, tarCorIndex])
    targetBinCorMat <- makeSymMatFromVec(as.double(tcVec), length(pAlleleMinor))
    
    
    # Simulate Using Madsen and Birkes method
    binVals <- GenerateSnpMadsenBirkes(pAlleleMinor, targetBinCorMat, nBinObs*nSamCorr, useEigAdjust = T)
    
    # Recover the correlation matrix used in the multivariate normal during the generation process
    sigz <- CalculateSigmaZ(nObs = 1, rep(2, length(pAlleleMinor)), pAlleleMinor, targetBinCorMat)
    
    isBinPSD <- is.positive.definite(targetBinCorMat)
    isNrmPSD <- is.positive.definite(sigz)
    
    # Which "sample" each value belongs to
    simRep <- rep(c(1:nSamCorr), nBinObs)
    
    # Adjusting the p Alleles and target correlation vectors so they can easily combine to 
    # get the return values
    pAlleleMat <- makeMatForCombine(pAlleleMinor, nBinObs*nSamCorr)
    tcMat <- makeMatForCombine(tcVec, nBinObs*nSamCorr)
    
    # Include other identifying information with the generated data
    dfRet <- cbind(binVals, 
                   tcMat,
                   simRep, nBinObs,
                   isBinPSD, isNrmPSD,
                   pAlleleMat)
    
    return(dfRet)
  } %>% as.data.frame()

colnames(dfGeneratedBinVals) <- getMadsenName(length(pAlleleIndex))
# dfGeneratedBinVals %>% head

endParT <- Sys.time()
startParT - endParT

######################
startParT <- Sys.time()

# Get summary dataframe
DfCorreAlleleSummary <- extractSummaryDataFrameInParallel(p123TarCorSim, dfGeneratedBinVals, pAlleleIndex, nBinObs)

endParT <- Sys.time()
startParT - endParT
stopCluster(cl)

browseURL('https://www.youtube.com/watch?v=SMX17bf_qr0&ab_channel=ShadowMusic')
# browseURL('https://www.youtube.com/watch?v=QH2-TGUlwu4')

DfCorreAlleleSummary %>% head
DfCorreAlleleSummary %>% dim
DfCorreAlleleSummary %>% data.frame()
DfCorreAlleleSummaryMadsen <- DfCorreAlleleSummary


zz <- apply(DfCorreAlleleSummary, 2, function(x) unlist(x))
DfCorreAlleleSummary[!complete.cases(DfCorreAlleleSummary),] %>% nrow

# getwd()

# plot(y = DfCorreAlleleSummary$GeneratedCorrelation, x = DfCorreAlleleSummary$SampleNum)
# DfCorreAlleleSummary$GeneratedCorrelation %>% summary

# SampleNum <- 1

# dfFilter <- dfGeneratedBinVals[
#   which(
#     dfGeneratedBinVals$p1Allele == p1 &
#       dfGeneratedBinVals$p2Allele == p2 &
#       dfGeneratedBinVals$TargetCorrelation == TargetCor &
#       dfGeneratedBinVals$SampleNum == SampleNum
#   ),]

################
# Saving Simulation Result
# Uncomment to use. 

# setwd("D:/Dropbox/Dropbox/School/SS_Thes/Exp/Results")
# setwd("F:/Dropbox/Dropbox/School/SS_Thes/Exp/Results")
# saveRDS(dfGeneratedBinVals, "SimMadsen3SnpBinValFisher10kSampleVal.rds")
# saveRDS(DfCorreAlleleSummary, "SimMadsen3SnpSumValFisher10kSampleVal.rds")


# Load in file. 
# DfCorreAlleleSummaryLoad <- readRDS("SimMadsen3SnpJuly22SumVal.rds")
# colnames(DfCorreAlleleSummaryLoad)



# DfCorreAlleleSummaryLoad$TargetCorrelation %>% unlist %>% unique
# Incase of overwriting 
# DfCorreAlleleSummary <- DfCorreAlleleSummaryLoad
