
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

source("../SourceCode/GenerateSnpCopula.R")
source("../SourceCode/TransNormMatToBin.R")

source("../SourceCode/NamingForBinomialAndBernoulliPmfVec.R")
source("../SourceCode/ExtractSnpMarginalFreq.R")
source("../SourceCode/ExtractSnpMultivarFreq.R")
source("../SourceCode/NamingForBinomialAndBernoulliPmfVec.R")
source("../SourceCode/ExtractSummaryDataframe.R")

source("../SourceCode/GetColNamesSimData.R")

p1Choice <- c(0.1)
p2Choice <- c(0.05)
p3Choice <- c(0.05, 0.1, 0.2)

nBinObs <- 50 # number of observations to generate in a single sample/generation
nSamCorr <- 50 # Number of samples/correlations to generate for each target correlation
tc12Vec <- c(1,2,5,10,25,50,100,200,400,800)
# tc13Vec <- 0.60
# tc23Vec <- c(0, 0.125)

# All Combinations of p1,p2 and target correlations we're interested in
p123TarCor.full <- expand.grid(p1Choice, p2Choice, 
                               # p3Choice,
                               tc12Vec
                               # tc13Vec, tc23Vec
)

# Filter out redundant choices
p123TarCor <- p123TarCor.full[p123TarCor.full$Var2 <= p123TarCor.full$Var1, ]
p123TarCor

pAlleleIndex <- c(1:2)
tarCorIndex <- c(3)

###################
# Parallellized Loop for simulation

cores = detectCores()
cl <- makeCluster(cores[1] - 1) #not to overload your computer. Leaves for regular use
registerDoParallel(cl)
set.seed(102030)

startParT <- Sys.time()
dfGeneratedBinVals <-
  foreach(
    index = 1:nrow(p123TarCor),
    .combine = rbind,
    .multicombine = TRUE,
    .packages = c("MASS", "mvtnorm", "foreach", "reshape2", "dplyr", "gdata", "copula")
  ) %do% {

    # Initialize values
    # index <- 1
    pAlleleMinor <- as.double(p123TarCor[index, pAlleleIndex])
    theta <- as.double(p123TarCor[index, tarCorIndex])

    # Simulate Using Madsen and Birkes method
    binVals <- GenerateSnpByCopulaAMH(pAlleleMinor, theta, nBinObs*nSamCorr)

    # Which "sample" each value belongs to
    simRep <- rep(c(1:nSamCorr), nBinObs)
    
    # Adjusting the p Alleles vectors so they can easily combine to 
    # get the return values
    pAlleleMat <- makeMatForCombine(pAlleleMinor, nBinObs*nSamCorr)

    # Include other identifying information with the generated data
    dfRet <- cbind(binVals, 
                   theta,
                   simRep, nBinObs,
                   pAlleleMat)
    
    return(dfRet)
  } %>% as.data.frame()
  
colnames(dfGeneratedBinVals) <- getAMHCopulaName(length(pAlleleIndex))
endParT <- Sys.time()
startParT - endParT


# Using a loop to summarize the raw data
# Need to do calculations for each "sample" and combination of input parameters (p1, p2, target cor)
# All Combinations of p1,p2 and target correlations we're interested in
p123TarCor.full <- expand.grid(p1Choice, p2Choice, 
                               # p3Choice,
                               tc12Vec, 
                               # tc13Vec, tc23Vec,
                               c(1:nSamCorr))

# Filter out redundant choices
p123TarCorSim <- p123TarCor.full[p123TarCor.full$Var2 <= p123TarCor.full$Var1, ]
p123TarCorSim

startParT <- Sys.time()

# Get summary dataframe
DfCorreAlleleSummary <- extractSummaryDataFrameInParallel(p123TarCorSim, dfGeneratedBinVals, pAlleleIndex, nBinObs)

endParT <- Sys.time()
startParT - endParT
stopCluster(cl)
browseURL('https://www.youtube.com/watch?v=SMX17bf_qr0&ab_channel=ShadowMusic')


DfCorreAlleleSummary %>% head
DfCorreAlleleSummary %>% dim

zz <- apply(DfCorreAlleleSummary, 2, function(x) unlist(x))
DfCorreAlleleSummary[!complete.cases(DfCorreAlleleSummary),] %>% nrow
zz

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
# saveRDS(dfGeneratedBinVals, "JoeCopulabinvals.rds")
# saveRDS(DfCorreAlleleSummary, "JoeCopulasumm.rds")


# Load in file. 
# DfCorreAlleleSummaryLoad <- readRDS("MadsenT0T1Size1000by1000p015q01SummVals.rds")
# DfCorreAlleleSummaryLoad$TargetCorrelation %>% unlist %>% unique
# DfCorreAlleleSummaryLoad$p1Allele %>% unlist %>% unique
# Incase of overwriting 
# DfCorreAlleleSummary <- DfCorreAlleleSummaryLoad
