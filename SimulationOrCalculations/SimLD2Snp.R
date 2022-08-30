
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
# library(tidyr)

# rm(list=ls())
# setwd("D:/Dropbox/Dropbox/School/SS_Thes/Exp/Simulation")
# setwd("E:/Dropbox/Dropbox/School/SS_Thes/Exp/Simulation")

source("../SourceCode/ExtractSnpFreq.R")
source("../SourceCode/ExtractSnpBivarFreq.R")
source("../SourceCode/ExtractLD.R")
source("../SourceCode/generateMultinomial.R")

# #######
# # Histogram of Simulaiton results
# snpToShow <- 1
# 
# simVal$BinObs
# hist(simVal$BinObs[,snpToShow]) # Histogram of binomial values
# hist(simVal$NormObs[,snpToShow]) # Histogram of normal values
# simVal$sigmaZ # Shows calculated multivariate normal correlation matrix
# 
# #########
# BinMat <- simVal$BinObs


# pAlleleChoices <- c(0.05, 0.1, 0.15, 0.3, 0.45)
p1Choice <- c(0.01, 0.1, 0.15)
p2Choice <- c(0.01, 0.1)
nBinObs <- 50 # number of observations to generate in a single sample/generation
nSamCorr <- 100 # Number of samples/correlations to generate for each target correlation
Dvec <- seq(0, 0.95, length.out = 20)# Target D we're interested in

# All Combinations of p1,p2 and target correlations we're interested in
p12TarCor.full <- expand.grid(p1Choice, p2Choice, Dvec)

# Filter out redundant choices
p12TarCor <- p12TarCor.full[p12TarCor.full$Var2 <= p12TarCor.full$Var1, ]
p12TarCor

###################
# Parallellized Loop for simulation

cores = detectCores()
cl <- makeCluster(cores[1] - 2) #not to overload your computer. Leaves for regular use
registerDoParallel(cl)
set.seed(102030)

startParT <- Sys.time()
dfGeneratedBinVals <-
  foreach(
    index = 1:nrow(p12TarCor), # Loops over the choice of alleles
    .combine = rbind,
    .multicombine = TRUE,
    .packages = c("MASS", "mvtnorm", "foreach", "reshape2", "dplyr")
  ) %dopar% {
    # Target Binomial correlations we're interested in
      
    # Initialize values
    # index.targetBinCor <- 2
    # index <- 2

    # Major Frequency from minor allele frequency (the input)
    pA <- 1 - p12TarCor[index, 1]
    pB <- 1 - p12TarCor[index, 2]
    D  <- p12TarCor[index, 3]
    
    # 2x2 table of frequencies (Haplotype)
    pAB <- pA*pB + D
    paB <- (1-pA)*pB - D
    pAb <- pA*(1-pB) - D
    pab <- (1-pA)*(1-pB) + D
    
    # 3x3 table of frequencies (Genotype)
    pAABB <- pAB^2
    pAaBB <- 2*pAB*paB
    paaBB <- paB^2
    pAABb <- 2*pAB*pAb
    pAaBb <- 2*pAB*pab + 2*pAb*paB
    paaBb <- 2*paB*pab
    pAAbb <- pAb^2
    pAabb <- 2*pAb*pab
    paabb <- pab^2
    
    # Assemble into single vector
    # Column is AA, Aa, aa
    # Row is BB, Bb, bb
    # Numbered 1-9 in multinomial in order
    pGeno <- c(
      pAABB, pAaBB, paaBB,
      pAABb, pAaBb, paaBb,
      pAAbb, pAabb, paabb
      )
    
    pGenoMat <- matrix(pGeno, ncol = 3, nrow = 3, byrow = T)
    
    
    # Assemble matrix that describess what each value of the
    # multinomial corresponds to a single snp genotype value
    pABmarginalMap <- cbind(
      c(0,1,2,0,1,2,0,1,2), 
      c(0,0,0,1,1,1,2,2,2) 
    )
    
    # We'll partition the [0,1] interval to get values
    # for the multinomial
    genUnifVal <- runif(nBinObs*nSamCorr)  
    genMultinomVals <- genMultinomial(genUnifVal, pGeno)
    
    # Transform the multinomial values into SNP genotype values
    # by using the mapping above
    matSnpVals <- apply(as.matrix(genMultinomVals), 1, function(multinomVal) {
      
      # For each SNP, we apply the mapping to genotype
      retval <- 
        apply(pABmarginalMap, 2, function(snpMap){
            snpGeno <- snpMap[multinomVal]
            
          return(snpGeno)
        })
        
      return(retval)
    }) %>% t()
    
    # Construct dataframe for return values
    
    # Which "sample" each value belongs to
    simRep <- rep(c(1:nSamCorr), nBinObs)
    
    # Include other identifying information with the generated data
    dfRet <- cbind(matSnpVals, simRep, 1-pA, 1-pB, nBinObs, D)

    return(dfRet)
  } %>% as.data.frame()
  
colnames(dfGeneratedBinVals) <- c("Snp1", "Snp2", "SampleNum", "p1", "p2", "NumSampleObs", "TargetD")
endParT <- Sys.time()
startParT - endParT


# Using a loop to summarize the raw data
# Need to do calculations for each "sample" and combination of input parameters (p1, p2, target cor)
p12TarCorSim.full <- expand.grid(p1Choice, p2Choice, Dvec, c(1:nSamCorr))
p12TarCorSim <- p12TarCorSim.full[p12TarCorSim.full$Var2 <= p12TarCorSim.full$Var1, ]
p12TarCorSim %>% dim

startParT <- Sys.time()
DfCorreAlleleSummary <-
  foreach(
    index = 1:nrow(p12TarCorSim), # Loops over the choice of alleles
    .combine = rbind,
    .multicombine = TRUE,
    .packages = c("MASS", "mvtnorm", "foreach", "reshape2", "dplyr")
  ) %dopar% {

    # 
    # index <- 1
    
    # For each input parameter (p1, p2, target Cor)
    # Calculate on SNP frequencies and linkage for each generated sample
    p1 <- p12TarCorSim[index, 1]
    p2 <- p12TarCorSim[index, 2]
    TargetD <- p12TarCorSim[index, 3]
    SampleNum <- p12TarCorSim[index, 4]

    tol <- 0.005
    
    # Filter out
    dfFilter <- dfGeneratedBinVals[
      which(
        dfGeneratedBinVals$p1 >= p1 - tol &
        dfGeneratedBinVals$p1 <= p1 + tol &
        dfGeneratedBinVals$p2 >= p2 - tol &
        dfGeneratedBinVals$p2 <= p2 + tol &
        dfGeneratedBinVals$TargetD >= TargetD - tol &
        dfGeneratedBinVals$TargetD <= TargetD + tol &
        dfGeneratedBinVals$SampleNum <= SampleNum + tol &
        dfGeneratedBinVals$SampleNum >= SampleNum - tol 
      ),]
    
    matBinVals <- as.matrix(dfFilter[,c(1:2)])           # Binomial values only
    
    # SNP Frequencies
    freqsDf <- extractSnpFreqDfVec(matBinVals, p1, p2)
    freqsBivarDf <- extractSnpBivarDfVec(matBinVals, p1, p2)
    
    # Calculate D
    LdDDp <- getLD2Snp(matBinVals, p1, p2)
    
    SampleCor <- cor(matBinVals[,1], matBinVals[,2])    # Sample Correlation
    NumSampleObs <- dfGeneratedBinVals$NumSampleObs[1]  # Parameter: Number of observations in sample 
    

    
    dfRet <- cbind(SampleCor, TargetD, 
                   freqsDf, freqsBivarDf, LdDDp, 
                   p1, p2, NumSampleObs, SampleNum)
    colnames(dfRet) <- c("GeneratedCorrelation", "TargetD", 
                         colnames(freqsDf), colnames(freqsBivarDf), colnames(LdDDp), 
                         "p1Allele", "p2Allele", "NumSampleObs", "SampleNum")
    
    return(dfRet)
  }

endParT <- Sys.time()
startParT - endParT
stopCluster(cl)

zz <- DfCorreAlleleSummary[!complete.cases(DfCorreAlleleSummary), ]
zz$TargetD %>% unique

zz %>% View
zz %>% dim

DfCorreAlleleSummary %>% head
DfCorreAlleleSummary %>% filter(TargetD > 0, TargetD < 0.06) %>% View
DfCorreAlleleSummary %>% dim
# getwd()

# SampleNum <- 1

dfFilter <- dfGeneratedBinVals[
  which(
    dfGeneratedBinVals$p1 >= p1s - tol &
      dfGeneratedBinVals$p1 <= p1s + tol &
      dfGeneratedBinVals$p2 >= p1s - tol &
      dfGeneratedBinVals$p2 <= p2s + tol &
      dfGeneratedBinVals$TargetD >= TargetDs - tol &
      dfGeneratedBinVals$TargetD <= TargetDs + tol &
      dfGeneratedBinVals$SampleNum == SampleNums
  ),]

DfCorreAlleleSummary %>% filter(D > 0)

################
# Saving Simulation Result
# Uncomment to use. 

# setwd("D:/Dropbox/Dropbox/School/SS_Thes/Exp/Results")
# setwd("E:/Dropbox/Dropbox/School/SS_Thes/Exp/Results")
# saveRDS(dfGeneratedBinVals, "SimMadsen2SnpBinVals0To81by10qtxc.rds")
# saveRDS(DfCorreAlleleSummary, "SimMadsen2SnpSumm0To81by10qtxc.rds")


# Load in file. 
# DfCorreAlleleSummaryLoad <- readRDS("SimMadsenP01P15FocusCorr785To810.rds")
# DfCorreAlleleSummaryLoad$TargetCorrelation %>% unlist %>% unique
# Incase of overwriting 
# DfCorreAlleleSummary <- DfCorreAlleleSummaryLoad