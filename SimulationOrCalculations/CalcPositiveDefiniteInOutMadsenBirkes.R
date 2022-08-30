library(dplyr)
library(matrixcalc)
library(doParallel)
library(foreach)
library(mvtnorm)

# rm(list=ls())
# setwd("D:/Dropbox/Dropbox/School/SS_Thes/Exp/Simulation")
# setwd("E:/Dropbox/Dropbox/School/SS_Thes/Exp/Simulation")

source("../SourceCode/CalcMadsenSigmaZMat.R")
source("../SourceCode/SolveMadsenSigmaEqn.R")

# Calculating for 3 SNPs
# Whether a positive definite target correlation matrix
# becomes a non-positive definite correlation matrix after transformation by 
# Madsen and Birkes equation

# Set up target correlations and p alleles
# gridCorSeq <- seq(0.05, 0.5, length.out = 10) %>% round(2)
gridCorSeq <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99)
gridCorSeq <- seq(0, 0.99, length.out = 50) %>% round(2)

pAlleChoi <- seq(0.1, 0.2)
pAlleChoiLong <- seq(0.01, 0.5, length.out = 50) %>% round(2)

matCorPvals <- expand.grid(c(0.1, 0.4, 0.6), gridCorSeq, gridCorSeq,
                           0.1, 0.1, 0.1)
matCorPvals

# Indicies needed to identify input parameters
lenTarcorEndIndex <- 3
lenPEndIndex <- 6


cores = detectCores()
cl <- makeCluster(cores[1] - 2) # not to overload your computer. Leaves for regular use
registerDoParallel(cl)

startParT <- Sys.time()
dfPDInOut3Snp <-
  foreach(
    index = 1:nrow(matCorPvals), # Loops over the choice of alleles
    .combine = rbind,
    .multicombine = TRUE,
    .packages = c("mvtnorm", "dplyr", "matrixcalc")
  ) %dopar% {
    
    # Set up input parameters
    # index <- 166
    n <- lenTarcorEndIndex
    vecT <- unlist(matCorPvals[index, c(1:lenTarcorEndIndex)])
    
    # quickly construct the symmetric correlation matrix from a 
    # vector containing the upper triangular values
    matR <- matrix(0, n, n) 
    matR[upper.tri(matR)] <- vecT
    matBinCor <- matR + diag(n) + t(matR) 
    
    # Vector of minor allele frequencies
    vecPmin <- unlist(matCorPvals[index, c((lenTarcorEndIndex+1):lenPEndIndex)])
    nBin <- rep(2, n)
    
    # Correlation matrix for the multivariate normal, calculated
    # from Madsen and Birkes
    matNrmCor <- CalculateSigmaZ(nObs = 1, nBin, vecPmin, matBinCor)
    vecN <- matNrmCor[upper.tri(matNrmCor)] 
    
    # Check positive definite and record values
    nrmPSD <- is.positive.definite(matNrmCor)
    binPSD <- is.positive.definite(matBinCor)
    
    retval <- data.frame(t(vecT), t(vecN), t(vecPmin), binPSD, nrmPSD)
    return(retval)
  }

endParT <- Sys.time()
startParT - endParT
stopCluster(cl)

colnames(dfPDInOut3Snp) <- c("Tar12", "Tar13", "Tar23", "Norm12", "Norm13","Norm23","p1", "p2", "p3","isbinPSD", "isnrmPSD")

dfDiff3Snp <- dfPDInOut3Snp %>% filter(isbinPSD == T & isnrmPSD == F)
dfDiff3Snp
dfDiff3Snp %>% filter(Tar12 == 0.2, Tar13 == 0.4, Tar23 == 0.5)
dfPDInOut3Snp %>% filter(Tar12 == 0.2, Tar13 == 0.4, Tar23 == 0.5, p2 == 0.45)

dfDiff3Snp %>% nrow

dfDiff3Snp[dfDiff3Snp$isbinPSD == F,]
ddf <- dfDiff3Snp[dfDiff3Snp$isbinPSD == T,]
ddf

# setwd("D:/Dropbox/Dropbox/School/SS_Thes/Exp/Results")
# setwd("E:/Dropbox/Dropbox/School/SS_Thes/Exp/Results")
# saveRDS(dfPDInOut3Snp, "calcPositiveDefiniteInOut3SnpMadsenT010406TVaryTVaryP010101.rds")
# saveRDS(dfPDInOut3Snp, "calcPositiveDefiniteInOut3SnpMadsen.rds")

