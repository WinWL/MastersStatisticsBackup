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
# setwd("D:/Dropbox/Dropbox/School/SS_Thes/Exp/Simulation")
# setwd("E:/Dropbox/Dropbox/School/SS_Thes/Exp/Simulation")

source("../SourceCode/ExtractSnpFreq.R")
source("../SourceCode/CalcMadsenSigmaZMat.R")
source("../SourceCode/SolveMadsenSigmaEqn.R")
source("../SourceCode/TransNormMatToBin.R")
source("../SourceCode/ExtractSnpBivarFreq.R")
source("../SourceCode/ExtractLD.R")
source("../SourceCode/GenerateSnpMadsenBirkes.R")
source("../SourceCode/CalcAdjustMadsenInverseEigen.R")

# Get data files and work with the ones that go from semi-positive definite
# to not semi-positive definite
dfRaw <- readRDS("../Results/calcPositiveDefiniteInOutData/calcPositiveDefiniteInOut3SnpMadsenT020406P01PvaryPVary.rds")
dfClean1 <- apply(dfRaw,2, function(x) unlist(x))

dfRaw <- readRDS("../Results/calcPositiveDefiniteInOutData/calcPositiveDefiniteInOut3SnpMadsenT020406P0203PvaryPVary.rds")
dfClean2 <- apply(dfRaw,2, function(x) unlist(x))

dfRaw <- readRDS("../Results/calcPositiveDefiniteInOutData/calcPositiveDefiniteInOut3SnpMadsen.rds")
dfClean3 <- apply(dfRaw,2, function(x) unlist(x))

dfRaw <- readRDS("../Results/calcPositiveDefiniteInOutData/calcPositiveDefiniteInOut3SnpMadsenT060606P010203PvaryPVary.rds")
dfClean4 <- apply(dfRaw,2, function(x) unlist(x))

dfRaw <- readRDS("../Results/calcPositiveDefiniteInOutData/calcPositiveDefiniteInOut3SnpMadsenT010406TVaryTVaryP010101.rds")
dfClean5 <- apply(dfRaw,2, function(x) unlist(x))

dfClean <- rbind(dfClean1, dfClean2, dfClean3, dfClean4, dfClean5) %>% as.data.frame() %>% distinct()

dfClean$isbinPSD <- as.logical(dfClean$isbinPSD)
dfClean$isnrmPSD <- as.logical(dfClean$isnrmPSD)
dfClean$Labelp1 <- paste("Corre 13, 12, 23:", 
                         dfClean$Tar12, dfClean$Tar13, dfClean$Tar23,
                         "and p1 Allele:", dfClean$p1)
dfClean$Labelt12 <- paste("p Alleles:", 
                          dfClean$p1, dfClean$p2, dfClean$p3,
                          "and Corre 12:", dfClean$Tar12)
dfWorkWith <- dfClean %>% filter(isbinPSD == T & isnrmPSD == F)

# 
MakeSymMat <- function(upperTriVec){
  mat <- diag(length(upperTriVec))
  mat[upper.tri(mat)] <- upperTriVec
  mat[lower.tri(mat)] <- upperTriVec
  return(mat)
}
dfWorkWith %>% head
dfone <- dfWorkWith[1,]

pAll <- c(dfone$p1, dfone$p2, dfone$p3)
tarCor <- c(dfone$Tar12, dfone$Tar13, dfone$Tar23)
tarCorMat <- MakeSymMat(tarCor)

genData <- GenerateSnpMadsenBirkes(pAlleleMinor = pAll, targetBinCorMat = tarCorMat,
                        nBinObs = 10000, useEigAdjust = T)
genTarCor <- c( 
  cor(genData[,1], genData[,2]),
  cor(genData[,1], genData[,3]),
  cor(genData[,2], genData[,3])
  )
genTarCor
tarCor

tarNCor <- c(dfone$Norm12, dfone$Norm13, dfone$Norm23)
tarNCorMat <- MakeSymMat(tarNCor)
CalcAdjustMadsenInverseEigen(tarNCorMat)


# Repeat for all rows of dataframe
# Save Data

# Toss into results file to check on the target correlations matching the samples. Same with the minor allele frequencies
# Clean up 3 SNPs results file. 
# Implement that 5 number summary. 
# Workout what it means to take repeated 5 number summaries. 

# Work out what to say during meeting. 
# Next steps.









# Trying to fix sleep schedule for some time now
# Lack of physical exercise and schedule based on day/night

# Next Steps:
# 3 Snps Madsen and Birkes Method -> Straightforward 
# 3 SNPs, figuring out what's feasible. Correlation and D wise
# 3 SNPs, adding Bernoulli -> ? Some free parameters exist
# Madsen Inverse Eigen Decomp

# Characterize variability around target parameters for each method. 
# Look Into New Coupla?
# 



