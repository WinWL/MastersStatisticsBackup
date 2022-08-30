
rm(list=ls())
library(hapsim)
library(dplyr)
# setwd("F:/Dropbox/Dropbox/School/SS_Thes/Exp/SimulationOrCalculations")
# setwd("D:/Dropbox/Dropbox/School/SS_Thes/Exp/SimulationOrCalculations")
source("../SourceCode/Generate2SnpBernoulli.R")










set.seed(100)
data(ACEdata)

# Take a subset of a dataset from hapsim package
subDat <- cbind(ACEdata[,1], ACEdata[,5])
subDat
cor(subDat[,1], subDat[,2])

# Construct haplotype object
hapSubDat <- haplodata(subDat)
hapSubDat

# Simulate two independent haplotypes and add to get genotype
nBinObs <- 100000
simHapDatBern1 <- haplosim(nBinObs, hapSubDat)
simHapDatBern2 <- haplosim(nBinObs, hapSubDat)
simHapDatBernAdded <- simHapDatBern1$data + simHapDatBern2$data
simHapDatBernAdded
# 2 - xxdat - yydat so that the 0,1,2 are aligned with the generate2SnpBernoulli
# 

# Simulate using adding bernoulli method with the same sample values as hapsim
p1 <- 1 - hapSubDat$freqs[1]; p2 <- 1 - hapSubDat$freqs[2]; rho <- hapSubDat$cor[1,2]
sim2SnpBern <- Generate2SnpBernoulli(p1, p2, rho, nBinObs)
sim2SnpBern

apply(simHapDatBernAdded, 2, function(col){
  table(col)/nBinObs
})
apply(sim2SnpBern, 2, function(col){
  table(col)/nBinObs
})

(1-p1)^2; 2*(1-p1)*p1; p1^2
(1-p2)^2; 2*(1-p2)*p2; p2^2

