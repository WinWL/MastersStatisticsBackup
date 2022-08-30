# rm(list=ls())
# setwd("D:/Dropbox/Dropbox/School/SS_Thes/Exp/SimulationOrCalculations")
# setwd("F:/Dropbox/Dropbox/School/SS_Thes/Exp/SimulationOrCalculations")





source("../SourceCode/GenerateAnySnpOptParamPmfBinBern.R")
source("../SourceCode/calcOptParamPmfBinBernAnySnp.R")

# Pmf from the optimizing bernoulli pmf method
p1 <- 0.1; p2 <- 0.12; rho <- 0.1
pAllele <- c(p1, p2)
matTarCor <- diag(2)
matTarCor[upper.tri(matTarCor)] <- c(rho)
matTarCor[lower.tri(matTarCor)] <- c(rho)
optBernPmf <- calcOptPmfAnySnpBernoulli(pAllele, matTarCor)


# Pmf from the solving bernoulli pmf method
a0 <-  rho * sqrt(p1*p2*(1-p1)*(1-p2)) + (1-p1)*(1-p2)
addBernPmf <-
  c(
    `(0,0)` = a0,
    `(0,1)` = 1 - p1 - a0,
    `(1,0)` = 1 - p2 - a0,
    `(1,1)` = a0 + p1 + p2 - 1
  )

addBernPmf
optBernPmf
