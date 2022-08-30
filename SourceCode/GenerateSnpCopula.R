# Input: Theta, Minor Allele frequencies and number of observations
# Output: Binomial Matrix and column indicating Kendall's Tau

# source("../SourceCode/TransNormMatToBin.R")
GenerateSnpByCopulaAMH <- function(pAllele, theta, nObs){
  require(copula)
  # Define copula and simulate values
  nSnp <- as.integer(length(pAllele))
  
  unifMat <- sapply(c(1:nObs), function(x){
    # Vsamp <- rgeom(1, 1 - theta) + 1
    Vsamp <- rSibuya (1, 1/theta)
    Rall <- rexp(nSnp)
    # Rall <- rSibuya (nSnp)
    
    uAll <- psiJoe(theta, Rall/Vsamp)
    as.vector(uAll)
  })
  unifMat <- t(unifMat)
  
  # unifMat
  # unifMat %>% typeof
  # unifMat %>% t
  # copulaA <- onacopula("A", C(theta, 1:3))
  # unifMat <- rnacopula(n = nObs, copulaA)
  
  # Transform uniform marginals to binomial
  binMat <- sapply(c(1:ncol(unifMat)), function(indexCol){
    pA <- pAllele[indexCol]
    unifVal <- unifMat[,indexCol]
    
    TransformUnifToBin(pA, unifVal)
  })
  
  # Return binomial values
  return(binMat)
}

psiAMH <- function(theta, t){
  (1 - theta)/(exp(t) - theta)
}

psiJoe <- function(theta, t){
  1 - (1 - exp(-t))^(1/theta)
}

# # library(copula)
# theta <- 0.99
# nObs <- 1000
# pAllele <- c(0.1, 0.2, 0.1)
# GenerateSnpByCopulaAMH(pAllele, theta, nObs)
# sum(binMat[,1])/(2*nrow(binMat))
# sum(binMat[,2])/(2*nrow(binMat))
# sum(binMat[,3])/(2*nrow(binMat))

# 
# dfFilter <- cbind(binMat, pAlleleAll[1], pAlleleAll[2], theta)
# dfFilter %>% head
# 
# dfAllVal <- dfFilter
# binValCol <- c(1:3)
# extValCol <- c(4:6)
# pAlValCol <- c(4:5)
# rowVec <- extractSummaryRowVec(dfAllVal, binValCol, extValCol, pAlValCol)
# rowVec
