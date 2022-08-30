# A single file that provides the column names to various different simulation methods
# This keeps the names all visible in one file and easy to modify when needed 


# Same functions to return names
# Used to keep naming consistent and easy to configure
getSnpName <- function(nSnp = c(1:3)){
  return(paste("Snp", nSnp, sep=""))
}

getpAlleleExp <- function(nSnp = c(1:3)){
  return(paste("p", nSnp, "e", sep=""))
}

getpAlleleSam <- function(nSnp = c(1:3)){
  return(paste("p", nSnp, "s", sep=""))
}

getSampleNumAndObsInSample <- function(){
  return(c("SampleIndexNumber","nObsInSample"))
}

getNamingPair <- function(nSnp = c(1:3)){
  grd <- expand.grid(nSnp, nSnp)
  grd <- grd[grd$Var1 < grd$Var2, ]
  return(grd)
}

getTargetCor <- function(nSnp = c(1:3)){
  grd <- getNamingPair(nSnp)
  return(paste("TargetCor", grd[,1], grd[,2], sep=""))
}

# getSnpName()
# getpAlleleExp()
# getpAlleleSam()
# getSampleNumAndObsInSample()
# getNamingPair()
# getTargetCor()

makeMatForCombine <- function(infoVector = c(1:3),
                              numrow = 100) {
  # Takes a vector of information (i.e. minor allele frequencies) and
  # converts it to a matrix. The vector is repeated numrow times to
  # create the matrix
  # Used for cbinding matrices and vectors
  matrix(
    infoVector,
    nrow = numrow,
    ncol = length(infoVector),
    byrow = T
  )
} 

makeSymMatFromVec <- function(vecNum = c(1:3), diagNum = 3){
  # Creates a symmetric matrix given the dimensions and a
  # vector describing the upper triangular values. 
  # The values given are assumed to correspond to the 
  # first row of the upper triangular matrix, then the second row and so on
  
  # vecNum <- c(1:10)
  # diagNum <- 5
  # # 
  # vecNum <- tcVec
  # diagNum <- 5
  
  matSym <- diag(diagNum)
  matSym[lower.tri(matSym)] <- vecNum
  matSym[upper.tri(matSym)] <- t(matSym)[upper.tri(matSym)]
  return(matSym)
}


# nSnp <- 3
getMadsenName <- function(nSnp = 1){
  # Creates the column names for the output of the madsen and birkes 
  # simulation method used
  snpVec <- c(1:nSnp)
  
  snpNam <- getSnpName(snpVec)
  pexNam <- getpAlleleExp(snpVec)
  # psaNam <- getpAlleleSam(snpVec)
  samNam <- getSampleNumAndObsInSample()
  
  tarNam <- getTargetCor(snpVec)
  grd <- getNamingPair(snpVec)
  mulNam <- paste("MultiCor", grd[,1], grd[,2], sep="")
  exInfo <- c("isBinPSD", "isNrmPSD")
  
  return(
    c(snpNam, mulNam, tarNam, samNam, exInfo, pexNam)
    )
}

# nSnp <- 3
getNortaName <- function(nSnp = 1){
  # Creates the column names for the output of the madsen and birkes 
  # simulation method used
  snpVec <- c(1:nSnp)
  
  snpNam <- getSnpName(snpVec)
  pexNam <- getpAlleleExp(snpVec)
  # psaNam <- getpAlleleSam(snpVec)
  samNam <- getSampleNumAndObsInSample()
  
  tarNam <- getTargetCor(snpVec)
  grd <- getNamingPair(snpVec)
  mulNam <- paste("AdjCor", grd[,1], grd[,2], sep="")
  exInfo <- c("isBinPSD", "isNrmPSD")
  
  return(
    c(snpNam, mulNam, tarNam, samNam, exInfo, pexNam)
  )
}
# getNortaName(2)
# getMadsenName(2)

getOptBernName <- function(nSnp = 1){
  # Creates the column names for the output of the madsen and birkes 
  # simulation method used
  snpVec <- c(1:nSnp)
  
  snpNam <- getSnpName(snpVec)
  pexNam <- getpAlleleExp(snpVec)
  # psaNam <- getpAlleleSam(snpVec)
  samNam <- getSampleNumAndObsInSample()
  
  tarNam <- getTargetCor(snpVec)
  grd <- getNamingPair(snpVec)

  return(
    c(snpNam, tarNam, samNam, pexNam)
  )
}

getOptBinName <- function(nSnp = 1){
  # Creates the column names for the output of the madsen and birkes 
  # simulation method used
  getOptBernName(nSnp)
}

getAMHCopulaName <- function(nSnp = 1){
  # Creates the column names for the output of the madsen and birkes 
  # simulation method used
  snpVec <- c(1:nSnp)
  
  snpNam <- getSnpName(snpVec)
  pexNam <- getpAlleleExp(snpVec)
  samNam <- getSampleNumAndObsInSample()
  
  theName <- "Theta"
  grd <- getNamingPair(snpVec)
  
  return(
    c(snpNam, theName, samNam, pexNam)
  )
}


# pAlleleMat <- makeMatForCombine(pAlleleMinor, nBinObs*nSamCorr)
# tcMat <- makeMatForCombine(tcVec, nBinObs*nSamCorr)
# 
# 
# pAlleleMat
# 
# zi1 <- c(p1Allele, p2Allele, p3Allele, 
# tc12, tc13, tc23,
# isBinPSD, isNrmPSD, nBinObs)
# ziM <- makeMatForCombine(zi1, nrow(binVals))
# 
# cbind(binVals, pAlleleMat, tcMat, isBinPSD) %>% head

# makeMatForCombine()

