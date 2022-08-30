

# Maps values from the normal distribution to the binomial distribution
# based on cut-offs of the multivariate normal
TransformNormToBin <- function(pAlleleMinor, normalGeneratedMat) {
  # pAlleleMinor = vector of minor allele frequencies
  # normalGeneratedMat = matrix of values from multivariate 
  # normal distribution
  
  transBinMat <- sapply(c(1:length(pAlleleMinor)),
                        function(indexBin) {
                          x <- normalGeneratedMat[, indexBin]
                          minAllele <-
                            pAlleleMinor[indexBin]
                          
                          # Transform normal marginal to standard uniform
                          unifVal <-
                            pnorm(x, mean = 0, sd = 1)
                          
                          # Map uniform to binomial using the cutoffs
                          binValue <-
                            TransformUnifToBin(minAllele, unifVal)
                          
                          return(binValue)
                        })
  
  return(data.frame(transBinMat))
}

# Transforms standard uniform values to binomial(2,p)
# given that p is the minor allele frequency
TransformUnifToBin <- function(minAllele, unifVal) {
  # pAlleleMinor = minor Allele Frequency
  # unifMat = Vector of uniform values
  
  unifVal <- as.double(unifVal)
  minAllele <- as.double(minAllele)
  
  # Transform to specified marginal binomial
  # Cutoffs for mapping uniform to binomial
  binCutHomoMinor <- minAllele ^ 2
  binCutHomoMajor <- (1 - minAllele) ^ 2
  binCutHetero <- 2 * minAllele * (1 - minAllele)
  
  # Map uniform to binomial using the cutoffs
  binValue <-
    ifelse(unifVal < binCutHomoMinor,
           2,
           ifelse(unifVal < binCutHetero + binCutHomoMinor,
                  1,
                  0))
  
  return(binValue)
}

#
# pAlleleMinor <- c(0.1, 0.1)
# corMat <-  matrix(c(1, 0, 0, 1), nrow = 2, ncol = 2)
# normalGeneratedMat <- mvrnorm(100000, c(0,0), corMat)
#
# TransformNormToBin(pAlleleMinor, normalGeneratedMat)
