# Purpose: Calculate D and D' for a 2SNP situation

# Input: matBinVals: 2 column matrix, entries are 0,1,2
# Input: p1: probability of minor allele
# Input: p2: probability of minor allele



# Output: Vector containing D, D', r and R^2
getLD2Snp <- function(matBinVals, p1, p2){
  require(genetics)
  require(dplyr)
  genVals <- ifelse(matBinVals == 0, "A/A",
                    ifelse(matBinVals == 1, "A/a",
                           "a/a")) %>% as.data.frame()
  
  # D and D' are only defined if there are different alleles present
  retval <- cbind(NA, NA, NA, NA)
  
  if (length(unique(matBinVals[,1])) > 1 &  
      length(unique(matBinVals[,2])) > 1){
    ldGeno <- LD(makeGenotypes(genVals))
    
    # By default D and D' are matrices
    retval <- cbind(ldGeno$D[1,2], ldGeno$`D'`[1,2], ldGeno$r[1,2], ldGeno$`R^2`[1,2])
  }
  
  
  colnames(retval) <- c("D", "DPrime", "LDr", "LDR2")
  return(retval)
}

# matBinVals <- cbind(
#   rbinom(10000,2,0.05),
#   rbinom(10000,2,0.05)
# )
# matBinVals <- cbind(
#   rbinom(10000,2,0.05),
#   rbinom(10000,2,0.05),
#   rbinom(10000,2,0.05)
# )
# 
# length(unique(matBinVals[,1]))
# 
# getLD2Snp(matBinVals, 1, 0.05)
# zzz <- ldGeno$`D'`
# upper.tri(zzz)
# zzz[upper.tri(zzz)]

