# Returns the probability mass function (pmf) that respects specified binomial margins with specified target correlation.
# Calculates the pmf by solving linear constratins derived from the target correlation and binomial margins and pinning down (by randomization) any free parameters.
# rm(list = ls())


#########
# Based on correlated binomials
# Calculates the equations associated with the marginal constraints
calcMatMarginalBinomial <- function(numSnp = 3,
                            pAlleles = c(0.1, 0.1, 0.1)) {
    # Calculates the constraint matrix associated to the marginal equations
    # numSnp = the number of snps
    # pAlleles = the minor allele frequency for each snp
    # Returns a typical constraint matrix where the last column 
    # is the RHS of all the equations and the LHS
    # denotes the coefficient on each variable
    # The variables are the frequencies of the multivariate SNPs. 
    # See below for the naming convention.
    
    # First, construct the Left hand side (LHS) of the equations
    snpMatLHS <- matrix(double(), 0, (3 ^ numSnp))
    snpVecRHS <- vector("double", 0)
    
    for (indSnp in c(1:numSnp)) {
      # For each snp,we separately work on the values they take on
      # Note that we let the SNP value take on (1,2,3) instead of (0,1,2)
      # so we don't need to deal with some index issues
      
      # construct the three rows associated to a particular SNP
      snpRows <- t(sapply(c(1:3), function(valSnp) {
        # Example of naming convention for 3 SNPs together
        # 111, 112, 113, 121, 122, 123, 131, 132, 133,
        # 211, 212, 213, 221, 222, 223, 231, 232, 233
        # 311, 312, 313, 321, 322, 323, 331, 332, 333
        
        # If we're interested in finding the SNPs that represent the marginals
        # i.e. if we want the marginal for SNP1 = 1,
        # then we want all of the above naming convention that 
        # have a 1 in the first place (1??)
        # i.e. if we want the marginal for SNP2 = 3,
        # then we want all of the above naming convention that 
        # have a 3 in the second place (?3?)
        
        # We can automate the the way we find the marginals.
        # Note that the SNPs are ordered lowest to highest 
        # (i.e. 111, 112, ..., 332, 333)
        
        # There will be a specific number of blocks
        # This number is related to the index of the SNP we're looking at
        numBlocks <- 3 ^ (indSnp - 1)
        
        # we will have a slew/block of 1's
        blockOneWidth <- 3 ^ (numSnp - indSnp)
        
        # The space between the first block of 1's is 
        # dependent on the value/snp of interest
        blockStartZeros <- (valSnp - 1) * blockOneWidth
        
        # The block of 1's will always be followed by a block of 0's
        blockEndZeros <- (3 - valSnp) * blockOneWidth
        
        ValBlock <- c(rep(0, blockStartZeros),
                      rep(1, blockOneWidth),
                      rep(0, blockEndZeros))
        ValRow <- rep(ValBlock, numBlocks)
        
        return(ValRow)
      }))
      
      # Append the 3 rows for a particular SNP to the matrix of constraints
      snpMatLHS <- rbind(snpMatLHS, snpRows)
      
      # Append the RHS of what the above 3 rows should equal
      pSnp <- pAlleles[indSnp]
      snpVecRHS <- c(snpVecRHS, #Previous values
                     (1 - pSnp) ^ 2, #1
                     2 * (1 - pSnp) * pSnp, #2
                     pSnp ^ 2) #3
    }
    
    # Combine into a matrix and return
    snpMatConstraints <- as.matrix(cbind(snpMatLHS, snpVecRHS))
    return(snpMatConstraints)
  }


# Calculates the equations associated to the target correlation constraints
calcMatTargetCorBinomial <- function(numSnp = 3,
                             matTargetCor = diag(3),
                             pAlleles = c(0.1, 0.1, 0.1)) {
  # Calculates the constraint matrix associated to the Target Correlation equations
  # numSnp = the number of snps
  # matTargetCor = the target correlation matrix
  # pAlleles = the alleles associated to each snp
  # Returns A typical constraint matrix where the last 
  # column is the RHS of all the equations and the LHS
  # denotes the coefficient on each variable
  # The variables are the frequencies of the multivariate SNPs. 
  # See below for the naming convention.
  
  # The snp values that are used in target correlation calculation
  valSnpXYmat <- matrix(c(2, 2, 3, 3,
                          2, 3, 2, 3), 2, 4, T)
  
  # The coefficients to be placed on the variables
  valXYTargetCor <- c(1, 2, 2, 4)

  # We will be appending to construct the matrix
  snpLHSMatTarCor <- matrix(double(), 0, (3 ^ numSnp))
  snpRHSVecTarCor <- vector("double", 0)
  
  for (indSnpX in c(1:numSnp)) {
    for (indSnpY in c((indSnpX):numSnp)) {
      if (isTRUE(all.equal(indSnpX, indSnpY)))
        next #Skip if we look a SNP paired with itself
      # Looking at all pairs of SNPs
      # SnpX is the on a lower index (more left) than SnpY
      
      # For each pair of SNPs
      # we will assemble the row for each value used in the target correlation
      # separately
      
      
      # Construct a matrix where each row indicates where to place
      # The appropriate value a single pair of values 
      # involved in the target correlation
      # i.e. the snp values of 2-2, 2-3, 3-2 and 3-3 are 
      # all used in the target correlation
      # The contribution of each pair is calculated separately and 
      # then added to get a separate row
      separateConstraint <- t(sapply(c(1:4), function(indXY) {
        valSnpX <- valSnpXYmat[1, indXY]
        valSnpY <- valSnpXYmat[2, indXY]
        valPlace <- valXYTargetCor[indXY]
        
        bigBlockNum <- 3 ^ (indSnpX - 1)
        bigBlockWidth <- 3 ^ (numSnp - indSnpX)
        
        bigBlockStartZeros <- (valSnpX - 1) * bigBlockWidth
        
        # Constructing smaller block inside big block
        smlBlockNum <- 3 ^ (indSnpY - indSnpX - 1)
        smlBlockoneWidth <- 3 ^ (numSnp - indSnpY)
        
        smlBlockStartZeros <- (valSnpY - 1) * smlBlockoneWidth
        smlBlockOnes <- smlBlockoneWidth
        smlBlockEndZeros <- (3 - valSnpY) * smlBlockoneWidth
        
        smlblock <- c(
          rep(0, smlBlockStartZeros),
          rep(valPlace, smlBlockOnes),
          rep(0, smlBlockEndZeros)
        )
        smlBlockRow <- rep(smlblock, smlBlockNum)
        
        bigBlockOnes <- smlBlockRow
        
        bigBlockEndZeros <- (3 - valSnpX) * bigBlockWidth
        
        # Constructing final big block
        bigBlock <- c(rep(0, bigBlockStartZeros),
                      rep(smlBlockRow, 1),
                      rep(0, bigBlockEndZeros))
        
        bigBlockRow <- rep(bigBlock, bigBlockNum)
        
        return(bigBlockRow)
      }))
      
      # collapse the above matrix into a row vector to be 
      # used in linear optimization
      singleRowConstraint <- apply(separateConstraint, 2, sum)
      snpLHSMatTarCor <- rbind(snpLHSMatTarCor, singleRowConstraint)
      
      
      # RHS of target correlation is targetCor * sqrt(VarX*VarY) + eX*EY
      pX <- pAlleles[indSnpX]
      pY <- pAlleles[indSnpY]
      nXY <- 2
      rhseXYTarCor <-
        matTargetCor[indSnpX, indSnpY] * 
        sqrt(nXY * pX * (1 - pX) * nXY * pY * (1 - pY)) +
        (nXY * pX * nXY * pY)
      
      snpRHSVecTarCor <-
        c(snpRHSVecTarCor, rhseXYTarCor)
    }
  }
  
  snpMatTarCor <- as.matrix(cbind(snpLHSMatTarCor, snpRHSVecTarCor))
  # snpMatTarCor
  return(snpMatTarCor)
}


# Combines the marginal and target correlation constraints into a single matrix
calcMatPmfConstaintBinomial <-function(numSnp = 3,
                              matTargetCor = diag(3),
                              pAlleles = c(0.1, 0.1, 0.1)) {
  # Calculates the constraint matrix associated to 
  # the Target Correlation equations
  # numSnp = the number of snps
  # matTargetCor = the target correlation matrix
  # pAlleles = the alleles associated to each snp
  # Returns A typical constraint matrix where the 
  # last column is the RHS of all the equations and the LHS
  # denotes the coefficient on each variable
  # The variables are the frequencies of the multivariate SNPs. 
  # See below for the naming convention.
  
  matMarg <- calcMatMarginalBinomial(numSnp, pAlleles)
  matTarC <- calcMatTargetCorBinomial(numSnp, matTargetCor, pAlleles)
  
  snpMat <- as.matrix(rbind(matMarg, matTarC))
  return(snpMat)
}

# Calculates the pmf based on optimizing the set of constraints as given by the marginals and target correlation
calcOptPmfAnySnpBinomial <- function(pAllele, matTarCor){
  require(lpSolve)
  # pAllele = vector in minor allele frequencies
  # matTarCor = correlation matrix
  
  numSnp <- length(pAllele)
  oriMatSnp <- calcMatPmfConstaintBinomial(length(pAllele), matTarCor, pAllele)
  
  # Need to separate matrix to LHS and RHS and add the unity constraint 
  LHSOriMatSnp <- oriMatSnp[,-ncol(oriMatSnp)]
  RHSOriMatSnp <- as.vector(oriMatSnp[,ncol(oriMatSnp)])
  
  # Add unity constraint to the above.
  # RHS only needs a 1 added
  uniRHSOriMatSnp <- c(RHSOriMatSnp, 1)
  
  # LHS needs an extra row and column, last diagonal entry needs to be 1
  uniLHSOriMatSnp1 <- rbind(LHSOriMatSnp, 0)
  uniLHSOriMatSnp2 <- cbind(uniLHSOriMatSnp1, 0)
  uniLHSOriMatSnp2[nrow(uniLHSOriMatSnp2), ncol(uniLHSOriMatSnp2)] <- 1
  uniLHSOriMatSnp <- uniLHSOriMatSnp2
  
  # Objective function
  # All of the variables added should equal the unity (last) variable
  obj <- c(
    rep(1,ncol(LHSOriMatSnp))
    ,-1)
  
  constranints_direction  <- rep("=", length(uniRHSOriMatSnp))
  
  optimum <-  lp(direction = "max",
                 objective.in = obj,
                 const.mat = uniLHSOriMatSnp,
                 const.dir = constranints_direction,
                 const.rhs = uniRHSOriMatSnp,
                 all.int = F)
  sol <- optimum$solution
  
  # Remove last entry since that's the total sum 
  # Remaining vector contains the entries in the PMF of the target
  # random variable
  pmfSolve <- sol[-length(sol)]
  
  # Name the vector according to the naming convention 
  names(pmfSolve) <- namingForBinomialPmf(numSnp)
  return(pmfSolve)
}


# p1 <- 0.1; p2 <- 0.3; p3 <- 0.1
# t12 <- 0.2; t13 <- 0.2; t23 <- 0
# pAllele <- c(p1, p2, p3)
# matTarCor <- diag(3)
# matTarCor[upper.tri(matTarCor)] <- c(t12,t13,t23)
# matTarCor[lower.tri(matTarCor)] <- c(t12,t13,t23)
# zzpmf <- calcOptPmfAnySnpBinomial(pAllele, matTarCor)
# 
# zz <- zzpmf
# zz
# 
# # p_1
# zz[c(1:9)] %>% sum()
# zz[c(10:18)] %>% sum()
# zz[c(19:27)] %>% sum()
# 
# # p_2
# zz[c(1,2,3,10,11,12,19,20,21)] %>% sum()
# zz[c(4,5,6,13,14,15,22,23,24)] %>% sum()
# zz[c(7,8,9,16,17,18,25,26,27)] %>% sum()
# 
# 
# 
# # p_3
# zz[c(1,4,7,10,13,16,19,22,25)] %>% sum()
# zz[c(2,5,8,11,14,17,20,23,26)] %>% sum()
# zz[c(3,6,9,12,15,18,21,24,27)] %>% sum()
# 
# # Cor12
# (zz[c(13, 14, 15)] %>% sum()) +
#   (2*zz[c(16, 17, 18)] %>% sum()) +
#   (2*zz[c(22, 23, 24)]  %>% sum()) +
#   (4*zz[c(25,26, 27)] %>% sum()) %>% sum()
# (0.12 - 4*0.1*0.3)/(2*sqrt(0.1*0.9*0.3*0.7))
# 
# # Cor13
# (zz[c(11, 14, 17)] %>% sum()) +
#   (2*zz[c(12, 15, 18)] %>% sum()) +
#   (2*zz[c(20, 23, 26)]  %>% sum()) +
#   (4*zz[c(21,24, 27)] %>% sum()) %>% sum()
# (0.076 - 4*0.1*0.1)/(2*sqrt(0.1*0.9*0.1*0.9))
# 
# 
# # Cor23
# (zz[c(5, 14, 23)] %>% sum()) +
#   (2*zz[c(6, 15, 24)] %>% sum()) +
#   (2*zz[c(8, 17, 26)]  %>% sum()) +
#   (4*zz[c(9,18, 27)] %>% sum()) %>% sum()
# (0.12 - 4*0.1*0.3)/(2*sqrt(0.1*0.9*0.3*0.7))
# 
# zz
# 
# genBagVals <- sample.int(27, 50 , replace=TRUE, prob=zz)
# genBagVals
# 
# sample.int(27, 50 , replace=TRUE, prob=zz)
# 
# samples <- matrix(sample(27, 5000, replace = TRUE, prob = zz), ncol = 50)
# samples
# 
# sample.int()
# 
# p1 <- 0.1; p2 <- 0.3; p3 <- 0.1
# t12 <- 0; t13 <- 0.2; t23 <- 0
# pAllele <- c(p1, p2, p3)
# matTarCor <- diag(3)
# matTarCor[upper.tri(matTarCor)] <- c(t12,t13,t23)
# matTarCor[lower.tri(matTarCor)] <- c(t12,t13,t23)
# zzpmfbern <- calcOptPmfAnySnpBernoulli(pAllele, matTarCor)
# 
# zzpmfbern
# 
# p <- 0.9480
# sqrt(p*(1-p)/1000)

#######
# Based on Adding two bernoulli's 
# Calculates the equations associated with the marginal constraints
calcMatMarginalBernoulli <- function(numSnp = 3,
                                    pAlleles = c(0.1, 0.1, 0.1)) {
  # Calculates the constraint matrix associated to the marginal equations
  # numSnp = the number of snps
  # pAlleles = the minor allele frequency for each snp
  # Returns A typical constraint matrix where the 
  # last column is the RHS of all the equations and the LHS
  # denotes the coefficient on each variable
  # The variables are the frequencies of the multivariate SNPs. 
  # See below for the naming convention.
  
  # First, construct the Left hand side (LHS) of the equations
  snpMatLHS <- matrix(double(), 0, (2 ^ numSnp))
  snpVecRHS <- vector("double", 0)
  
  for (indSnp in c(1:numSnp)) {
    # For each snp,we separately work on the values they take on
    # Note that we let the SNP value take on (1,2,3) instead of (0,1,2)
    # so we don't need to deal with some index issues
    
    # construct the three rows associated to a particular SNP
    snpRows <- t(sapply(c(1:2), function(valSnp) {
      # Example of naming convention for 3 SNPs together
      # 111, 112, 121, 122, 211, 212, 221, 222
      
      # If we're interested in finding the SNPs that represent the marginals
      # i.e. if we want the marginal for SNP1 = 0,
      # then we want all of the above naming convention that have a 1 
      # in the first place (1)
      # i.e. if we want the marginal for SNP2 = 2,
      # then we want all of the above naming convention that 
      # have a 3 in the second place
      
      # We can automate the the way we find the marginals.
      # Note that the SNPs are ordered lowest to highest 
      # (i.e. 111, 112, ..., 221, 222)
      
      # There will be a specific number of blocks
      # This number is related to the index of the SNP we're looking at
      numBlocks <- 2 ^ (indSnp - 1)
      
      # we will have a slew/block of 1's
      blockOneWidth <- 2 ^ (numSnp - indSnp)
      
      # The space between the first block of 1's is dependent 
      # on the value/snp of interest
      blockStartZeros <- (valSnp - 1) * blockOneWidth
      
      # The block of 1's will always be followed by a block of 0's
      blockEndZeros <- (2 - valSnp) * blockOneWidth
      
      ValBlock <- c(rep(0, blockStartZeros),
                    rep(1, blockOneWidth),
                    rep(0, blockEndZeros))
      ValRow <- rep(ValBlock, numBlocks)
      
      return(ValRow)
    }))
    
    # Append the 2 rows for a particular SNP to the matrix of constraints
    snpMatLHS <- rbind(snpMatLHS, snpRows)
    
    # Append the RHS of what the above 2 rows should equal
    pSnp <- pAlleles[indSnp]
    snpVecRHS <- c(snpVecRHS, #Previous values
                   (1 - pSnp), #1
                   pSnp) #2
  }
  
  # Combine into a matrix and return
  snpMatConstraints <- as.matrix(cbind(snpMatLHS, snpVecRHS))
  return(snpMatConstraints)
}


# Calculates the equations associated to the target correlation constraints
calcMatTargetCorBernoulli <- function(numSnp = 3,
                                     matTargetCor = diag(3),
                                     pAlleles = c(0.1, 0.1, 0.1)) {
  # Calculates the constraint matrix associated to the Target Correlation equations
  # numSnp = the number of snps
  # matTargetCor = the target correlation matrix
  # pAlleles = the alleles associated to each snp
  # Returns a typical constraint matrix where the last 
  # column is the RHS of all the equations and the LHS
  # denotes the coefficient on each variable
  # The variables are the frequencies of the multivariate 
  # SNPs. See below for the naming convention.
  
  # The snp values that are used in target correlation calculation
  valSnpXYmat <- matrix(c(2,
                          2), 2, 1, T)
  
  # The coefficients to be placed on the variables
  valXYTargetCor <- c(1)
  
  # We will be appending to construct the matrix
  snpLHSMatTarCor <- matrix(double(), 0, (2 ^ numSnp))
  snpRHSVecTarCor <- vector("double", 0)
  
  for (indSnpX in c(1:numSnp)) {
    for (indSnpY in c((indSnpX):numSnp)) {
      if (isTRUE(all.equal(indSnpX, indSnpY)))
        next #Skip if we look a SNP paired with itself
      # Looking at all pairs of SNPs
      # SnpX is the on a lower index (more left) than SnpY
      
      
      # Construct a matrix where each row indicates where to place
      # The appropriate value a single pair of values involved 
      # in the target correlation
      # i.e. the snp values of 2-2, 2-3, 3-2 and 3-3 are 
      # all used in the target correlation
      # The contribution of each pair is calculated separately 
      # and then added to get a separate row
      separateConstraint <- t(sapply(c(1), function(indXY) {
        valSnpX <- valSnpXYmat[1, indXY]
        valSnpY <- valSnpXYmat[2, indXY]
        valPlace <- valXYTargetCor[indXY]
        
        bigBlockNum <- 2 ^ (indSnpX - 1)
        bigBlockWidth <- 2 ^ (numSnp - indSnpX)
        
        bigBlockStartZeros <- (valSnpX - 1) * bigBlockWidth
        
        # Constructing smaller block inside big block
        smlBlockNum <- 2 ^ (indSnpY - indSnpX - 1)
        smlBlockoneWidth <- 2 ^ (numSnp - indSnpY)
        
        smlBlockStartZeros <- (valSnpY - 1) * smlBlockoneWidth
        smlBlockOnes <- smlBlockoneWidth
        smlBlockEndZeros <- (2 - valSnpY) * smlBlockoneWidth
        
        smlblock <- c(
          rep(0, smlBlockStartZeros),
          rep(valPlace, smlBlockOnes),
          rep(0, smlBlockEndZeros)
        )
        smlBlockRow <- rep(smlblock, smlBlockNum)
        
        bigBlockOnes <- smlBlockRow
        
        bigBlockEndZeros <- (2 - valSnpX) * bigBlockWidth
        
        # Constructing final big block
        bigBlock <- c(rep(0, bigBlockStartZeros),
                      rep(smlBlockRow, 1),
                      rep(0, bigBlockEndZeros))
        
        bigBlockRow <- rep(bigBlock, bigBlockNum)
        
        return(bigBlockRow)
      }))
      
      # collapse the above matrix into a row vector to be 
      # used in linear optimization
      singleRowConstraint <- apply(separateConstraint, 2, sum)
      snpLHSMatTarCor <- rbind(snpLHSMatTarCor, singleRowConstraint)
      
      
      # RHS of target correlation is targetCor * sqrt(VarX*VarY) + eX*EY
      pX <- pAlleles[indSnpX]
      pY <- pAlleles[indSnpY]
      nXY <- 1
      rhseXYTarCor <-
        matTargetCor[indSnpX, indSnpY] * 
        sqrt(nXY * pX * (1 - pX) * nXY * pY * (1 - pY)) +
        (nXY * pX * nXY * pY)
      
      snpRHSVecTarCor <-
        c(snpRHSVecTarCor, rhseXYTarCor)
    }
  }
  
  snpMatTarCor <- as.matrix(cbind(snpLHSMatTarCor, snpRHSVecTarCor))
  return(snpMatTarCor)
}


# Combines the marginal and target correlation constraints into a single matrix
calcMatPmfConstaintBernoulli <- function(numSnp = 3,
                                       matTargetCor = diag(3),
                                       pAlleles = c(0.1, 0.1, 0.1)) {
  # Calculates the constraint matrix associated to the Target Correlation equations
  # numSnp = the number of snps
  # matTargetCor = the target correlation matrix
  # pAlleles = the alleles associated to each snp
  # Returns A typical constraint matrix where the last column 
  # is the RHS of all the equations and the LHS
  # denotes the coefficient on each variable
  # The variables are the frequencies of the multivariate SNPs. 

  
  matMarg <- calcMatMarginalBernoulli(numSnp, pAlleles)
  matTarC <- calcMatTargetCorBernoulli(numSnp, matTargetCor, pAlleles)
  
  snpMat <- as.matrix(rbind(matMarg, matTarC))
  return(snpMat)
}



# Calculates the pmf based on optimizing the set of constraints as given by the marginals and target correlation
calcOptPmfAnySnpBernoulli <- function(pAllele, matTarCor){
  require(lpSolve)
  # pAllele = vector of minor allele frequencies 
  # matTarCor = correlation matrix
  
  numSnp <- length(pAllele)
  oriMatSnp <- calcMatPmfConstaintBernoulli(length(pAllele), matTarCor, pAllele)
  
  # Need to separate matrix to LHS and RHS and add the unity constraint 
  
  LHSOriMatSnp <- oriMatSnp[,-ncol(oriMatSnp)]
  RHSOriMatSnp <- as.vector(oriMatSnp[,ncol(oriMatSnp)])
  
  # Add unity constraint to the above.
  # RHS only needs a 1 added
  uniRHSOriMatSnp <- c(RHSOriMatSnp, 1)
  
  # LHS needs an extra row and column, last diagonal entry needs to be 1
  uniLHSOriMatSnp1 <- rbind(LHSOriMatSnp, 0)
  uniLHSOriMatSnp2 <- cbind(uniLHSOriMatSnp1, 0)
  uniLHSOriMatSnp2[nrow(uniLHSOriMatSnp2), ncol(uniLHSOriMatSnp2)] <- 1
  uniLHSOriMatSnp <- uniLHSOriMatSnp2
  
  # Objective function
  # All of the variables added should equal the unity (last) variable
  obj <- c(
    rep(1,ncol(LHSOriMatSnp))
    ,-1)
  
  constranints_direction  <- rep("=", length(uniRHSOriMatSnp))
  
  optimum <-  lp(direction = "max",
                 objective.in = obj,
                 const.mat = uniLHSOriMatSnp,
                 const.dir = constranints_direction,
                 const.rhs = uniRHSOriMatSnp,
                 all.int = F)
  sol <- optimum$solution
  
  # Remove last entry since that's the total sum 
  pmfSolve <- sol[-length(sol)] 
  
  # Name the vector according to naming convention
  names(pmfSolve) <- namingForBernoullipmf(numSnp)
  return(pmfSolve)
}

# # ###############
p1 <- 0.1; p2 <- 0.15; p3 <- 0.05
t12 <- 0.1;
pAllele <- c(p1, p2, p3)
matTarCor <- diag(3)
matTarCor[upper.tri(matTarCor)] <- c(t12, t12 ,t12)
matTarCor[lower.tri(matTarCor)] <- c(t12, t12, t12)
# 
# eigen(matTarCor)
# 
# # #
# # # pAllele
# # # matTarCor
# cc <- calcOptPmfAnySnpBinomial(pAllele, matTarCor)
# cc
# matrix(cc, 3,3)
# 
# 
# p1c <- c((1-p1)^2, 2*(p1)*(1-p1), p1^2)
# p2c<- c((1-p2)^2, 2*(p2)*(1-p2), p2^2)
# xt <- matrix(p1c,3,1) %*% matrix(p2c,1,3)
# 
dd <- calcOptPmfAnySnpBernoulli(pAllele, matTarCor)
dd


# 
# matrix(dd,2,2)
# matrix(c(1-p1,p1),2,1) %*% matrix(c(1-p2,p2),1,2)

# namingForBinomialPmf <- function(numSnp){
#   # Naming convention follows a specific pattern
#   # they are the ordered (lowest to highest) of what value (1,2,3) each snp (indicated by index) takes on 
#   # For example, if snp1 is 1, snp2 is 3 and snp3 is 2, then the corresponding name is 132.
#   # This is a simple pattern
#   
#   # Construct the names from right to left.
#   # The pattern of 1,2,3 is a repetitive one.
#   nmedoub <- rep(0, 3^numSnp)
#   for (snpInd in c(1:numSnp)) {
#     # snpInd <- 1
#     inner <- 3 ^ (snpInd - 1)
#     outer <- 3 ^ (numSnp - snpInd)
#     eexp <- snpInd - 1
#     
#     toAdd <- (10^eexp)*rep(
#       c(rep(1, inner),
#         rep(2, inner),
#         rep(3, inner)),
#       outer
#     )
#     toAdd
#     nmedoub <- nmedoub + toAdd
#   }
#   
#   nmestri <- paste("x", nmedoub, sep = "")
#   nmestri
# }
# 
# namingForBernoullipmf <- function(numSnp){
#   # Naming convention follows a specific pattern
#   # they are the ordered (lowest to highest) of what value (1,2,3) each snp (indicated by index) takes on 
#   # For example, if snp1 is 1, snp2 is 1 and snp3 is 2, then the corresponding name is 112.
#   # This is a simple pattern
#   
#   # Construct the names from right to left.
#   # The pattern of 1,2 is a repetitive one. 
#   nmedoub <- rep(0, 2^numSnp)
#   for (snpInd in c(1:numSnp)) {
#     # snpInd <- 1
#     inner <- 2 ^ (snpInd - 1)
#     outer <- 2 ^ (numSnp - snpInd)
#     eexp <- snpInd - 1
#     
#     toAdd <- (10^eexp)*rep(
#       c(rep(1, inner),
#         rep(2, inner)),
#       outer
#     )
#     toAdd
#     nmedoub <- nmedoub + toAdd
#   }
#   
#   nmestri <- paste("x", nmedoub, sep = "")
#   nmestri
# }
# 
# # mmat <- matrix(cc, 3,3)
# # mmat
# # rowSums(mmat)
# # colSums(mmat)
# 
