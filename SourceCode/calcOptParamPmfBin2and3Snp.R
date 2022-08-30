# Returns the probability mass function (pmf) that respects specified binomial margins with specified target correlation.
# Calculates the pmf by solving linear constraints derived from the target correlation and binomial margins and pinning down (by randomization) any free parameters.
# rm(list = ls())

calcOptPmf2SnpBinomial <- function(pCol, pRow, targetCor) {
  require(lpSolve)

  # pCol is the Minor allele frequency (MAF) corresponding to one of the snps
  # pRow is the MAF for the other snp
  # It does not matter which snp is pCol or which snp is pRow
  # TargetCor is the target correlation.

  # Check target correlation is within bounds
  # Note Update this with the Frechet Hoeffding bounds

  # Check MAF are between 0 and 0.5 exclusive

  # Calculate the expected value and variances
  pColvar <- 2 * pCol * (1 - pCol)
  pColexp <- 2 * pCol
  pRowvar <- 2 * pRow * (1 - pRow)
  pRowexp <- 2 * pRow

  lhsmat <- matrix(c(
    1, 0, 0, 1, 0, 0, 1, 0, 0, 0,
    0, 1, 0, 0, 1, 0, 0, 1, 0, 0,
    0, 0, 1, 0, 0, 1, 0, 0, 1, 0,
    0, 0, 0, 1, 1, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 1, 1, 1, 0,
    0, 0, 0, 0, 1, 2, 0, 2, 4, 0,
    # 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    # 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1
  ), nrow = 7, byrow = T)


  rhsOfeXYTarCor <- targetCor * sqrt(pColvar * pRowvar) + pColexp * pRowexp

  # rx22 <- 0.1*min(2*pCol*(1-pCol), 2*pRow*(1-pRow), rhsOfeXYTarCor)
  # rx23 <- 0.5*min(pCol^2, 2*pRow*(1-pRow) - rx22, rhsOfeXYTarCor - rx22)

  rhs <- c(
    (1 - pCol)^2, 2 * (1 - pCol) * pCol, pCol^2,
    2 * (1 - pRow) * pRow, pRow^2, rhsOfeXYTarCor,
    # rx22,
    # rx23,
    1
  )


  obj <- c(
    1, 1, 1,
    1, 1, 1,
    1, 1, 1,
    -1
  )
  constranints_direction <- c(
    "=", "=", "=",
    "=", "=", "=",
    "="
  )

  optimum <- lp(
    direction = "max",
    objective.in = obj,
    const.mat = lhsmat,
    const.dir = constranints_direction,
    const.rhs = rhs,
    all.int = F
  )
  optimum
  sol <- optimum$solution
  # sol
  # sol %>% t %>% t
  #
  # sol[-10] %>% sum
  # (lhsmat %*% sol)
  # ((lhsmat %*% sol) - rhs) %>% t %>% t
  # ((lhsmat %*% sol) - rhs) %>% t %>% t %>% sum
  #
  # ss <- optimum$objval
  # ss

  pmfSolve <- sol[-10]

  # pmf <- pmfSolve %>% matrix(nrow = 3, ncol = 3, T)
  # pmf
  # colSums(pmf)
  # rowSums(pmf)
  # EXY <- lhsmat[6,] %*% sol
  #
  # (EXY - pColexp*pRowexp)/sqrt(pRowvar*pColvar)
  # targetCor; pCol; pRow

  names(pmfSolve) <- c(
    "x11", "x12", "x13",
    "x21", "x22", "x23",
    "x31", "x32", "x33"
  )
  return(pmfSolve)
}
# library(dplyr)
# pCol <- 0.1; pRow <- 0.1
# targetCor <- 0.1;
# pmfMat <- calcPmfRndParam2(pCol, pRow, targetCor) %>% matrix(3,3,T)
# pmfMat
# #
# # # Marginals
# colSums(pmfMat)
# c((1 - pCol)^2, 2*(1 - pCol)*pCol, pCol^2)
#
# rowSums(pmfMat)
# c((1 - pRow)^2, 2*(1 - pRow)*pRow, pRow^2)
#
# # pmfMat
# # # Calculate the expected value and variances to check target correlation
# pColvar <- 2*pCol*(1-pCol); pColexp <- 2*pCol; pRowvar <- 2*pRow*(1-pRow); pRowexp <- 2*pRow
# EXY <- 4*pmfMat[3,3] + 2*pmfMat[3,2] + 2*pmfMat[2,3] + pmfMat[2,2]
# (EXY - pColexp*pRowexp)/sqrt(pRowvar*pColvar)
# targetCor; pCol; pRow

# Note that it currently gives negative probabilities outside of the range

# CalcFrecHoefBoundsBinBin(0.1,as.integer(2) ,0.5, as.integer(2))

# Calculates the pmf using the pmf constraint optimization method

calcOptPmf3SnpBinomial <- function(p1,p2,p3, t12, t13, t23){
  require(lpSolve)
  
  # p1,p2 and p3 are the Minor allele frequency (MAF) corresponding to one of the snps
  # t12, t13, t23 are the target correlation between SNPs 1,2 and 3 
  
  # p1 <- 0.1; p2 <- 0.12; p3 <- 0.15
  # t12 <- 0.1; t13 <- 0.13; t23 <- 0.15
  
  # Calculate the expected value and variances
  p1var <- 2*p1*(1-p1)
  p1exp <- 2*p1
  p2var <- 2*p2*(1-p2)
  p2exp <- 2*p2
  p3var <- 2*p3*(1-p3)
  p3exp <- 2*p3
  
  lhsmat <- matrix(c(
    1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, # SNP1 = 0
    0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0, # SNP1 = 1
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,0, # SNP1 = 2
    1,1,1,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0, # SNP2 = 0
    0,0,0,1,1,1,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,1,1,1,0,0,0,0, # SNP2 = 1
    0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,1,1,1,0, # SNP2 = 2
    1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,0, # SNP3 = 0
    0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0, # SNP3 = 1
    0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0, # SNP3 = 2
    0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,2,2,2,0,0,0,2,2,2,4,4,4,0, # t12
    0,0,0,0,0,0,0,0,0,0,1,2,0,1,2,0,1,2,0,2,4,0,2,4,0,2,4,0, # t13
    0,0,0,0,1,2,0,2,4,0,0,0,0,1,2,0,2,4,0,0,0,0,1,2,0,2,4,0, # t23
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1  # sum of all values = 1
  ), nrow = 13, byrow = T)
  
  # c(1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) %>% length()
  
  rhsOfeXYt12 <- t12*sqrt(p1var*p2var) + p1exp*p2exp
  rhsOfeXYt13 <- t13*sqrt(p1var*p3var) + p1exp*p3exp
  rhsOfeXYt23 <- t23*sqrt(p2var*p3var) + p2exp*p3exp
  
  rhs <- c(
    (1 - p1)^2, 2*(1 - p1)*p1, p1^2,
    (1 - p2)^2, 2*(1 - p2)*p2, p2^2, 
    (1 - p3)^2, 2*(1 - p3)*p3, p3^2, 
    rhsOfeXYt12,
    rhsOfeXYt13, 
    rhsOfeXYt23,
    1
  )
  
  
  obj <- c(
    rep(1,27)
    ,-1)
  
  constranints_direction  <- rep("=", 13)
  
  optimum <-  lp(direction = "max",
                 objective.in = obj,
                 const.mat = lhsmat,
                 const.dir = constranints_direction,
                 const.rhs = rhs,
                 all.int = F)
  # optimum
  sol <- optimum$solution
  # sol
  # sol %>% t %>% t
  # 
  # sol[-28] %>% sum
  # (lhsmat %*% sol)
  # ((lhsmat %*% sol) - rhs) %>% t %>% t
  # ((lhsmat %*% sol) - rhs) %>% t %>% t %>% sum
  # 
  # ss <- optimum$objval
  # ss
  
  # Remove last entry since that's the total sum 
  pmfSolve <- sol[-length(sol)] 
  
  # pmf <- pmfSolve %>% matrix(nrow = 3, ncol = 3, T)
  # pmf
  # colSums(pmf)
  # rowSums(pmf)
  # EXY <- lhsmat[6,] %*% sol
  # 
  # (EXY - pColexp*pRowexp)/sqrt(pRowvar*pColvar)
  # targetCor; pCol; pRow
  
  names(pmfSolve) <- c(
    "x111", "x112", "x113", "x121", "x122", "x123", "x131", "x132", "x133",
    "x211", "x212", "x213", "x221", "x222", "x223", "x231", "x232", "x233",
    "x311", "x312", "x313", "x321", "x322", "x323", "x331", "x332", "x333"
  )
  return(pmfSolve)
}

# p1 <- 0.1; p2 <- 0.12; p3 <- 0.15
# t12 <- 0.1; t13 <- 0.13; t23 <- 0.15

# zz <- calcPmfParam3Snp(p1, p2, p3, t12, t13, t23)
# 
# zz
# (zz * 
#     c(0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,1)
#   ) %>% sum
# 
# (1 - p3)^2; 2*(1 - p3)*p3; p3^2
