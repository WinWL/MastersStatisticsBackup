# Calculates the Fisher Z COnfidence Interval (C.I.)  and whether the target correlation is within the C.I.

# r <- corVec
# n <- nObs
# alpha <- rep(0.05, 1000)
# target <- rep(0.2, 1000)
calcFisherZCI <- function(r, n, alpha, target = NULL){
  # r is the sample pearson correlation
  # n is the sample size
  # alpha is the desired significance level
  # target are values to check are within C.I.
  # Calculates an equal two-tailed confidence interval 
  # of the Fisher Z tranformation
  
  # Outputs a dataframe of the lower bound, upper bound and whether the target was within the bounds

  # The transformation, variance and z_alpha value
  z <- atanh(r)
  sdz <- 1/(sqrt(n-3))
  za <- qnorm(alpha/2, lower.tail = F)
  
  # Calculating the upper and lower bounds
  lower <- tanh(z - za*sdz)
  upper <- tanh(z + za*sdz) 
  ret <- data.frame(lower, upper)
  colnames(ret) <- c("FisherZLower", "FisherZUpper")
  
  # Whether the target value is within the confidence interval bounds
  # Only checks if a target is provided. 
  if (!is.null(target)){
    within <- (lower <= target) & (target <= upper)
    ret <- data.frame(lower, upper, within)
    colnames(ret) <- c("FisherZLower", "FisherZUpper", "isTargetWithinFisherZCI")
  }
    
  return(ret)
}

# r <- c(0.714); n <- c(27)
# alpha <- c(0.05)
# target <- c(0.5)
# zz <- calcFisherZCI(r, n, alpha, target)
# zz
# 
# sum(zz$isTargetWithinCI)/nrow(zz)
# 
# # Comparing to another library that has the confidence intervals
# library(DescTools)
# CorCI(r, n, conf.level = 0.95, alternative = c("two.sided", "less", "greater"))
# 
# atanh(0.5)
# FisherZ(target)