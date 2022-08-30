

calcPercentileCI <- function(obs, target = NULL, alpha = 0.05){
  # obs are the sample observations of the target
  # Target is the value that we will check to see if it's in the CI bounds
  # alpha is the desired significance level. Automatically does a two-tailed test
  # Calculates confidence interval based on percentile. 
  # Checks whether the target value is within the confidence interval
  
  
  sortobs <- sort(obs)
  n <- length(obs)
  
  lowerInd <- floor(n*alpha/2)
  upperInd <- ceiling(n*(1- alpha/2))
  
  percentileLower <- sortobs[lowerInd]  
  percentileUpper <- sortobs[upperInd]
  
  ret <- data.frame(percentileLower, percentileUpper)
  
  if(!is.null(target)){
    isTargetWithinPercentileCI <- (percentileLower < target) & (target < percentileUpper)
    ret <- data.frame(percentileLower, percentileUpper, isTargetWithinPercentileCI)
  }
  
  return(ret)
}

# ret
# obs <- zd$GeneratedCorrelation
# target <- zd$TargetCorrelation[1]
# alpha <- 0.05
# calcPercentileCI(obs, target, alpha)
