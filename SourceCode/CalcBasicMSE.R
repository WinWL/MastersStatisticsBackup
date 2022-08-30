# Calculates a basic mean square error

calcBasicMSE <- function(obsVal,expVal){
  # obsVal: The observed values
  # expVal: The expected value
  
  sum((obsVal - expVal)^2)/length(obsVal)
}

# expVal <- 1;
# obsVal <- rnorm(100, expVal, 3)
# calcBasicMSE(obsVal, expVal)
# 
# # # Comparing to another library that calculates a basic MSE
# # library(DescTools)
# MSE(expVal,obsVal)
