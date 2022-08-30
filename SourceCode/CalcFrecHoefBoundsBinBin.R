# Calculates the Frechet-Hoeffding bounds on correlation for two binomial r.v.
# px = success rate for first binomial r.v.
# nx = number of trials for first binomial r.v.
# py = success rate for second binomial r.v.
# ny = number of trials for second binomial r.v.
# Outputs the lower bounds, upper bounds and the pmf's that lead to each bound

CalcFrecHoefBoundsBinBin <- function(px = 0.1, nx = as.integer(2), 
                                     py = 0.1, ny = as.integer(2)){
 
  # Check basic conditions
  if (!is.double(px) | px < 0 | 1 < px) {
    stop("Error In initial values, px;
         should be double between 0 and 1 (inclusive)")
  } 
  if (!is.double(py) | py < 0 | 1 < py) {
    stop("Error In initial values, py; 
         should be double between 0 and 1 (inclusive)")
  } 
  if (!is.integer(nx) | nx < 0 ) {
    stop("Error In initial values, nx; 
         should be non-zero integer")
  } 
  if (!is.integer(ny) | ny < 0 ) {
    stop("Error In initial values, ny; 
         should be non-zero integer")
  } 
  
  # Support of X and Y
  ysup <- c(0:ny)
  xsup <- c(0:nx)
  
  # Matrices to store pmf/cdf values
  lowcdfmat <- matrix(0, 
                   nrow = length(ysup), 
                   ncol = length(xsup))
  lowpmfmat <- matrix(0, 
                   nrow = length(ysup), 
                   ncol = length(xsup))
  highcdfmat <- matrix(0, 
                      nrow = length(ysup), 
                      ncol = length(xsup))
  highpmfmat <- matrix(0, 
                      nrow = length(ysup), 
                      ncol = length(xsup))
  
  # Calculate pmf/cdf
  for (yrow in ysup) {
    for (xcol in  xsup) {
     xval <- xcol
     yval <- yrow
     
     xind <- xcol + 1;
     yind <- yrow + 1
     
     # Calculate cdf/pmf for lower bound
     # Definition of cdf for lower bound
      lowcdfmat[yind, xind] <- 
        max(
          pbinom(xval, nx, px) +
            pbinom(yval, ny, py) -
            1,
          0)
      
      # f(x,y) = F(X <= x, Y <= y) - 
      #   F(X <= x - 1, Y <= y) - 
      #   F(X <= x , Y <= y - 1) + 
      #   F(X <= x - 1, Y <= y - 1) 
      lowcdftop <- ifelse(yind > 1, 
                          lowcdfmat[yind - 1, xind],
                          0)
      lowcdfleft <- ifelse(xind > 1, 
                           lowcdfmat[yind, xind - 1],
                           0)
      lowcdftopleft <- ifelse(yind > 1 & xind > 1, 
                              lowcdfmat[yind - 1, xind - 1],
                              0)
      
      lowpmfmat[yind, xind] <- 
        lowcdfmat[yind, xind] - 
        lowcdfleft -
        lowcdftop +
        lowcdftopleft
      
      # Calculate pmf/cdf of upper bound
      # definition of cdf for upper bound
      highcdfmat[yind, xind] <- 
        min(
          pbinom(xval, nx, px),
          pbinom(yval, ny, py)
          )
      
      # f(x,y) = F(X <= x, Y <= y) - 
      #   F(X <= x - 1, Y <= y) - 
      #   F(X <= x , Y <= y - 1) + 
      #   F(X <= x - 1, Y <= y - 1) 
      highcdftop <- ifelse(yind > 1, 
                          highcdfmat[yind - 1, xind],
                          0)
      highcdfleft <- ifelse(xind > 1, 
                           highcdfmat[yind, xind - 1],
                           0)
      highcdftopleft <- ifelse(yind > 1 & xind > 1, 
                              highcdfmat[yind - 1, xind - 1],
                              0)
      
      highpmfmat[yind, xind] <- 
        highcdfmat[yind, xind] - 
        highcdfleft -
        highcdftop +
        highcdftopleft
      
    }
  }
  
  # Calculate the correlation
  # E[X], Var[X]; E[Y], Var[Y]
  eX <- nx*px; varX <- nx*px*(1-px)
  eY <- ny*py; varY <- ny*py*(1-py)
  
  # E[XY]
  xyvalmat <- ysup %*% t(xsup)
  loweXY <- sum(lowpmfmat * xyvalmat)
  higheXY <- sum(highpmfmat * xyvalmat)

  # Correlation
  upperCor <- (higheXY - eX*eY)/sqrt(varX*varY)
  lowerCor <- (loweXY - eX*eY)/sqrt(varX*varY)
  
  # Return value
  ret <- list(
    LowerBound = lowerCor,
    UpperBound = upperCor,
    pmfLower = lowpmfmat,
    pmfUpper = highpmfmat
  )
  
  return(ret)
}

# rm(list = ls())
# px <- 0.1; nx <- as.integer(1)
# py <- 0.3; ny <- as.integer(1)
# testRet <- CalcFrecHoefBoundsBinBin(px,nx,py,ny)
# testRet
# c(testRet$LowerBound, testRet$UpperBound)
# 
# testRet$UpperBound

# Checks that all of the pmf's created have the correct margins
#####
# nxtest <- c(1:100)
# nytest <- c(1:100)
# pxtest <- seq(0, 1, 0.01)
# pytest <- seq(0, 1, 0.01)
# params <- expand.grid(nxtest,pxtest,nytest,pytest)
# 
# retvec <- double(0)
# for (i in nrow(params)) {
#   # Set up parameters
#   # i <- 2
#   nx <- params[i,1]
#   px <- params[i,2]
#   ny <- params[i,1]
#   py <- params[i,2]
#   
#   # Calculate bounds
#   testRet <- CalcFrecHoefBoundsBinBin(px,nx,py,ny)
#   
#   
#   # Checking Margins match intended
#   s1 <- sapply(c(0:nx), function(x){
#     dbinom(x, nx, px)
#   }) - colSums(testRet$pmfUpper)
#   
#   s2 <- sapply(c(0:nx), function(x){
#     dbinom(x, nx, px)
#   }) - colSums(testRet$pmfLower)
#   
#   s3 <- sapply(c(0:ny), function(y){
#     dbinom(y, ny, py)
#   }) - rowSums(testRet$pmfLower)
#   
#   s4 <- sapply(c(0:ny), function(y){
#     dbinom(y, ny, py)
#   }) - rowSums(testRet$pmfUpper)
#   
#   # Append to return vector
#   retvec <- c(retvec, s1, s2, s3, s4)
# }
# summary(retvec)
# Should all be 0
#####