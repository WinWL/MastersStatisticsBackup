
# Given p1 and p2 minor allele frequencies,
# Gives the lower and upper bound for Rho, the correlation coefficient
# Derivation is seen on page 23 of sept 9 book.
# Can be derived from https://stats.stackexchange.com/questions/284996/generating-correlated-binomial-random-variables
# And bounding all frequencies from 0 to 1. 

getRhoRangeBernEqs <- function(p1, p2) {
  z <- sqrt(p1*p2*(1-p1)*(1-p2))
  Lr <- max(
    -(1-p1)*(1-p2),
    -p1*p2)/z
  
  Ur <- min(
    1 - p1 - (1-p1)*(1-p2),
    1 - p2 - (1-p1)*(1-p2))/z
  
  lrUr <- c(Lr, Ur)
  names(lrUr) <- c("Lower", "Upper")
  return(lrUr)
}

# getRhoRangeBern(0.05,0.1)
# getRhoRangeBern(0.1,0.1)
# getRhoRangeBern(0.1,0.15)


# Calculates the maximum rho based on optimizing a set of linear equations
# Linear equations derived from the definition of correlation given bivariate frequencies
# And the binomial marginals holding. Uses the Frequency Table given by Bernoulli variables
getRhoRangeBernOpt <- function(p1,p2){
  require(lpSolve)
  
  # Minor allele frequencies testing values
  # p1 <- 0.01
  # p2 <- 0.015
  # p2 <- 0.99
  # p1 <- 0.59
  
  # Format of 00, 01, 10, 11,
  # A = 0, a = 1,  where 0,1 are form the binomial(1, p)
  ## Set the coefficients of the decision variables -> C
  C <- c(0, 0, 0, 1)

  # Create constraint martix 
  A <- matrix(c(1,0,1,0,
                0,1,0,1,
                0,0,1,1),
              nrow = 3, byrow = T)

  # Right hand side for the constraints
  B <- c(1-p1, p1, p2)

  # Direction of the constraints
  constranints_direction  <- c("=", "=", "=")

  # Find the optimal solution
  optimum <-  lp(direction = "max",
                 objective.in = C,
                 const.mat = A,
                 const.dir = constranints_direction,
                 const.rhs = B,
                 all.int = F)
  maxCorre <- (optimum$objval - p1*p2)/(sqrt(p1*(1-p1)*p2*(1-p2)))
  Ur <- maxCorre
  
  optimum <-  lp(direction = "min",
                 objective.in = C,
                 const.mat = A,
                 const.dir = constranints_direction,
                 const.rhs = B,
                 all.int = F)
  minCorre <- (optimum$objval - p1*p2)/(sqrt(p1*(1-p1)*p2*(1-p2)))
  Lr <- minCorre
  
  # Testing values
  # Print status: 0 = success, 2 = no feasible solution
  # optimum$status
  
  # bivarFreq <- matrix(optimum$solution, nrow=3, ncol=3)
  # maxCorre <- (optimum$objval - 2*p1*2*p2)/(sqrt(2*p1*(1-p1)*2*p2*(1-p2)))
  
  # bivarFreq
  # maxCorre
  # 
  # rowSums(bivarFreq)
  # c((1-p1)^2, 2*p1*(1-p1), p1^2)
  # 
  # 
  # colSums(bivarFreq)
  # c((1-p2)^2, 2*p2*(1-p2), p2^2)
  
  lrUr <- c(Lr, Ur)
  names(lrUr) <- c("Lower", "Upper")
  return(lrUr)
}

# lrUr
# lrUr