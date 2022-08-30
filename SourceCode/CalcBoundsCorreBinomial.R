# Calculates the maximum rho based on optimizing a set of linear equations
# Linear equations derived from the definition of correlation given bivariate frequencies
# And the binomial marginals holding. Uses the table given by Binomial variables


getRhoRangeBinOpt <- function(p1,p2){
  # Input: 2 minor allele frequencies
  # output: The lower and upper bound in a named vector ("Lower", "Upper")
  
  require(lpSolve)
  
  # Testing values
  # Minor allele frequencies
  # p1 <- 0.01
  # p2 <- 0.015
  
  # Format of 00, 01, 02, 10, 11, 12, 20, 21, 22
  # AA = 0, Aa = 1, aa = 2 where 0,1,2 are form the binomial(n, p)
  ## Set the coefficients of the decision variables -> C
  C <- c(0, 0, 0, 0, 1, 2, 0, 2, 4)
  # C
  
  # Create constraint martix 
  A <- matrix(c(1, 0, 0, 1, 0, 0, 1, 0, 0,
                0, 1, 0, 0, 1, 0, 0, 1, 0,
                0, 0, 1, 0, 0, 1, 0, 0, 1,
                0, 0, 0, 1, 1, 1, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 1, 1, 1),
              nrow = 5, byrow = T)

  # Right hand side for the constraints
  B <- c((1 - p1)^2, 
         2*(1 - p1)*p1, 
         p1^2,
         2*(1 - p2)*p2, 
         p2^2)

  # Direction of the constraints
  constranints_direction  <- c("=", "=", "=", "=", "=")

  # Find the optimal solution for maximum
  optimum <-  lp(direction = "max",
                 objective.in = C,
                 const.mat = A,
                 const.dir = constranints_direction,
                 const.rhs = B,
                 all.int = F)
  maxCorre <- (optimum$objval - 2*p1*2*p2)/(sqrt(2*p1*(1-p1)*2*p2*(1-p2)))
  Ur <- maxCorre
  
  # Find optimal solution for minimum 
  optimum <-  lp(direction = "min",
                 objective.in = C,
                 const.mat = A,
                 const.dir = constranints_direction,
                 const.rhs = B,
                 all.int = F)
  minCorre <- (optimum$objval - 2*p1*2*p2)/(sqrt(2*p1*(1-p1)*2*p2*(1-p2)))
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
  # colSums(bivarFreq)
  # c((1-p2)^2, 2*p2*(1-p2), p2^2)
  
  lrUr <- c(Lr, Ur)
  names(lrUr) <- c("Lower", "Upper")
  return(lrUr)
}


# getRhoRangeBinOpt(0.1, 0.15)
# getRhoRangeBinOpt(0.1, 0.5)

