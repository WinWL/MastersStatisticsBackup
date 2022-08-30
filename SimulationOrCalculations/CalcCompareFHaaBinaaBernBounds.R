# Constructs a dataframe to compare the correlation
# Frechet Hoeffding bounds,
# THe bounds calculated by optimizing bernoulli pmf
# And the bounds calculated by optimizing Binomial PMf

library(dplyr)
library(lpSolve)
library(ggplot2)

# rm(list=ls())
# setwd("F:/Dropbox/Dropbox/School/SS_Thes/Exp/SimulationOrCalculations")
# setwd("D:/Dropbox/Dropbox/School/SS_Thes/Exp/SimulationOrCalculations")
source("../SourceCode/CalcBoundsCorreBernoulli.R")
source("../SourceCode/CalcBoundsCorreBinomial.R")
source("../SourceCode/CalcFrecHoefBoundsBinBin.R")

# plist <- seq(0.01, 0.1, 0.01)
# param <- expand.grid(
#   plist, plist
#   ) %>% filter(Var1 > Var2)
# 
# row <- param[1,]
# row
# dfBounds <- apply(param, 1, function(row){
#   px <- as.double(row[1])
#   py <- as.double(row[2])
#   
#   nbin <- as.integer(2)
#   nber <- as.integer(1)
#   
#   FHBin <- CalcFrecHoefBoundsBinBin(px, nbin, py, nbin)
#   FHBer <- CalcFrecHoefBoundsBinBin(px, nber, py, nber)
#   
#   BinOpt <- getRhoRangeBinOpt(px, py)
#   BernOpt <- getRhoRangeBernOpt(px, py)
#   
#   FHBinBound <- c(FHBin$LowerBound, FHBin$UpperBound)
#   names(FHBinBound) <- c("BinFrHLower", "BinFrHUpper")
#   
#   FHBerBound <- c(FHBer$LowerBound, FHBer$UpperBound)
#   names(FHBerBound) <- c("BerFrHLower", "BerFrHUpper")
#   
#   BinBound <- BinOpt
#   names(BinBound) <- c("BinOptLower", "BinOptUpper")
#   
#   BernBound <- BernOpt
#   names(BernBound) <- c("BerOptLower", "BerOptUpper")
#   
#   pa <- c(px, py)
#   names(pa) <- c("p1", "p2")
#   
#   ret <- c(pa, FHBinBound, BinBound, BernBound, FHBerBound)
#   ret
#   # ret %>% as.matrix() %>% t
#   return(ret)  
# }) %>% t() %>% as.data.frame()
# 
# dfBounds %>% head()
# 
# # Looking at the data
# dfUpperOnly <- dfBounds[-c(4,6,8)]
# dfUpperOnly
# 
# dfLowerOnly <- dfBounds[-c(5,7,9)]
# dfLowerOnly
# 
# 
# (dfBounds$FrHLower - dfBounds$BinLower) %>% summary
# (dfBounds$FrHUpper - dfBounds$BinUpper) %>% summary
# 
# (dfBounds$FrHLower - dfBounds$BerLower) %>% summary
# (dfBounds$FrHUpper - dfBounds$BerUpper) %>% summary
# dfBounds %>% colnames
# # Constructing boxplot
# library(reshape2)
# library(ggplot2)
# dfPlot <- data.frame(
#   BinFHLowerMinusBinOptLower =
#     (dfBounds$BinOptLower - dfBounds$BinOptLower),
#   BinFHUpperMinusBinOptUpper =
#     (dfBounds$BinFrHUpper - dfBounds$BinOptUpper),
#   BinFHLowerMinusBerOptLower =
#     (dfBounds$BinFrHLower - dfBounds$BerOptLower),
#   BinFHUpperMinusBerOptUpper =
#     (dfBounds$BinFrHUpper - dfBounds$BerOptUpper),
#   BerFHLowerMinusBinOptLower =
#     (dfBounds$BerFrHLower - dfBounds$BinOptLower),
#   BerFHUpperMinusBinOptUpper =
#     (dfBounds$BerFrHUpper - dfBounds$BinOptUpper),
#   BerFHLowerMinusBerOptLower =
#     (dfBounds$BerFrHLower - dfBounds$BerOptLower),
#   BerFHUpperMinusBerOptUpper =
#     (dfBounds$BerFrHUpper - dfBounds$BerOptUpper)
# ) %>% melt
# dfPlot

# pl <- 
#   ggplot(dfPlot, 
#          aes(x = variable, y = value, fill = variable)) +
#   geom_boxplot() +
#   theme(legend.position = "none") +
#   ggtitle("Difference between Frechet Hoeffding, Binomial and Bernoulli Opt Bounds") +
#   xlab("") +
#   ylab("Difference")
# pl
# ggsave("../Results/summimages/CorrelationBoundDiff.png",
       # plot = pl,
       # width = 7, height = 5)


# Caption: COrrelation bounds were calculated on all combination of two allele frequencies from 0.01 to 0.99 (increments of 0.01). The correlations were calculated using four method: The Frechet Hoeffding bounds between bernoulli random variables, the Frechet Hoeffding bounds between binomial random variables, optimizing the correlation from the binomial pmf and optimizing the correlation from the bernoulli pmf. The boxplots of the differences are displayed above. 
# 
# The bounds on optimization are the same as the Frechet Hoeffding bounds when the underlying random variables (both binomial or both bernoulli) are the same. Effectively, by optimizing the correlation by solving a set of equations involving the probability mass functions, we are getting the pmfs that also lead to the bounds on correlation. This is effectively how the Frechet Hoeffding bounds are calculated although they are much more elegant.
# 
# Note that the bernoulli optimization bounds originate from one of the methods of creating correlated binomials. This suggests that the method of adding correlation bernoullis are bounded by the FH bounds between bernoulli random variables rather than the FH bounds between binomial random variables. 

# Check the bernoulli method agaisnt the bounds derived for it via optimiation. 


#p1 <- c(0.1, 0.3)
p1 <- c(0.25)
p2 <- seq(0.01, 0.5, 0.02)
nbin <- as.integer(2)

# px <- 0.1; py <- 0.1
# FHBin <- CalcFrecHoefBoundsBinBin(px, nbin, py, nbin)
p12 <- expand.grid(p1,p2)
# p12
# p12r <- p12[1,]
# p12r

dfPlotBoundBin <- apply(p12, 1, function(p12r){
  px <- as.double(p12r[1]); py <- as.double(p12r[2]); nbin <- as.integer(2)
  FHBin <- CalcFrecHoefBoundsBinBin(px, nbin, py, nbin)
  
  ret <- c(px,py, FHBin$LowerBound, FHBin$UpperBound)
  names(ret) <- c("p1", "p2", "LowerBound", "UpperBound")
  ret
}) %>% t() %>% data.frame()

dfPlotBoundBin %>% head

dfPlotBoundBin$LabelNames <- ifelse(dfPlotBoundBin$p1 == p1[1], "(a)", "(b)" )
dfPlotBoundBin$LabelNames %>% unique

dfPlotBoundBin$ppp <- p12[,2]/p12[,1]

dfPlotBoundBin$VertLine <- ifelse(dfPlotBoundBin$p1 == p1[1], p1[1], p1[2])
dfPlotBoundBin$VertLine %>% unique


dfPlotBoundBern<- apply(p12, 1, function(p12r){
  px <- as.double(p12r[1]); py <- as.double(p12r[2]); nbin <- as.integer(1)
  FHBin <- CalcFrecHoefBoundsBinBin(px, nbin, py, nbin)
  
  ret <- c(px,py, FHBin$LowerBound, FHBin$UpperBound)
  names(ret) <- c("p1", "p2", "LowerBound", "UpperBound")
  ret
}) %>% t() %>% data.frame()

dfPlotBoundBern %>% head

dfPlotBoundBern$LabelNames <- ifelse(dfPlotBoundBern$p1 == p1[1], "(a)", "(b)" )
dfPlotBoundBern$LabelNames %>% unique

dfPlotBoundBern$ppp <- p12[,2]/p12[,1]

dfPlotBoundBern$VertLine <- ifelse(dfPlotBoundBern$p1 == p1[1], p1[1], p1[2])
dfPlotBoundBern$VertLine %>% unique

dfPlotBoundBin$BB <- "Bin"
dfPlotBoundBern$BB <- "Bern"
dfPlotBoundBin %>% head()
dfPlotBoundBern %>% head()

dfPlotBound <- rbind(dfPlotBoundBin, dfPlotBoundBern)
pl <- ggplot(dfPlotBound, aes(x = p2, y = UpperBound)) + 
  # geom_point(aes(shape = BB), colour="gray70", size = 5) +
  # geom_point(aes(color = BB, shape = BB), size = 3) +
  geom_point(aes(color = BB, shape = BB), size = 6) +
  facet_wrap( ~ LabelNames, scales = "free") +
  # ggtitle("Upper Bound on Correlation Between Binomial Random Variables") +
  scale_x_continuous(name ="P2", 
                     breaks = seq(0,0.5,0.1),
                     limits = c(0, 0.51)) +
  scale_y_continuous(name ="Correlation", 
                     breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
                     limits = c(0, 1.05)) +
  geom_vline(aes(xintercept = VertLine), linetype = "dashed", alpha = 0.5) +
  theme(text = element_text(size = 15),
        legend.position = "none")
pl

# ggsave("../Results/summimages/ttCh2FHCorBound2.pdf", plot = pl, width = 11, height = 6)
# ggsave("../PaperFigs/ttCh2FHCorBound2.pdf", plot = pl, width = 11, height = 6)


pl <- ggplot(dfPlotBound, aes(x = p2, y = LowerBound)) + 
  geom_point() +
  facet_wrap( ~ LabelNames, scales = "free") +
  ggtitle("Upper Bound on Correlation Between Binomial Random Variables") 


pl
