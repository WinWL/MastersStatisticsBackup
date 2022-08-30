# Input: binMat: 2 column matrix, entries are 0,1,2
# Input: p1A: probability of minor allele
# Input: p2A: probability of minor allele

# Output: single row df whose columns are the sample 
# and expected frequncies of bivariate binomial 0,1,2 values
# (0,0), (0,1), (0,2), 
# (1,0), (1,1), (1,2)
# (2,0), (2,1), (2,2)

extractSnpBivarDfVec <- function(binMat, p1A, p2A) {
  require(dplyr)
  
  # Bivariate Snp Totals of generated Sample
  b <- binMat %>% as.data.frame()
  btotalnum <- nrow(b)
  colnames(b) <- c("S1", "S2")
  
  
  Freq00Samp <- b %>% filter(S1 == 0 & S2 == 0) %>% nrow()
  Freq01Samp <- b %>% filter(S1 == 0 & S2 == 1) %>% nrow()
  Freq02Samp <- b %>% filter(S1 == 0 & S2 == 2) %>% nrow()
  
  Freq10Samp <- b %>% filter(S1 == 1 & S2 == 0) %>% nrow()
  Freq11Samp <- b %>% filter(S1 == 1 & S2 == 1) %>% nrow()
  Freq12Samp <- b %>% filter(S1 == 1 & S2 == 2) %>% nrow()
  
  Freq20Samp <- b %>% filter(S1 == 2 & S2 == 0) %>% nrow()
  Freq21Samp <- b %>% filter(S1 == 2 & S2 == 1) %>% nrow()
  Freq22Samp <- b %>% filter(S1 == 2 & S2 == 2) %>% nrow()
  
  # Converting to total to frequenceis
  Freq00Samp <- Freq00Samp/btotalnum
  Freq01Samp <- Freq01Samp/btotalnum
  Freq02Samp <- Freq02Samp/btotalnum
  
  Freq10Samp <- Freq10Samp/btotalnum
  Freq11Samp <- Freq11Samp/btotalnum
  Freq12Samp <- Freq12Samp/btotalnum
  
  Freq20Samp <- Freq20Samp/btotalnum
  Freq21Samp <- Freq21Samp/btotalnum
  Freq22Samp <- Freq22Samp/btotalnum
  
  
  # Expected SNP frequencies with 0 correlation
  # SNP1
  b1s0 <- (1-p1A)^2
  b1s1 <- 2*p1A*(1-p1A)
  b1s2 <- p1A^2
  
  # SNP2
  b2s0 <- (1-p2A)^2
  b2s1 <- 2*p2A*(1-p2A)
  b2s2 <- p2A^2
  
  # Expected bivaraite SNP frequencies
  Freq00Expe <- b1s0*b2s0
  Freq01Expe <- b1s0*b2s1
  Freq02Expe <- b1s0*b2s2
  
  Freq10Expe <- b1s1*b2s0
  Freq11Expe <- b1s1*b2s1
  Freq12Expe <- b1s1*b2s2
  
  Freq20Expe <- b1s2*b2s0
  Freq21Expe <- b1s2*b2s1
  Freq22Expe <- b1s2*b2s2
  
  singleRowFreq <- cbind(Freq00Samp, Freq01Samp, Freq02Samp,
                         Freq10Samp, Freq11Samp, Freq12Samp,
                         Freq20Samp, Freq21Samp, Freq22Samp,
                         Freq00Expe, Freq01Expe, Freq02Expe,
                         Freq10Expe, Freq11Expe, Freq12Expe,
                         Freq20Expe, Freq21Expe, Freq22Expe
  )
  
  return(singleRowFreq)
}

# p1A <- 0.1
# p2A <- 0.15
# bin1 <- rbinom(n = 1e8, size = 2, prob = p1A)
# bin2 <- rbinom(n = 1e8, size = 2, prob = p2A)
# df <- extractSnpBivarDfVec(cbind(bin1, bin2), p1A, p2A)
# df %>% t
# df

# z1 <- binMat[,1] %>% head(1e8)
# z2 <- binMat[,2] %>% head(1e8)
# 
# Paste all columns together with no separation
# zp <- paste(z1,z2,sep="")
# 
# Convert back to integer and use table to count
# zp %>% as.integer() %>% table()/1e8
# Rename table headers using the naming system for the pmf
# Easy way to get the multivariate pmf for more than 1 snp

