
###################################
# 2 SNP
p1Choice <- c(0.1)
# p2Choice <- c(0.01, 0.1, 0.3)
p2Choice <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5)
p2Choice <- c(0.15)


nBinObs <- 1000 # number of observations to generate in a single sample/generation
nSamCorr <- 1000 # Number of samples/correlations to generate for each target correlation
# tc12Vec <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
tc12Vec <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
tc12Vec <- c(0, 0.2, 0.4, 0.6, 0.8, 1)


# nBinObs <- 10000000 # number of observations to generate in a single sample/generation
# nSamCorr <- 1 # Number of samples/correlations to generate for each target correlation
# # tc12Vec <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
# tc12Vec <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
# tc12Vec <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
# tc12Vec <- c(0.1)


# All Combinations of p1,p2 and target correlations we're interested in
p123TarCor.full <- expand.grid(p1Choice, p2Choice, 
                               tc12Vec
)

# Filter out redundant choices - Unneeded right now
# p123TarCor <- p123TarCor.full[p123TarCor.full$Var1 <= p123TarCor.full$Var2, ]
p123TarCor  <- p123TarCor.full
p123TarCor
p123TarCor %>% dim

pAlleleIndex <- c(1:2)
tarCorIndex <- c(3:3)


# All Combinations of of parameters combined with sample number we're interested in
p123TarCor.full <- expand.grid(p1Choice, p2Choice, 
                               tc12Vec, 
                               c(1:nSamCorr))

# Filter out redundant choices - Unneeded right now
# p123TarCorSim <- p123TarCor.full[p123TarCor.full$Var1 <= p123TarCor.full$Var2, ]
p123TarCorSim <- p123TarCor.full
p123TarCorSim

###################################
# 3 SNP

# p1Choice <- c(0.1)
# p2Choice <- c(0.01, 0.1, 0.3)
# p3Choice <- c(0.01, 0.1, 0.3)
# 
# nBinObs <- 1000 # number of observations to generate in a single sample/generation
# nSamCorr <- 1000 # Number of samples/correlations to generate for each target correlation
# tc12Vec <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
# tc13Vec <- c(0, 0.2, 0.6)
# tc23Vec <- c(0, 0.2, 0.6)
# 486

p1Choice <- c(0.1)
p2Choice <- c(0.1)
p3Choice <- c(0.1)

p1Choice <- c(0.1)
p2Choice <- c(0.15)
p3Choice <- c(0.05)

nBinObs <- 1000 # number of observations to generate in a single sample/generation
nSamCorr <- 10 # Number of samples/correlations to generate for each target correlation
tc12Vec <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
tc13Vec <- c(0, 0.2, 0.8)
tc23Vec <- c(0, 0.2, 0.8)
# tc12Vec <- c(0.1)
# tc13Vec <- c(0.1)
# tc23Vec <- c(0.1)

# All Combinations of p1,p2 and target correlations we're interested in
p123TarCor.full <- expand.grid(p1Choice, p2Choice, 
                               p3Choice,
                               tc12Vec,
                               tc13Vec, tc23Vec
)

# Filter out redundant choices - Unneeded right now
# p123TarCor <- p123TarCor.full[p123TarCor.full$Var1 <= p123TarCor.full$Var2, ]
p123TarCor <- p123TarCor.full
p123TarCor
p123TarCor %>% dim

# Need to split up into 3 so vector isn't too large for R to handle



pAlleleIndex <- c(1:3)
tarCorIndex <- c(4:6)


# All Combinations of of parameters combined with sample number we're interested in
p123TarCor.full <- expand.grid(p1Choice, p2Choice, 
                               p3Choice,
                               tc12Vec, tc13Vec, tc23Vec,
                               c(1:nSamCorr))

# Filter out redundant choices - Unneeded right now
# p123TarCorSim <- p123TarCor.full[p123TarCor.full$Var1 <= p123TarCor.full$Var2, ]
p123TarCorSim <- p123TarCor.full
p123TarCorSim

# Might need to split up so R can handle large object





#############################
# 5 SNP
# Number of possibilities will be too large. Might fix SNp4-SNP5 so the possibilities are only double that of 3 SNPs
p1Choice <- c(0.1)
p2Choice <- c(0.1)
p3Choice <- c(0.1)
p4Choice <- c(0.1)
p5Choice <- c(0.1)

p1Choice <- c(0.1)
p2Choice <- c(0.15)
p3Choice <- c(0.05)
p4Choice <- c(0.3)
p5Choice <- c(0.025)

nBinObs <- 1000 # number of observations to generate in a single sample/generation
nSamCorr <- 1000 # Number of samples/correlations to generate for each target correlation
tc12Vec <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
# tc13Vec <- c(0, 0.2, 0.8)
# tc23Vec <- c(0, 0.8)
tc13Vec <- c(0.2)
tc23Vec <- c(0.2)
tc14Vec <- c(0.2)

tc14Vec <- tc14Vec; tc15Vec <- tc14Vec
tc24Vec <- tc14Vec; tc25Vec <- tc14Vec 
tc34Vec <- tc14Vec; tc35Vec <- tc14Vec
tc45Vec <- tc14Vec

# All Combinations of p1,p2 and target correlations we're interested in
p123TarCor.full <- expand.grid(p1Choice, p2Choice, 
                               p3Choice, p4Choice, p5Choice, 
                               tc12Vec, tc13Vec, tc14Vec, tc15Vec,
                               tc23Vec, tc24Vec, tc25Vec,
                               tc34Vec, tc35Vec, 
                               tc45Vec
)

# Filter out redundant choices - Unneeded right now
# p123TarCor <- p123TarCor.full[p123TarCor.full$Var1 <= p123TarCor.full$Var2, ]
p123TarCor <- p123TarCor.full
p123TarCor %>% head
p123TarCor %>% dim


pAlleleIndex <- c(1:5)
tarCorIndex <- c(6:15)


# All Combinations of of parameters combined with sample number we're interested in
p123TarCor.full <- expand.grid(p1Choice, p2Choice, 
                               p3Choice, p4Choice, p5Choice, 
                               tc12Vec, tc13Vec, tc14Vec, tc15Vec,
                               tc23Vec, tc24Vec, tc25Vec,
                               tc34Vec, tc35Vec, 
                               tc45Vec,
                               c(1:nSamCorr))

# Filter out redundant choices - Unneeded right now
# p123TarCorSim <- p123TarCor.full[p123TarCor.full$Var1 <= p123TarCor.full$Var2, ]
p123TarCorSim <- p123TarCor.full
p123TarCorSim



##############################
# 7 Snp

# Number of possibilities will be too large. Might fix parameters so the possibilities are only double that of 3 SNPs

p1Choice <- c(0.1)
p2Choice <- c(0.01, 0.1, 0.3)
p3Choice <- c(0.01, 0.1, 0.3)
p4Choice <- c(0.01, 0.1, 0.3)
p5Choice <- c(0.01, 0.1, 0.3)
p6Choice <- c(0.01, 0.1, 0.3)
p7Choice <- c(0.01, 0.1, 0.3)

nBinObs <- 1000 # number of observations to generate in a single sample/generation
nSamCorr <- 1000 # Number of samples/correlations to generate for each target correlation
tc12Vec <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
tc13Vec <- c(0, 0.2, 0.6)

# All Combinations of p1,p2 and target correlations we're interested in
# There are 21 pairs with 7 SNPs. 
p123TarCor.full <- expand.grid(p1Choice, p2Choice, 
                               p3Choice, p4Choice, p5Choice, 
                               p6choice, p7choice, 
                               tc12Vec, 
                               tc13Vec, tc13Vec, tc13Vec, tc13Vec, tc13Vec, 
                               tc13Vec, tc13Vec, tc13Vec, tc13Vec, tc13Vec,
                               tc13Vec, tc13Vec, tc13Vec, tc13Vec, tc13Vec, 
                               tc13Vec, tc13Vec, tc13Vec, tc13Vec, tc13Vec
)

# Filter out redundant choices - Unneeded right now
# p123TarCor <- p123TarCor.full[p123TarCor.full$Var1 <= p123TarCor.full$Var2, ]
p123TarCor <- p123TarCor.full
p123TarCor

pAlleleIndex <- c(1:7)
tarCorIndex <- c(8:28)


# All Combinations of of parameters combined with sample number we're interested in
# There are 21 pairs with 7 SNPs. 
p123TarCor.full <- expand.grid(p1Choice, p2Choice, 
                               p3Choice, p4Choice, p5Choice, 
                               p6choice, p7choice, 
                               tc12Vec, 
                               tc13Vec, tc13Vec, tc13Vec, tc13Vec, tc13Vec, 
                               tc13Vec, tc13Vec, tc13Vec, tc13Vec, tc13Vec,
                               tc13Vec, tc13Vec, tc13Vec, tc13Vec, tc13Vec, 
                               tc13Vec, tc13Vec, tc13Vec, tc13Vec, tc13Vec,
                               c(1:nSamCorr))

# Filter out redundant choices - Unneeded right now
# p123TarCorSim <- p123TarCor.full[p123TarCor.full$Var1 <= p123TarCor.full$Var2, ]
p123TarCorSim <- p123TarCor.full
p123TarCorSim

#######################################################
# 3 SNP Fisher Coverage Large Smaple Size

p1Choice <- c(0.1)
p2Choice <- c(0.3)
p3Choice <- c(0.1)

nBinObs <- 10000 # number of observations to generate in a single sample/generation
nSamCorr <- 1000 # Number of samples/correlations to generate for each target correlation
tc12Vec <- c(0.2)
tc13Vec <- c(0.2)
tc23Vec <- c(0)

# All Combinations of p1,p2 and target correlations we're interested in
p123TarCor.full <- expand.grid(p1Choice, p2Choice, 
                               p3Choice,
                               tc12Vec,
                               tc13Vec, tc23Vec
)

# Filter out redundant choices - Unneeded right now
# p123TarCor <- p123TarCor.full[p123TarCor.full$Var1 <= p123TarCor.full$Var2, ]
p123TarCor <- p123TarCor.full
p123TarCor
p123TarCor %>% dim

# Need to split up into 3 so vector isn't too large for R to handle



pAlleleIndex <- c(1:3)
tarCorIndex <- c(4:6)


# All Combinations of of parameters combined with sample number we're interested in
p123TarCor.full <- expand.grid(p1Choice, p2Choice, 
                               p3Choice,
                               tc12Vec, tc13Vec, tc23Vec,
                               c(1:nSamCorr))

# Filter out redundant choices - Unneeded right now
# p123TarCorSim <- p123TarCor.full[p123TarCor.full$Var1 <= p123TarCor.full$Var2, ]
p123TarCorSim <- p123TarCor.full
p123TarCorSim

# Might need to split up so R can handle large object

