
library(rmetasim)
library(dplyr)

lndScp <- landscape.new.example()
lndScp %>% View



exampleland <- landscape.simulate(lndScp, 4)

exampleland %>% View

exampleland$individuals %>% head()
exampleland

theta.s.mat <- landscape.theta.s(exampleland)
theta.s.mat


exampleland <- landscape.new.example()
exampleland <- landscape.simulate(exampleland, 4)
print("Allele frequencies at locus 1")
table(landscape.states(exampleland,1)[,c(-1:-landscape.democol())])
