rm(list = ls())
library(coala)
model <- coal_model(sample_size = 3, loci_number = 1)

model

model <- model + feat_mutation(rate = 1, model = "IFS")
model

?feat_mutation
sumstats <- simulate(model, seed = 123)
sumstats
