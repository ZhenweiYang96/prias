source("src/R/common.R")

library(MASS)
library(splines)



nSub <- 50 # number of subjects

nDataSets = 2

seeds = seq(7001, 7000 + nDataSets)
weibullScales = rep(1.5, nDataSets)
weibullShapes = rep(2.8, nDataSets)

