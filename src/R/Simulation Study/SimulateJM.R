library(MASS)
library(splines)

nSub <- 1000 # number of subjects

nDataSets = 1

seeds = seq(7001, 7000 + nDataSets)
weibullScales = rep(8.0, nDataSets)
weibullShapes = rep(5.5, nDataSets)

source("src/R/common.R")
source("src/R/Simulation Study/simCommon.R")
source("src/R/Simulation Study/personalizedDynamicCutoff.R")

simulatedDsList = vector("list", nDataSets)
for(i in 1:nDataSets){
  set.seed(seeds[i])
  simulatedDsList[[i]] = generateSimLongitudinalData(nSub)
  simulatedDsList[[i]]$weibullShape = weibullShapes[i]
  simulatedDsList[[i]]$weibullScale = weibullScales[i]
  simulatedDsList[[i]] = generateSimJointData(i, nSub)
  
  #ggplot(simulatedDsList[[i]]$simDs.id, aes(x=progression_time)) + geom_density()
  #ggplot(data=simulatedDsList[[1]]$trainingDs, aes(x=visitTimeYears, y=logpsa1, group=P_ID)) + geom_line()
  
  simulatedDsList[[i]]$models = fitJointModelSimDs(simulatedDsList[[i]]$trainingDs.id, simulatedDsList[[i]]$trainingDs)
  
  registerDoParallel(cores = 4)
  simulatedDsList[[i]]$cutoffValues = computeCutOffValues(i)
  simulatedDsList[[i]]$biopsyTimes = computeBiopsyTimes(minVisits = 5, dsId = i)
  
  temp = list(simulatedDsList[[i]])
  save(temp,file = paste("Rdata/Gleason as event/Sim Study/simDs",i,".Rdata", sep=""))
}
