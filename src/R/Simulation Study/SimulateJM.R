library(MASS)
library(splines)

nSub <- 1000# number of subjects

nDataSets = 1

seeds = seq(7001, 7000 + nDataSets)
weibullScales = rep(10.9, nDataSets)
weibullShapes = rep(9.5, nDataSets)

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
  range(simulatedDsList[[1]]$testDs.id$progression_time)
  range(simulatedDsList[[1]]$trainingDs.id$progression_time[simulatedDsList[[1]]$trainingDs.id$progressed==1])
  simulatedDsList[[i]]$models = fitJointModelSimDs(simulatedDsList[[i]]$trainingDs.id, simulatedDsList[[i]]$trainingDs)
  
  ct= makeCluster(4)
  registerDoParallel(ct)
  
  simulatedDsList[[i]]$cutoffValues = computeCutOffValues(i)
  tStart = Sys.time()
  simulatedDsList[[i]]$biopsyTimes = vector("list", nrow(simulatedDsList[[i]]$testDs.id))
  for(patientRowNum in 1:nrow(simulatedDsList[[i]]$testDs.id)){
    simulatedDsList[[i]]$biopsyTimes[[patientRowNum]] = computeBiopsyTimes(minVisits = 5, dsId = i, patientRowNum)
    print(paste("Subject", patientRowNum))
  }
  tEnd = Sys.time()
  tEnd-tStart
  
  stopCluster(ct)

  temp = list(simulatedDsList[[i]])
  save(temp,file = paste("Rdata/Gleason as event/Sim Study/simDs",2,".Rdata", sep=""))
}

#See the results
biopsyResults = getBiopsyResults(1)
biopsyResults = biopsyResults[!biopsyResults$methodName %in% "survTimeMaxNPV",]
incompleteRowNum = unique(biopsyResults[is.na(biopsyResults$biopsyTimeOffset) | biopsyResults$biopsyTimeOffset < 0, ]$patientRowNum)
biopsyResultsCC = biopsyResults[!(biopsyResults$patientRowNum %in% incompleteRowNum),]

methodNames = unique(biopsyResultsCC$methodName)

nbSummary = do.call(cbind, by(data = biopsyResultsCC[, "nb"], INDICES = biopsyResultsCC$methodName, summary))
nbSummaryMeanSorted = nbSummary[, names(sort(nbSummary[4,]))]
nbSummaryMedianSorted = nbSummary[, names(sort(nbSummary[3,]))]
write.csv2(nbSummaryMeanSorted, file = "Rdata/Gleason as event/Sim Study/nbSummaryMeanSorted2.csv")
write.csv2(nbSummaryMedianSorted, file = "Rdata/Gleason as event/Sim Study/nbSummaryMedianSorted2.csv")

offsetSummary = do.call(cbind, by(data = biopsyResultsCC[, "biopsyTimeOffset"], INDICES = biopsyResultsCC$methodName, summary))
offsetSummaryMeanSorted = offsetSummary[, names(sort(offsetSummary[4,]))]
offsetSummaryMedianSorted = offsetSummary[, names(sort(offsetSummary[3,]))]

write.csv2(offsetSummaryMeanSorted*12, file = "Rdata/Gleason as event/Sim Study/offsetSummaryMeanSorted2.csv")
write.csv2(offsetSummaryMedianSorted*12, file = "Rdata/Gleason as event/Sim Study/offsetSummaryMedianSorted2.csv")

ggplot(data = biopsyResultsCC) + geom_boxplot(aes(methodName, biopsyTimeOffset*12)) + scale_y_continuous(breaks = 12*seq(0,8, by = 0.25)) +  
  ylab("Biopsy offset (months)")
