library(MASS)
library(splines)

nSub <- 1000# number of subjects

source("src/R/common.R")
source("src/R/Simulation Study/simCommon.R")
source("src/R/Simulation Study/personalizedDynamicCutoff.R")
source("src/R/Simulation Study/rocJM_mod.R")

cores = detectCores()

nDataSets = 10
getNextSeed = function(lastSeed){
  lastSeed + 1
}

#2.2 and 6
weibullScales = rep(6, nDataSets)
weibullShapes = rep(2.5, nDataSets)

simulatedDsList = vector("list", nDataSets)
lastSeed = 4000
for(i in 1:8){
  load(paste("Rdata/Gleason as event/Sim Study/sc_6_sh2pt5/simDs",i,".Rdata", sep=""))
  simulatedDsList[[i]] = temp[[1]]
  rm(temp)
  
  simulatedDsList[[i]]$testDs = simulatedDsList[[i]]$testDs[,1:6]
  simulatedDsList[[i]]$dynamicCutOffTimes = seq(generateLongtiudinalTimeBySchedule()[2], max(simulatedDsList[[i]]$trainingDs$progression_time)-0.0001, 0.1)
  
  simulatedDsList[[i]]$rocList = computeRocDataDrivenDt(i)
  simulatedDsList[[i]]$cutoffValues = computeCutOffValues(i)

  #make new columns
  #This order should match the order of result from compute biopsy times
  simulatedDsList[[i]]$testDs$expectedFailureTime = NA
  simulatedDsList[[i]]$testDs$survTime85 = NA
  simulatedDsList[[i]]$testDs$survTimeYouden = NA
  simulatedDsList[[i]]$testDs$survTimeMaxTPR = NA
  simulatedDsList[[i]]$testDs$survTimeAccuracy = NA
  simulatedDsList[[i]]$testDs$survTimeF1Score = NA

  resultColNumbers = (ncol(simulatedDsList[[i]]$testDs) - 6 + 1):ncol(simulatedDsList[[i]]$testDs)

  ct= makeCluster(cores)
  registerDoParallel(ct)
  
  tStart = Sys.time()
  lastPossibleVisit = timesPerSubject - 1
  
  for(visitNumber in 1:lastPossibleVisit){
    biopsyTimesList = computeBiopsyTimes(i, visitNumber)
    for(j in 1:nrow(simulatedDsList[[i]]$testDs.id)){
      simulatedDsList[[i]]$testDs[(j-1)*timesPerSubject + visitNumber, resultColNumbers] = biopsyTimesList[[j]]  
    }
    print(paste("visitNumber", visitNumber))
  }
  
  # simulatedDsList[[i]]$biopsyTimes = vector("list", nrow(simulatedDsList[[i]]$testDs.id))
  # for(patientRowNum in 1:length(simulatedDsList[[i]]$biopsyTimes)){
  #   simulatedDsList[[i]]$biopsyTimes[[patientRowNum]] = computeBiopsyTimes(minVisits = 1, dsId = i, patientRowNum)
  #   print(paste("Subject", patientRowNum))
  # }
  tEnd = Sys.time()
  print(tEnd-tStart)

  stopCluster(ct)
   
  temp = list(simulatedDsList[[i]])
  save(temp,file = paste("Rdata/Gleason as event/Sim Study/sc_6_sh2pt5/simDs",i,".Rdata", sep=""))
  
  #Save RAM
  simulatedDsList[[i]] = NA
}

#See the results
biopsyResults = getBiopsyResults(2, biopsyIfLessThanTime = 1, biopsyEveryKYears = 3,  minVisits = 1)
biopsyResults = biopsyResults[biopsyResults$methodName %in% c("survTime85", "survTimeAccuracy", "survTimeMaxTPR","survTimeYouden","survTimeF1Score" ,"fixed", "johnsSummary"),]
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

plotBiopsy2DPlot(nbSummary, offsetSummary)
