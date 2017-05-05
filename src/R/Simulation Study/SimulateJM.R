library(MASS)
library(splines)

nSub <- 1000# number of subjects

source("src/R/common.R")
source("src/R/Simulation Study/simCommon.R")
source("src/R/Simulation Study/personalizedDynamicCutoff.R")
source("src/R/Simulation Study/rocJM_mod.R")

nDataSets = 10
getNextSeed = function(lastSeed){
  lastSeed + 1
}

#2.2 and 6
weibullScales = rep(6, nDataSets)
weibullShapes = rep(1.5, nDataSets)

simulatedDsList = vector("list", nDataSets)
lastSeed = 3000
for(i in 1:nDataSets){
  lastSeed = getNextSeed(lastSeed)
  repeat{
    set.seed(lastSeed)
    simulatedDsList[[i]] = generateSimLongitudinalData(nSub)
    simulatedDsList[[i]]$weibullShape = weibullShapes[i]
    simulatedDsList[[i]]$weibullScale = weibullScales[i]
  
    simulatedDsList[[i]] = try(generateSimJointData(i, nSub), T)
    if(inherits(simulatedDsList[[i]], "try-error")){
      print(simulatedDsList[[i]])
      lastSeed = getNextSeed(lastSeed)
    }else{
      simulatedDsList[[i]]$seed = lastSeed
      break
    }
  }
  
  failureTimeCompDs = data.frame(failuretime = c(simulatedDsList[[i]]$trainingDs.id$progression_time,
                                                 prias.id$progression_time),
                                 type = c(rep("Sim",length(simulatedDsList[[i]]$trainingDs.id$progression_time)),
                                          rep("Obs", length(prias.id$progression_time))))

  ggplot(data = failureTimeCompDs) + geom_density(aes(x=failuretime, color=type, fill=type), alpha=0.3)
  
  summary(prias.id$progression_time)
  summary(simulatedDsList[[1]]$testDs.id$progression_time)
  summary(simulatedDsList[[1]]$trainingDs.id$progression_time)
  simulatedDsList[[i]]$models = fitJointModelSimDs(simulatedDsList[[i]]$trainingDs.id, simulatedDsList[[i]]$trainingDs)
  
  print(summary(simulatedDsList[[i]]$models$mvJoint_psa_tdboth_training))
  
  simulatedDsList[[i]]$dynamicCutOffTimes = seq(generateLongtiudinalTimeBySchedule()[2], max(simulatedDsList[[1]]$trainingDs$progression_time)-0.0001, 0.1)
  
  simulatedDsList[[i]]$rocList = computeRoc(i)
  simulatedDsList[[i]]$cutoffValues = computeCutOffValues(i)

  tStart = Sys.time()

  ct= makeCluster(4)
  registerDoParallel(ct)
  simulatedDsList[[i]]$biopsyTimes = vector("list", nrow(simulatedDsList[[i]]$testDs.id))
  for(patientRowNum in 1:length(simulatedDsList[[i]]$biopsyTimes)){
    simulatedDsList[[i]]$biopsyTimes[[patientRowNum]] = computeBiopsyTimes(minVisits = 1, dsId = i, patientRowNum)
    print(paste("Subject", patientRowNum))
  }
  tEnd = Sys.time()
  print(tEnd-tStart)

  stopCluster(ct)
   
  temp = list(simulatedDsList[[i]])
  save(temp,file = paste("Rdata/Gleason as event/Sim Study/sc_6_sh1pt5/simDs",i,".Rdata", sep=""))
}

#See the results
biopsyResults = getBiopsyResults(1, biopsyIfLessThanTime = 1, minVisits = 1)
biopsyResults = biopsyResults[!biopsyResults$methodName %in% c("survTimeMinFPR", "expectedFailureTime"),]
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
