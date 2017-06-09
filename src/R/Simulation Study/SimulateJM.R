library(MASS)
library(splines)

nSub <- 1000# number of subjects

source("src/R/common.R")
source("src/R/Simulation Study/simCommon.R")
source("src/R/Simulation Study/personalizedDynamicCutoff.R")
source("src/R/rocJM_mod.R")

cores = detectCores()

nDataSets = 3
getNextSeed = function(lastSeed){
  lastSeed + 1
}

#2.2 and 6
weibullScales = rep(8, nDataSets)
weibullShapes = rep(4.5, nDataSets)

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
  summary(simulatedDsList[[i]]$testDs.id$progression_time)
  summary(simulatedDsList[[i]]$trainingDs.id$progression_time)
  simulatedDsList[[i]]$models = fitJointModelSimDs(simulatedDsList[[i]]$trainingDs.id, simulatedDsList[[i]]$trainingDs)
  
  print(summary(simulatedDsList[[i]]$models$mvJoint_psa_tdboth_training))
  
  simulatedDsList[[i]]$dynamicCutOffTimes = seq(generateLongtiudinalTimeBySchedule()[2], max(simulatedDsList[[i]]$trainingDs$progression_time)-0.0001, 0.1)
  
  simulatedDsList[[i]]$rocList = computeRoc(i)
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
  save(temp,file = paste("Rdata/Gleason as event/Sim Study/sc_8_sh_4pt5/simDs",i,".Rdata", sep=""))
  
  #Save RAM
  simulatedDsList[[i]] = NA
}
