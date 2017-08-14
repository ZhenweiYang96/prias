library(MASS)
library(splines)

nSub <- 1000# number of subjects

source("src/R/common.R")
source("src/R/Sim Study Mixture/simCommon.R")
source("src/R/Sim Study Mixture/rocAndCutoff.R")
source("src/R/Sim Study Mixture/nbAndOffset.R")
source("src/R/Sim Study Mixture/produceResults.R")
source("src/R/rocJM_mod.R")

cores = detectCores()

datasetNumbers = 338:340
simulatedDsList = vector("list", max(datasetNumbers))

for(i in datasetNumbers){
  load(paste("Rdata/Gleason as event/Sim Study/sc_mixed_sh_mixed/Dt_1/simDs",i,".Rdata", sep=""))
  simulatedDsList[[i]] = temp[[1]]
  rm(temp)
  
  print(paste("******** Started working on Data Set: ", i, "*******"))
  
  ##############  THE COMMENTED PART IS FOR THE DYN_DT ################
  #simulatedDsList[[i]]$testDs = simulatedDsList[[i]]$testDs[,1:6]
  #simulatedDsList[[i]]$dynamicCutOffTimes = seq(generateLongtiudinalTimeBySchedule()[2], max(simulatedDsList[[i]]$trainingDs$progression_time)-0.0001, 0.1)
  
  #simulatedDsList[[i]]$rocList = computeRocDataDrivenDt(i)
  #simulatedDsList[[i]]$cutoffValues = computeCutOffValues(i)
  
  #simulatedDsList[[i]]$testDs$expectedFailureTime = NA
  #simulatedDsList[[i]]$testDs$survTime85 = NA
  #simulatedDsList[[i]]$testDs$survTimeYouden = NA
  #simulatedDsList[[i]]$testDs$survTimeMaxTPR = NA
  #simulatedDsList[[i]]$testDs$survTimeAccuracy = NA
  #simulatedDsList[[i]]$testDs$survTimeF1Score = NA
  #####################################################################

  ct= makeCluster(cores)
  registerDoParallel(ct)
  
  tStart = Sys.time()
  lastPossibleVisit = timesPerSubject - 1
  
  #Combo of patientRowNum and Method
  mixed_res = data.frame(foreach(patientRowNum=1:nrow(simulatedDsList[[i]]$testDs.id),
                                 .combine = rbind,
                                 .packages =  c("splines", "JMbayes", "coda"),
                                 .export = c("timesPerSubject"))%dopar%{
                                   res = c(P_ID=patientRowNum,
                                           methodName="MixedF1Score",
                                           computeNbAndOffset_Mixed(dsId = i, patientRowNum=patientRowNum,
                                                                    1, timesPerSubject-1,
                                                                    alternative = "f1score"))
                                   return(res)
                                 })
  mixed_res[,"weibullScale"] = tail(simulatedDsList[[i]]$weibullScale, nrow(simulatedDsList[[i]]$testDs.id))
  mixed_res[,"weibullShape"] = tail(simulatedDsList[[i]]$weibullShape, nrow(simulatedDsList[[i]]$testDs.id))
  
  colnames(mixed_res) = colnames(simulatedDsList[[i]]$biopsyTimes)
  simulatedDsList[[i]]$biopsyTimes = rbind(simulatedDsList[[i]]$biopsyTimes, mixed_res)
  
  tEnd = Sys.time()
  print(tEnd-tStart)
  
  stopCluster(ct)
  
  temp = list(simulatedDsList[[i]])
  save(temp,file = paste("Rdata/Gleason as event/Sim Study/sc_mixed_sh_mixed/Dt_1/simDs",i,".Rdata", sep=""))
  
  print(paste("******** End working on Data Set: ", i, "*******"))
  
  #Save RAM
  rm(temp)
  simulatedDsList[[i]] = NA
}
