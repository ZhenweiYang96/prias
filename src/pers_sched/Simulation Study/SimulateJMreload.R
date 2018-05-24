library(MASS)
library(splines)

nSub <- 1000# number of subjects

source("src/R/common.R")
source("src/R/Simulation Study/simCommon.R")
source("src/R/Simulation Study/rocAndCutoff.R")
source("src/R/Simulation Study/nbAndOffset.R")
source("src/R/Simulation Study/produceResults.R")
source("src/R/rocJM_mod.R")

cores = detectCores()

nDataSets = 1
simulatedDsList = vector("list", nDataSets)

for(i in 1:nDataSets){
  load(paste("Rdata/Gleason as event/Sim Study/sc_8pt5_sh_3/Dt_1/simDs",i,".Rdata", sep=""))
  simulatedDsList[[i]] = temp[[1]]
  rm(temp)
  
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
                                           methodName="MixedYouden",
                                           computeNbAndOffset_Mixed(dsId = i, patientRowNum=patientRowNum,
                                                                    1, timesPerSubject-1,
                                                                    alternative = "youden"))
                                   return(res)
                                 })
  
  colnames(mixed_res) = colnames(simulatedDsList[[i]]$biopsyTimes)
  simulatedDsList[[i]]$biopsyTimes = rbind(simulatedDsList[[i]]$biopsyTimes, mixed_res)
  
  tEnd = Sys.time()
  print(tEnd-tStart)
  
  stopCluster(ct)
  
  temp = list(simulatedDsList[[i]])
  save(temp,file = paste("Rdata/Gleason as event/Sim Study/sc_8pt5_sh_3/Dt_1/simDs",i,".Rdata", sep=""))
  
  #Save RAM
  rm(temp)
  simulatedDsList[[i]] = NA
}
