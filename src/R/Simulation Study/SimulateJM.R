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

nDataSets = 10
getNextSeed = function(lastSeed){
  lastSeed + 1
}

#2.2 and 6
weibullScales = rep(8, nDataSets)
weibullShapes = rep(4.5, nDataSets)

methods = c("expectedFailureTime", "medianFailureTime","youden", 
            "accuracy", "f1score")

simulatedDsList = vector("list", nDataSets)
lastSeed = 3000 + 8
for(i in 9:nDataSets){
  print(paste("******** Started working on Data Set: ", i, "*******"))
  
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
  
  print(paste("Nr. of training:", nrow(simulatedDsList[[i]]$trainingDs.id), 
              "; Nr. of test:", nrow(simulatedDsList[[i]]$testDs.id)))
  
  print("Beginning to fit joint model")
  simulatedDsList[[i]]$models = fitJointModelSimDs(simulatedDsList[[i]]$trainingDs.id, simulatedDsList[[i]]$trainingDs)
  print("Joint model fitted")
  
  print(summary(simulatedDsList[[i]]$models$mvJoint_psa_tdboth_training))
  
  simulatedDsList[[i]]$dynamicCutOffTimes = seq(generateLongtiudinalTimeBySchedule()[2], max(simulatedDsList[[i]]$trainingDs$progression_time)-0.0001, 0.1)
  
  print("Computing ROC")
  simulatedDsList[[i]]$rocList = computeRoc(i, Dt = 1)
  #simulatedDsList[[i]]$rocList = computeRocDataDrivenDt(i)
  print("Done computing ROC, now computing cutoff values")
  simulatedDsList[[i]]$cutoffValues = computeCutOffValues(i)

  ct= makeCluster(cores)
  registerDoParallel(ct)
  
  tStart = Sys.time()
  lastPossibleVisit = timesPerSubject - 1
  
  #Combo of patientRowNum and Method
  nTasks = length(methods) * nrow(simulatedDsList[[i]]$testDs.id)
  simulatedDsList[[i]]$biopsyTimes = foreach(taskNumber=1:nTasks, .combine = rbind,
                                             .packages =  c("splines", "JMbayes", "coda"),
                                             .export = c("timesPerSubject"))%dopar%{
                                                patientRowNum = ceiling(taskNumber/length(methods))
                                                methodName = taskNumber%%length(methods)
                                                if(methodName == 0){
                                                  methodName = methods[length(methods)]
                                                }else{
                                                  methodName = methods[methodName]
                                                }

                                                res = c(P_ID=patientRowNum, methodName=methodName, computeNbAndOffset(dsId = i, patientRowNum=patientRowNum,
                                                                                                                      minVisits = 1,
                                                                                                                      methodName = methodName,
                                                                                                                      lastPossibleVisit = lastPossibleVisit))

                                                return(res)
                                             }

  
  prias_res = data.frame(foreach(patientRowNum=1:nrow(simulatedDsList[[i]]$testDs.id),
                                 .combine = rbind,
                                 .packages =  c("splines", "JMbayes", "coda"))%dopar%{
                                   res = c(P_ID=patientRowNum,
                                           methodName="PRIAS",
                                           computeNbAndOffset_PRIAS(dsId = i, patientRowNum=patientRowNum))
                                   return(res)
                                 })
  simulatedDsList[[i]]$biopsyTimes = rbind(simulatedDsList[[i]]$biopsyTimes, prias_res)
  
  jh_res = data.frame(foreach(patientRowNum=1:nrow(simulatedDsList[[i]]$testDs.id),
                              .combine = rbind,
                              .packages =  c("splines", "JMbayes", "coda"))%dopar%{
                                res = c(P_ID=patientRowNum,
                                        methodName="JH",
                                        computeNbAndOffset_JH(dsId = i, patientRowNum=patientRowNum))
                                return(res)
                              })
  simulatedDsList[[i]]$biopsyTimes = rbind(simulatedDsList[[i]]$biopsyTimes, jh_res)
  
  # mixed_res = data.frame(foreach(patientRowNum=1:nrow(simulatedDsList[[i]]$testDs.id), 
  #                                .combine = rbind, 
  #                                .packages =  c("splines", "JMbayes", "coda"),
  #                                .export = c("timesPerSubject", "dynamicCutOffTimes"))%dopar%{
  #                                  res = c(P_ID=patientRowNum, 
  #                                          methodName="MixedAccu", 
  #                                          computeNbAndOffset_Mixed(dsId = i, patientRowNum=patientRowNum,
  #                                                                   1, timesPerSubject-1, 
  #                                                                   alternative = "accuracy"))
  #                                  return(res)
  #                                })
  # 
  # 
  # simulatedDsList[[i]]$biopsyTimes = rbind(simulatedDsList[[i]]$biopsyTimes, mixed_res)
  
  tEnd = Sys.time()
  print(tEnd-tStart)

  stopCluster(ct)
   
  temp = list(simulatedDsList[[i]])
  save(temp,file = paste("Rdata/Gleason as event/Sim Study/sc_8_sh_4pt5/Dt_1/simDs",i,".Rdata", sep=""))
  
  #Save RAM
  rm(temp)
  simulatedDsList[[i]] = NA
}
