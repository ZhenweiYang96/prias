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

nDataSets = 366
getNextSeed = function(lastSeed){
  lastSeed + 1
}

#3 types of weibull scales and shapes
weibullScales = c(4,5,6)
weibullShapes = c(1.5,3,4.5)

methods = c("expectedFailureTime", "medianFailureTime","youden", "f1score") 
            #"accuracy", "f1score")

simulatedDsList = vector("list", nDataSets)
lastSeed = 19557
for(i in 317:nDataSets){
  print(paste("******** Started working on Data Set: ", i, "*******"))
  
  lastSeed = getNextSeed(lastSeed)
  repeat{
    set.seed(lastSeed)
    simulatedDsList[[i]] = generateSimLongitudinalData(nSub)
    indices = sample(1:length(weibullScales), size = nSub, replace = T, prob = rep(1/length(weibullScales),length(weibullScales)))
    simulatedDsList[[i]]$weibullShape = weibullShapes[indices]
    simulatedDsList[[i]]$weibullScale = weibullScales[indices]
    
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
  colnames(prias_res) = colnames(simulatedDsList[[i]]$biopsyTimes)
  simulatedDsList[[i]]$biopsyTimes = rbind(simulatedDsList[[i]]$biopsyTimes, prias_res)
  
  jh_res = data.frame(foreach(patientRowNum=1:nrow(simulatedDsList[[i]]$testDs.id),
                              .combine = rbind,
                              .packages =  c("splines", "JMbayes", "coda"))%dopar%{
                                res = c(P_ID=patientRowNum,
                                        methodName="JH",
                                        computeNbAndOffset_JH(dsId = i, patientRowNum=patientRowNum))
                                return(res)
                              })
  colnames(jh_res) = colnames(simulatedDsList[[i]]$biopsyTimes)
  simulatedDsList[[i]]$biopsyTimes = rbind(simulatedDsList[[i]]$biopsyTimes, jh_res)
  
  # mixed_res = data.frame(foreach(patientRowNum=1:nrow(simulatedDsList[[i]]$testDs.id),
  #                                .combine = rbind,
  #                                .packages =  c("splines", "JMbayes", "coda"),
  #                                .export = c("timesPerSubject"))%dopar%{
  #                                  res = c(P_ID=patientRowNum,
  #                                          methodName="MixedF1Score",
  #                                          computeNbAndOffset_Mixed(dsId = i, patientRowNum=patientRowNum,
  #                                                                   1, timesPerSubject-1,
  #                                                                   alternative = "f1score"))
  #                                  return(res)
  #                                })
  # 
  # colnames(mixed_res) = colnames(simulatedDsList[[i]]$biopsyTimes)
  # simulatedDsList[[i]]$biopsyTimes = rbind(simulatedDsList[[i]]$biopsyTimes, mixed_res)
  
  mixed_res2 = data.frame(foreach(patientRowNum=1:nrow(simulatedDsList[[i]]$testDs.id),
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
  
  colnames(mixed_res2) = colnames(simulatedDsList[[i]]$biopsyTimes)
  simulatedDsList[[i]]$biopsyTimes = rbind(simulatedDsList[[i]]$biopsyTimes, mixed_res2)
  
  tEnd = Sys.time()
  print(tEnd-tStart)
  
  stopCluster(ct)
  
  nSchedulingMethods = nrow(simulatedDsList[[i]]$biopsyTimes) / nrow(simulatedDsList[[i]]$testDs.id)
  simulatedDsList[[i]]$biopsyTimes$weibullScale = rep(tail(simulatedDsList[[i]]$weibullScale, nrow(simulatedDsList[[i]]$testDs.id)), nSchedulingMethods)
  simulatedDsList[[i]]$biopsyTimes$weibullShape = rep(tail(simulatedDsList[[i]]$weibullShape, nrow(simulatedDsList[[i]]$testDs.id)), nSchedulingMethods)
  
  temp = list(simulatedDsList[[i]])
  save(temp,file = paste("Rdata/Gleason as event/Sim Study/sc_mixed_sh_mixed/Dt_1/simDs",i,".Rdata", sep=""))
  
  #Save RAM
  rm(temp)
  simulatedDsList[[i]] = NA
}
