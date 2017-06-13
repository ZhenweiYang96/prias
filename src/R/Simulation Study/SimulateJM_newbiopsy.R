library(MASS)
library(splines)

nSub <- 1000# number of subjects

source("src/R/common.R")
source("src/R/Simulation Study/simCommon.R")
source("src/R/Simulation Study/personalizedDynamicCutoff.R")
source("src/R/Simulation Study/biopsyTimes.R")
source("src/R/rocJM_mod.R")

cores = detectCores()

nDataSets = 10
simulatedDsList = vector("list", nDataSets)

methods = c("expectedFailureTime", "survTime85","survTimeYouden", 
           "survTimeAccuracy", "survTimeF1Score", "survTimeMaxTPR")

for(i in 1:nDataSets){
  load(paste("Rdata/Gleason as event/Sim Study/sc_6_sh_1pt5/Dt_1/simDs",i,".Rdata", sep=""))
  simulatedDsList[[i]] = temp[[1]]
  rm(temp)
  
  simulatedDsList[[i]]$testDs = simulatedDsList[[i]]$testDs[,1:6]
  lastPossibleVisit = timesPerSubject - 1
  
  ct= makeCluster(cores)
  registerDoParallel(ct)
  
  tStart = Sys.time()
  
  #Combo of patientRowNum and Method
  nTasks = length(methods) * nrow(simulatedDsList[[i]]$testDs.id)
  simulatedDsList[[i]]$biopsyTimes = foreach(taskNumber=1:nTasks, .combine = rbind, 
                                            .packages =  c("splines", "JMbayes", "coda"),
                                            .export = c("timesPerSubject", "dynamicCutOffTimes"))%dopar%{
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
    
    #temp = read.csv(csvBackUpFile)
    #temp = rbind(temp, res)
    #write.csv(temp, csvBackUpFile)
    
    #Write it in a file as well
    return(res)
  }
  
  tEnd = Sys.time()
  print(tEnd-tStart)

  stopCluster(ct)
   
  temp = list(simulatedDsList[[i]])
  save(temp,file = paste("Rdata/Gleason as event/Sim Study/sc_6_sh_1pt5/new/simDs",i,".Rdata", sep=""))
  
  #Save RAM
  simulatedDsList[[i]] = NA
}

#########################################
#Adding Johns hopkins and PRIAS and MIXED
#########################################
for(i in 1:10){
  load(paste("Rdata/Gleason as event/Sim Study/sc_6_sh_1pt5/new/simDs",i,".Rdata", sep=""))
  simulatedDsList[[i]] = temp[[1]]
  rm(temp)
  
  lastPossibleVisit = timesPerSubject - 1
  
  ct= makeCluster(cores)
  registerDoParallel(ct)
  
  tStart = Sys.time()
  
  # prias_res = data.frame(foreach(patientRowNum=1:nrow(simulatedDsList[[i]]$testDs.id), 
  #                        .combine = rbind, 
  #                        .packages =  c("splines", "JMbayes", "coda"))%dopar%{
  #                           res = c(P_ID=patientRowNum, 
  #                                   methodName="PRIAS", 
  #                                   computeNbAndOffset_PRIAS(dsId = i, patientRowNum=patientRowNum))
  #                           return(res)
  # })
  # 
  # jh_res = data.frame(foreach(patientRowNum=1:nrow(simulatedDsList[[i]]$testDs.id), 
  #                                .combine = rbind, 
  #                                .packages =  c("splines", "JMbayes", "coda"))%dopar%{
  #                                  res = c(P_ID=patientRowNum, 
  #                                          methodName="JH", 
  #                                          computeNbAndOffset_JH(dsId = i, patientRowNum=patientRowNum))
  #                                  return(res)
  #                                })
  
  mixed_res = data.frame(foreach(patientRowNum=1:nrow(simulatedDsList[[i]]$testDs.id), 
                              .combine = rbind, 
                              .packages =  c("splines", "JMbayes", "coda"),
                              .export = c("timesPerSubject", "dynamicCutOffTimes"))%dopar%{
                                res = c(P_ID=patientRowNum, 
                                        methodName="MixedAccu", 
                                        computeNbAndOffset_Mixed(dsId = i, patientRowNum=patientRowNum,
                                                                1, timesPerSubject-1, 
                                                                alternative = "accuracy"))
                                return(res)
                              })
  
  
  simulatedDsList[[i]]$biopsyTimes = rbind(simulatedDsList[[i]]$biopsyTimes, mixed_res)
  
  tEnd = Sys.time()
  print(tEnd-tStart)
  
  stopCluster(ct)
  
  temp = list(simulatedDsList[[i]])
  save(temp,file = paste("Rdata/Gleason as event/Sim Study/sc_6_sh_1pt5/new/simDs",i,".Rdata", sep=""))
  
  #Save RAM
  simulatedDsList[[i]] = NA
}

