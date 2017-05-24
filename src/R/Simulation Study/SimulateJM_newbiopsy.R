library(MASS)
library(splines)

nSub <- 1000# number of subjects

source("src/R/common.R")
source("src/R/Simulation Study/simCommon.R")
source("src/R/Simulation Study/personalizedDynamicCutoff.R")
source("src/R/Simulation Study/biopsyTimes.R")
source("src/R/Simulation Study/rocJM_mod.R")

cores = detectCores()

nDataSets = 10
simulatedDsList = vector("list", nDataSets)

methods = c("expectedFailureTime", "survTime85","survTimeYouden", 
           "survTimeAccuracy", "survTimeF1Score", "survTimeMaxTPR")

for(i in 1:nDataSets){
  load(paste("Rdata/Gleason as event/Sim Study/sc_8_sh_4pt5/Dt_1/simDs",i,".Rdata", sep=""))
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
                                                    minVisits = 5, 
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
  save(temp,file = paste("Rdata/Gleason as event/Sim Study/sc_8_sh_4pt5/new/simDs",i,".Rdata", sep=""))
  
  #Save RAM
  simulatedDsList[[i]] = NA
}
