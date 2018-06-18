fitTimeDepCoxModelTrainingDs = function(simDs.id){
  simDs.id = simDs.id[order(simDs.id$progression_time, decreasing = F),]
  
  simDs_time_dep = simDs.id[0,]
  for(i in 1:nrow(simDs.id)){
    intervals = c(0, simDs.id$progression_time[1:i])
    
    truePSA_i = as.numeric(generateTruePSAProfile(simDs.id$Age[i], intervals, simDs.id[i,bNames[3:7]]))
    truePSASlope_i = as.numeric(generateTruePSASlope(intervals, simDs.id[i,bNames[4:7]]))
    trueDRELogOdds_i = as.numeric(generateTrueDRELogOdds(simDs.id$Age[i], intervals, simDs.id[i,bNames[1:2]]))
    
    simDs.id_i = simDs.id[i,]
    simDs.id_i$startTime = -1
    simDs.id_i$stopTime = -1
    simDs.id_i$progressed_in_period = 0
    simDs.id_i$truePSA_in_period = -1
    simDs.id_i$truePSASlope_in_period = -1
    simDs.id_i$trueDRELogOdds_in_period = -1
    
    newDs = simDs.id_i[rep(1, i),]
    for(j in 1:i){
      newDs$startTime[j] = intervals[j]
      newDs$stopTime[j] = intervals[j+1]
      newDs$progressed_in_period[j] = 0
      newDs$truePSA_in_period[j] = truePSA_i[j]
      newDs$truePSASlope_in_period[j] = truePSASlope_i[j]
      newDs$trueDRELogOdds_in_period[j] = trueDRELogOdds_i[j]
    }
    newDs$progressed_in_period[j] = simDs.id_i$progressed
    
    simDs_time_dep = rbind(simDs_time_dep, newDs)
  }
  
  relative_risk_fitted  = coxph(Surv(startTime, stopTime, progressed_in_period) ~ I(Age - 70) +  I((Age - 70)^2) + 
                                  truePSA_in_period + truePSASlope_in_period + trueDRELogOdds_in_period + 
                                  cluster(P_ID), data=simDs_time_dep)
  print(relative_risk_fitted)
  return(relative_risk_fitted)
}

load("Rdata/decision_analytic/DRE_PSA/mvJoint_dre_psa_dre_value_superlight.Rdata")

savedFiles = list.files(path = "D:/", full.names = T)

timeDepCoxModels = vector("list", length(savedFiles))
for(i in 1:length(savedFiles)){
  rm(jointModelData)
  print(paste("Reading the file:", savedFiles[i]))
  load(savedFiles[i])
  
  print("Creating time dependent DS")
  b = jointModelData$mvJoint_dre_psa_simDs$statistics$postMeans$b
  #b = jointModelData$mvglmer_dre_psa_simDs$postMeans$b
  #b = jointModelData$b[1:1000,]
  colnames(b) = bNames
  simDs.id = cbind(jointModelData$trainingData$trainingDs.id, b)
  
  timeDepCoxModels[[i]] = fitTimeDepCoxModelTrainingDs(simDs.id)
}
