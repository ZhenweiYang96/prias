computeNbAndOffset_PRIAS = function(dsId, patientRowNum){
  progressionTime = simulatedDsList[[dsId]]$testDs.id$progression_time[patientRowNum]
  patientId = simulatedDsList[[dsId]]$testDs.id$P_ID[patientRowNum]
  patientDs = simulatedDsList[[dsId]]$testDs[simulatedDsList[[dsId]]$testDs$P_ID == patientId,]
  
  nb = 0
  offset = NA
  
  fixedSchedule = c(1, 4, 7, 10, 15, 20, 25, 30, 35, 40, 45, 50)
  
  lastBiopsyTime = 0
  proposedBiopsyTime = Inf
  
  #Min number of measurements before which PSA-DT can't be used
  for(j in 4:nrow(patientDs)){
    curVisitTime = patientDs$visitTimeYears[j]
    
    if(curVisitTime > proposedBiopsyTime){
      nb = nb + 1
      offset = proposedBiopsyTime - progressionTime
      lastBiopsyTime = proposedBiopsyTime
      proposedBiopsyTime = Inf
    }
    
    if(!is.na(offset) & offset > 0){
      break
    }
    
    betaReg = lm(logpsa1~visitTimeYears, data = patientDs[1:j, ])$coefficients[2]
    
    #in the paper we submitted the 1/betaReg >= 0 condition was not present
    if((1/betaReg) <= 10 & (1/betaReg)>=0){
      if((curVisitTime - lastBiopsyTime) > 1){
        proposedBiopsyTime = curVisitTime
      }else{
        proposedBiopsyTime = lastBiopsyTime + 1
      }
    }else{
      fixedScheduleIndex = 1
      proposedBiopsyTime = fixedSchedule[fixedScheduleIndex]
      fixedScheduleIndex = fixedScheduleIndex + 1
      
      while((proposedBiopsyTime - lastBiopsyTime) <= 0){
        proposedBiopsyTime = fixedSchedule[fixedScheduleIndex]
        fixedScheduleIndex = fixedScheduleIndex + 1
      }
      
      if(proposedBiopsyTime - lastBiopsyTime < 1){
        proposedBiopsyTime = lastBiopsyTime + 1
      }
    }
  }
  
  if(proposedBiopsyTime < Inf){
    nb = nb + 1
    offset = proposedBiopsyTime - progressionTime
  }
  
  return(c(nb = nb, offset = offset))
}

computeNbAndOffset_JH = function(dsId, patientRowNum){
  progressionTime = simulatedDsList[[dsId]]$testDs.id$progression_time[patientRowNum]
  detectionTime = ceiling(progressionTime)
  
  return(c(nb = detectionTime, offset = detectionTime-progressionTime))
}

pDynSurvTimeOptimal = function(dsId, patientDs_j, cutoff, lastBiopsyTime){
  lastVisitTime = max(patientDs_j$visitTimeYears)
  
  round1SurvTimes = seq(lastBiopsyTime + 0.0001, lastBiopsyTime + 15, 1)
  
  round1 = survfitJM(simulatedDsList[[dsId]]$models$simJointModel_replaced, 
                     newdata = patientDs_j, idVar = "P_ID", last.time=lastBiopsyTime,
                     survTimes = round1SurvTimes)
  
  times = round1$summaries[[1]][, "times"]
  probs = round1$summaries[[1]][, "Median"]
  
  probs[probs < 0.009] = 0
  
  nearest_time = times[which(abs(probs-cutoff)==min(abs(probs-cutoff)))[1]]
  upper = nearest_time + 1
  lower = nearest_time - 1
 
  if(lower <= lastBiopsyTime){
    lower = lastBiopsyTime + 0.0001  
  }
  
  round2SurvTimes = seq(lower, upper, 0.05)
  
  round2 = survfitJM(simulatedDsList[[dsId]]$models$simJointModel_replaced, 
                     newdata = patientDs_j, idVar = "P_ID", last.time=lastBiopsyTime,
                     survTimes = round2SurvTimes)
  
  times = round2$summaries[[1]][, "times"]
  probs = round2$summaries[[1]][, "Median"]
  
  times[which(abs(probs-cutoff)==min(abs(probs-cutoff)))[1]]
}

computeNbAndOffset = function(dsId, patientRowNum, methodName="expectedFailureTime", 
                              minVisits=1, lastPossibleVisit = timesPerSubject-1, biopsyEveryKYears=NA){
  
  trueProgressionTime = simulatedDsList[[dsId]]$testDs.id$progression_time[patientRowNum]
  
  patientId = simulatedDsList[[dsId]]$testDs.id$P_ID[patientRowNum]
  patientDs_i = simulatedDsList[[dsId]]$testDs[simulatedDsList[[dsId]]$testDs$P_ID == patientId,]
  
  visitTimeYears = patientDs_i$visitTimeYears
  dynamicCutOffTimes = simulatedDsList[[dsId]]$dynamicCutOffTimes
  
  nb = 0
  biopsyTimeOffset = 0 - trueProgressionTime
  
  lastBiopsyTime = 0
  proposedBiopsyTime = Inf
  scheduledBiopsyTime = Inf
  
  nextVisitTime = visitTimeYears[1]
  while(!is.na(nextVisitTime) & biopsyTimeOffset < 0){
    curVisitTime = nextVisitTime
    
    #print(paste(curVisitTime, "----", scheduledBiopsyTime))
    
    if(curVisitTime >= scheduledBiopsyTime){
      lastBiopsyTime = scheduledBiopsyTime
      nextVisitTime = scheduledBiopsyTime
      
      nb = nb + 1
      biopsyTimeOffset = scheduledBiopsyTime - trueProgressionTime
      scheduledBiopsyTime = Inf
    }else{
      patientDsSubset = patientDs_i[patientDs_i$visitTimeYears <= curVisitTime, ]
      
      if(methodName=="expectedFailureTime"){
        proposedBiopsyTime = expectedCondFailureTime(dsId, patientDsSubset, lastBiopsyTime)
      }else if(methodName=="medianFailureTime"){
        proposedBiopsyTime = pDynSurvTimeOptimal(dsId, patientDsSubset, 0.5, lastBiopsyTime)
      }else if(methodName %in% c("youden", "accuracy", "f1score")){
        #nearest_time_index = which(abs(dynamicCutOffTimes-curVisitTime)==min(abs(dynamicCutOffTimes-curVisitTime)))[1]
        nearest_time_index = which(abs(dynamicCutOffTimes-lastBiopsyTime)==min(abs(dynamicCutOffTimes-lastBiopsyTime)))[1]
        
        cutoff =  simulatedDsList[[dsId]]$cutoffValues[[nearest_time_index]][methodName]
        if(!is.na(cutoff)){
          proposedBiopsyTime = pDynSurvTimeOptimal(dsId, patientDsSubset, cutoff, lastBiopsyTime)
        }else{
          proposedBiopsyTime = NA
        }
      }else if(methodName %in% c("Cut95")){
        proposedBiopsyTime = pDynSurvTimeOptimal(dsId, patientDsSubset, 0.95, lastBiopsyTime)
      }
      
      #print("Iteration Next")
      #print(paste("Cur:",curVisitTime))
      #print(paste("Last:",lastBiopsyTime))
      #print(paste("Prop:",proposedBiopsyTime))
      #print(paste("Enfor:",scheduledBiopsyTime))
      
      if(!is.na(proposedBiopsyTime)){
        if(proposedBiopsyTime < scheduledBiopsyTime){
          
          if(proposedBiopsyTime < curVisitTime){
            proposedBiopsyTime = curVisitTime
          }
          
          if((proposedBiopsyTime - lastBiopsyTime) > 1){
            scheduledBiopsyTime = proposedBiopsyTime
            nextVisitTime = patientDs_i$visitTimeYears[patientDs_i$visitTimeYears > curVisitTime][1]
          }else{
            scheduledBiopsyTime = lastBiopsyTime + 1
            nextVisitTime = scheduledBiopsyTime
          }
        }else{
          #scheduledBiopsyTime = proposedBiopsyTime
          nextVisitTime = patientDs_i$visitTimeYears[patientDs_i$visitTimeYears > curVisitTime][1]
        }
      }else{
        nextVisitTime = patientDs_i$visitTimeYears[patientDs_i$visitTimeYears > curVisitTime][1]
      }
    }
  }
  
  if(scheduledBiopsyTime < Inf){
    nb = nb + 1
    biopsyTimeOffset = scheduledBiopsyTime - trueProgressionTime
  }
  
  return(c(nb=nb, offset=biopsyTimeOffset))
}

computeNbAndOffset_Mixed=function(dsId, patientRowNum, minVisits=1, lastPossibleVisit, alternative="youden"){
  
  trueProgressionTime = simulatedDsList[[dsId]]$testDs.id$progression_time[patientRowNum]
  
  patientId = simulatedDsList[[dsId]]$testDs.id$P_ID[patientRowNum]
  patientDs_i = simulatedDsList[[dsId]]$testDs[simulatedDsList[[dsId]]$testDs$P_ID == patientId,]
  
  visitTimeYears = patientDs_i$visitTimeYears
  dynamicCutOffTimes = simulatedDsList[[dsId]]$dynamicCutOffTimes
  
  nb = 0
  biopsyTimeOffset = 0 - trueProgressionTime
  
  lastBiopsyTime = 0
  proposedBiopsyTime = Inf
  scheduledBiopsyTime = Inf
  
  nextVisitTime = visitTimeYears[1]
  while(!is.na(nextVisitTime) & biopsyTimeOffset < 0){
    curVisitTime = nextVisitTime
    
    #print(paste(curVisitTime, "----", scheduledBiopsyTime))
    
    if(curVisitTime >= scheduledBiopsyTime){
      lastBiopsyTime = scheduledBiopsyTime
      nextVisitTime = scheduledBiopsyTime
      
      nb = nb + 1
      biopsyTimeOffset = scheduledBiopsyTime - trueProgressionTime
      scheduledBiopsyTime = Inf
    }else{
      patientDsSubset = patientDs_i[patientDs_i$visitTimeYears <= curVisitTime, ]
      
      pt025quantile = pDynSurvTimeOptimal(dsId, patientDsSubset, 0.975, lastBiopsyTime)
      medianFailureTime = pDynSurvTimeOptimal(dsId, patientDsSubset, 0.5, lastBiopsyTime)
      #medianFailureTime = expectedCondFailureTime(dsId, patientDsSubset, lastBiopsyTime)
      
      #nearest_time_index = which(abs(dynamicCutOffTimes-curVisitTime)==min(abs(dynamicCutOffTimes-curVisitTime)))[1]
      nearest_time_index = which(abs(dynamicCutOffTimes-lastBiopsyTime)==min(abs(dynamicCutOffTimes-lastBiopsyTime)))[1]
      
      cutoff =  simulatedDsList[[dsId]]$cutoffValues[[nearest_time_index]][alternative]
      if((medianFailureTime - pt025quantile > 3) & !is.na(cutoff)){
        proposedBiopsyTime = pDynSurvTimeOptimal(dsId, patientDsSubset, cutoff, lastBiopsyTime)
      }else{
        proposedBiopsyTime = medianFailureTime
      }
      
      if(!is.na(proposedBiopsyTime)){
        if(proposedBiopsyTime < scheduledBiopsyTime){
          
          if(proposedBiopsyTime < curVisitTime){
            proposedBiopsyTime = curVisitTime
          }
          
          if((proposedBiopsyTime - lastBiopsyTime) > 1){
            scheduledBiopsyTime = proposedBiopsyTime
            nextVisitTime = patientDs_i$visitTimeYears[patientDs_i$visitTimeYears > curVisitTime][1]
          }else{
            scheduledBiopsyTime = lastBiopsyTime + 1
            nextVisitTime = scheduledBiopsyTime
          }
        }else{
          #scheduledBiopsyTime = proposedBiopsyTime
          nextVisitTime = patientDs_i$visitTimeYears[patientDs_i$visitTimeYears > curVisitTime][1]
        }
      }else{
        nextVisitTime = patientDs_i$visitTimeYears[patientDs_i$visitTimeYears > curVisitTime][1]
      }
    }
  }
  
  if(scheduledBiopsyTime < Inf){
    nb = nb + 1
    biopsyTimeOffset = scheduledBiopsyTime - trueProgressionTime
  }
  
  return(c(nb=nb, offset=biopsyTimeOffset))
}