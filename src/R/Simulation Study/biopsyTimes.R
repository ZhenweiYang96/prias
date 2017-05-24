computeNbAndOffset = function(dsId, patientRowNum, methodName="expectedFailureTime", 
                              minVisits, lastPossibleVisit, biopsyEveryKYears=NA){
  
  trueProgressionTime = simulatedDsList[[dsId]]$testDs.id$progression_time[patientRowNum]
  
  patientId = simulatedDsList[[dsId]]$testDs.id$P_ID[patientRowNum]
  patientDs_i = simulatedDsList[[dsId]]$testDs[simulatedDsList[[dsId]]$testDs$P_ID == patientId,]
  
  visitTimeYears = patientDs_i$visitTimeYears
  dynamicCutOffTimes = simulatedDsList[[dsId]]$dynamicCutOffTimes
  
  nb = 0
  biopsyTimeOffset = NA
  
  lastBiopsyTime = 0
  proposedBiopsyTime = Inf
  enforcedBiopsyTime = Inf
  
  skipProposals = FALSE
  
  for(j in minVisits:lastPossibleVisit){
    curVisitTime = visitTimeYears[j]
    if(curVisitTime > enforcedBiopsyTime){
      lastBiopsyTime = enforcedBiopsyTime
      nb = nb + 1
      biopsyTimeOffset = enforcedBiopsyTime - trueProgressionTime
      enforcedBiopsyTime = Inf
      skipProposals = FALSE
    }
    
    if(!is.na(biopsyTimeOffset) & biopsyTimeOffset > 0){
      break
    }
    
    #If the K years have passed since last biopsy
    if(!is.na(biopsyEveryKYears)){
      if((curVisitTime - lastBiopsyTime) >= biopsyEveryKYears){
        enforcedBiopsyTime = curVisitTime
      }
    }
    
    #Now deal with the proposed biopsy time
    if(skipProposals == FALSE){
      nearest_time_index = which(abs(dynamicCutOffTimes-curVisitTime)==min(abs(dynamicCutOffTimes-curVisitTime)))[1]
      
      if(methodName=="expectedFailureTime"){
        proposedBiopsyTime = expectedCondFailureTime(dsId, patientDs_i[1:j,], lastBiopsyTime)
      }else if(methodName=="survTime85"){
        proposedBiopsyTime = pDynSurvTimeOptimal(dsId, patientDs_i[1:j,], 0.85, lastBiopsyTime)
      }else if(methodName=="survTimeYouden"){
        proposedBiopsyTime = pDynSurvTimeOptimal(dsId, patientDs_i[1:j,], simulatedDsList[[dsId]]$cutoffValues[[nearest_time_index]]["youden"], lastBiopsyTime)
      }else if(methodName=="survTimeAccuracy"){
        proposedBiopsyTime = pDynSurvTimeOptimal(dsId, patientDs_i[1:j,], simulatedDsList[[dsId]]$cutoffValues[[nearest_time_index]]["accuracy"], lastBiopsyTime)
      }else if(methodName=="survTimeF1Score"){
        proposedBiopsyTime = pDynSurvTimeOptimal(dsId, patientDs_i[1:j,], simulatedDsList[[dsId]]$cutoffValues[[nearest_time_index]]["f1score"], lastBiopsyTime)
      }else if(methodName=="survTimeMaxTPR"){
        proposedBiopsyTime = pDynSurvTimeOptimal(dsId, patientDs_i[1:j,], simulatedDsList[[dsId]]$cutoffValues[[nearest_time_index]]["maxTPR"], lastBiopsyTime)
      }
      
      if(!is.na(proposedBiopsyTime)){
        if(proposedBiopsyTime < enforcedBiopsyTime){
          if((proposedBiopsyTime - lastBiopsyTime) > 1){
            enforcedBiopsyTime = proposedBiopsyTime
          }else{
            enforcedBiopsyTime = lastBiopsyTime + 1
            skipProposals = TRUE
          }
        }
      }
    }
    
  }
  
  if(enforcedBiopsyTime < Inf){
    nb = nb + 1
    biopsyTimeOffset = enforcedBiopsyTime - trueProgressionTime
  }
  
  return(c(nb=nb, offset=biopsyTimeOffset))
}

pDynSurvTimeOptimal = function(dsId, patientDs_j, cutoff, lastBiopsyTime){
  lastVisitTime = max(patientDs_j$visitTimeYears)
  
  round1Granularity = 1
  if(15-lastVisitTime < round1Granularity){
    round1Granularity = 15-lastVisitTime
  }
  
  round1 = survfitJM(simulatedDsList[[dsId]]$models$simJointModel_replaced, 
                     newdata = patientDs_j, idVar = "P_ID", last.time=lastBiopsyTime,
                     survTimes = seq(lastVisitTime, 15, round1Granularity))
  
  times = round1$summaries[[1]][, "times"]
  probs = round1$summaries[[1]][, "Median"]
  
  nearest_time = times[which(abs(probs-cutoff)==min(abs(probs-cutoff)))[1]]
  upper = nearest_time + 1
  lower = nearest_time - 1
  if(upper > 15){
    upper = 15
  }
  
  if(lower < 0){
    lower = 0
  }
  newTimes = seq(lower, upper, 0.05)
  
  round2 = survfitJM(simulatedDsList[[dsId]]$models$simJointModel_replaced, 
                     newdata = patientDs_j, idVar = "P_ID", last.time=lastBiopsyTime,
                     survTimes = newTimes)
  
  times = round2$summaries[[1]][, "times"]
  probs = round2$summaries[[1]][, "Median"]
  
  times[which(abs(probs-cutoff)==min(abs(probs-cutoff)))[1]]
}
