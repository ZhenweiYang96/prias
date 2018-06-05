runFixedSchedule = function(testDs.id, biopsyTimes){
  progressionTimes = testDs.id$progression_time
  
  nb = sapply(progressionTime, function(x){which(biopsyTimes>=x)[1]})
  detectionTimes = biopsyTimes[nb]
  offset = detectionTimes - progressionTimes
  
  return(cbind(nb, offset))
}

runPRIASSchedule = function(testDs.id, testDs){
  fixedSchedule = c(1, 4, 7, 10, 15, 20, 25, 30, 35, 40, 45, 50)
  
  psa = 2^testDs$log2psaplus1 - 1
  psa[psa<=0] = 0.01
  testDs$log2psa = log(psa, base = 2)
  
  switchToAnnual = function(ds){
    psaDt = 1/(lm(log2psa~visitTimeYears, data = ds)$coefficients[2])
    return(psaDt>=0 & psaDt<=10)
  }
  
  getNbOffset = function(pid){
    patientDs = testDs[testDs$P_ID == pid,]
    progressionTime = testDs.id$progression_time[testDs.id$P_ID == pid]
    
    nb = 0
    offset = NA
    
    #4 is the minimum number of measurements before which PSA-DT can't be used
    #The biopsy at year 1 will be able to detect the prostate cancer
    if(nrow(patientDs)<4){
      nb = 1
      offset = 1 - progressionTime
    }else{
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
        
        if(switchToAnnual(patientDs[1:j, ])){
          if((curVisitTime - lastBiopsyTime) >= 1){
            proposedBiopsyTime = curVisitTime
          }else{
            proposedBiopsyTime = lastBiopsyTime + 1
          }
        }else{
          proposedBiopsyTime = fixedSchedule[which(fixedSchedule >= curVisitTime)[1]]
          
          if(proposedBiopsyTime - lastBiopsyTime < 1){
            proposedBiopsyTime = lastBiopsyTime + 1
          }
        }
      }
      
      #If we iterated through the entire vector of visit times
      if(proposedBiopsyTime < Inf){
        nb = nb + 1
        offset = proposedBiopsyTime - progressionTime
      }
    }
    
    return(c(nb, offset))
  }
  
  return(t(sapply(testDs.id$P_ID, getNbOffset)))
}


runDynRiskGRSchedule = function(jointModelData, riskMethodName){

  getRiskThreshold = function(riskMethodName, lastBiopsyTime, curVisitTime){
    #if F1 then use last biopsy time and curvisit time
    #else return the fixed threshold value correponding to the name
  }
    
}