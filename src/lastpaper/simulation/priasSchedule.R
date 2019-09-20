runPRIASSchedule = function(testDs.id, testDs){
  fixedSchedule = c(1, 4, 7, 10)
  
  psa = 2^testDs$log2psaplus1 - 1
  psa[psa<=0] = 0.01
  testDs$log2psa = log(psa, base = 2)
  
  switchToAnnual = function(ds){
    psaDt = 1/(lm(log2psa~year_visit, data = ds)$coefficients[2])
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
        curVisitTime = patientDs$year_visit[j]
        
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
