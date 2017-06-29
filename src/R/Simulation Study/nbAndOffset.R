computeNbAndOffset_PRIAS = function(dsId, patientRowNum){
  progressionTime = simulatedDsList[[dsId]]$testDs.id$progression_time[patientRowNum]
  patientId = simulatedDsList[[dsId]]$testDs.id$P_ID[patientRowNum]
  patientDs = simulatedDsList[[dsId]]$testDs[simulatedDsList[[dsId]]$testDs$P_ID == patientId,]
  
  nb = 0
  offset = NA
  
  fixedSchedule = c(1, 4, 7, 10, 15, 20, 25, 30, 35, 40, 45, 50)
  fixedScheduleIndex = 1
  
  skipFixedBiopsies = FALSE
  lastBiopsyTime = 0
  proposedBiopsyTime = Inf
  
  #Min number of measurements before which PSA-DT cant be used
  for(j in 4:nrow(patientDs)){
    curVisitTime = patientDs$visitTimeYears[j]
    
    if(curVisitTime > proposedBiopsyTime){
      nb = nb + 1
      offset = proposedBiopsyTime - progressionTime
      lastBiopsyTime = proposedBiopsyTime
      proposedBiopsyTime = Inf
      skipFixedBiopsies = FALSE
    }
    
    if(!is.na(offset) & offset > 0){
      break
    }
    
    betaReg = lm(logpsa1~visitTimeYears, data = patientDs[1:j, ])$coefficients[2]
    
    if((1/betaReg) <= 10){
      if((curVisitTime - lastBiopsyTime) > 1){
        proposedBiopsyTime = curVisitTime
      }else{
        proposedBiopsyTime = lastBiopsyTime + 1
      }
    }else if(skipFixedBiopsies==FALSE){
      proposedBiopsyTime = fixedSchedule[fixedScheduleIndex]
      fixedScheduleIndex = fixedScheduleIndex + 1
      
      while((proposedBiopsyTime - lastBiopsyTime) < 1){
        proposedBiopsyTime = fixedSchedule[fixedScheduleIndex]
        fixedScheduleIndex = fixedScheduleIndex + 1
      }
      
      skipFixedBiopsies = TRUE
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
  
  return(c(nb = detectionTime-1, offset = detectionTime-progressionTime))
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

computeNbAndOffset = function(dsId, patientRowNum, methodName="expectedFailureTime", 
                              minVisits, lastPossibleVisit, biopsyEveryKYears=NA){
  
  trueProgressionTime = simulatedDsList[[dsId]]$testDs.id$progression_time[patientRowNum]
  
  patientId = simulatedDsList[[dsId]]$testDs.id$P_ID[patientRowNum]
  patientDs_i = simulatedDsList[[dsId]]$testDs[simulatedDsList[[dsId]]$testDs$P_ID == patientId,]
  
  visitTimeYears = patientDs_i$visitTimeYears
  dynamicCutOffTimes = simulatedDsList[[dsId]]$dynamicCutOffTimes
  
  nb = 0
  biopsyTimeOffset = 0 - trueProgressionTime
  
  lastBiopsyTime = 0
  proposedBiopsyTime = Inf
  enforcedBiopsyTime = Inf
  
  nextVisitTime = visitTimeYears[1]
  while(!is.na(nextVisitTime) & biopsyTimeOffset < 0){
    curVisitTime = nextVisitTime
    
    if(curVisitTime == enforcedBiopsyTime){
      lastBiopsyTime = enforcedBiopsyTime
      nb = nb + 1
      biopsyTimeOffset = enforcedBiopsyTime - trueProgressionTime
      enforcedBiopsyTime = Inf

      nextVisitTime = curVisitTime
    }else{
      nearest_time_index = which(abs(dynamicCutOffTimes-curVisitTime)==min(abs(dynamicCutOffTimes-curVisitTime)))[1]
      
      patientDsSubset = patientDs_i[patientDs_i$visitTimeYears <= curVisitTime, ]
      
      if(methodName=="expectedFailureTime"){
        proposedBiopsyTime = expectedCondFailureTime(dsId, patientDsSubset, lastBiopsyTime)
      }else if(methodName=="medianFailureTime"){
        proposedBiopsyTime = pDynSurvTimeOptimal(dsId, patientDsSubset, 0.5, lastBiopsyTime)
      }else if(methodName=="survTimeYouden"){
        proposedBiopsyTime = pDynSurvTimeOptimal(dsId, patientDsSubset, simulatedDsList[[dsId]]$cutoffValues[[nearest_time_index]]["youden"], lastBiopsyTime)
      }else if(methodName=="survTimeAccuracy"){
        proposedBiopsyTime = pDynSurvTimeOptimal(dsId, patientDsSubset, simulatedDsList[[dsId]]$cutoffValues[[nearest_time_index]]["accuracy"], lastBiopsyTime)
      }else if(methodName=="survTimeF1Score"){
        proposedBiopsyTime = pDynSurvTimeOptimal(dsId, patientDsSubset, simulatedDsList[[dsId]]$cutoffValues[[nearest_time_index]]["f1score"], lastBiopsyTime)
      }
      
      #print("Iteration Next")
      #print(paste("Cur:",curVisitTime))
      #print(paste("Last:",lastBiopsyTime))
      #print(paste("Prop:",proposedBiopsyTime))
      #print(paste("Enfor:",enforcedBiopsyTime))
      
      if(!is.na(proposedBiopsyTime)){
        if(proposedBiopsyTime < enforcedBiopsyTime){
          
          if(proposedBiopsyTime < curVisitTime){
            proposedBiopsyTime = curVisitTime
          }
          
          if((proposedBiopsyTime - lastBiopsyTime) > 1){
            enforcedBiopsyTime = proposedBiopsyTime
            nextVisitTime = patientDs_i$visitTimeYears[patientDs_i$visitTimeYears > curVisitTime][1]
          }else{
            enforcedBiopsyTime = lastBiopsyTime + 1
            nextVisitTime = enforcedBiopsyTime
          }
        }#else{
        #  enforcedBiopsyTime = proposedBiopsyTime
        #}
      }else{
        nextVisitTime = patientDs_i$visitTimeYears[patientDs_i$visitTimeYears > curVisitTime][1]
      }
    }
  }
  
  if(enforcedBiopsyTime < Inf){
    nb = nb + 1
    biopsyTimeOffset = enforcedBiopsyTime - trueProgressionTime
  }
  
  return(c(nb=nb, offset=biopsyTimeOffset))
}

computeNbAndOffset_Mixed=function(dsId, patientRowNum, minVisits=1, lastPossibleVisit, 
                                  thresholdSd = 2, alternative="youden", upperLimitIntegral = 15){
  
  #Very specific to this function
  dynamicPredProb = function(futureTimes){
    survProbNearest = sapply(futureTimes, function(futureTime){
      dynamic_probs[which(abs(prob_times-futureTime)==min(abs(prob_times-futureTime)))[1]]
    })
    
    return(survProbNearest)
  }
  
  dynamicPredProbTimesTdiff = function(futureTimes, last.time){
    survProbNearest = sapply(futureTimes, function(futureTime){
      dynamic_probs[which(abs(prob_times-futureTime)==min(abs(prob_times-futureTime)))[1]]
    })
    temp = (futureTimes - last.time) * survProbNearest
    return(temp)
  }
  
  calcSurvProbs = function(futureTimes, dsId, patientDs, last.time=NULL){
    dynamic_probs = c(1, survfitJM(simulatedDsList[[dsId]]$models$simJointModel_replaced, patientDs, 
                          idVar="P_ID",last.time=last.time, 
                          survTimes = prob_times)$summaries[[1]][, "Median"])
    return(dynamic_probs)
  }
  
  expectedCondFailureTime = function(last.time){
    last.time + integrate(dynamicPredProb, lower=last.time,
                                              upper=upperLimitIntegral, abs.tol = 0.1)$value
  }
  
  varianceCondFailureTime = function(last.time){
    
    lhs = 2 * integrate(dynamicPredProbTimesTdiff, last.time, 
                        upperLimitIntegral, last.time=last.time, rel.tol = 0.05)$value
    
    rhs = (integrate(dynamicPredProb,lower=last.time, upperLimitIntegral, abs.tol = 0.1)$value)^2
    
    return(lhs-rhs)
  }
  
  #Now doing the mixed approach
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
  
  useExpectedFailureTime = FALSE
  
  for(j in minVisits:lastPossibleVisit){
    print(j)
    
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
    
    #Now deal with the proposed biopsy time
    if(skipProposals == FALSE){
      nearest_time_index = which(abs(dynamicCutOffTimes-curVisitTime)==min(abs(dynamicCutOffTimes-curVisitTime)))[1]
      
      prob_times = seq(lastBiopsyTime, upperLimitIntegral, 0.1)
      dynamic_probs = calcSurvProbs(prob_times, dsId = dsId, patientDs_i[1:j,], lastBiopsyTime)
      
      if(useExpectedFailureTime==F){
        
        varFailTime = varianceCondFailureTime(lastBiopsyTime)
        
        if(varFailTime < thresholdSd^2){
          useExpectedFailureTime = T
        }else{
          probInterest = simulatedDsList[[dsId]]$cutoffValues[[nearest_time_index]][alternative]
          proposedBiopsyTime = prob_times[which(abs(dynamic_probs-probInterest)==min(abs(dynamic_probs-probInterest)))[1]]
        }
      }
      
      if(useExpectedFailureTime==T){
        proposedBiopsyTime = expectedCondFailureTime(lastBiopsyTime)  
      }

      if(!is.na(proposedBiopsyTime)){
        if(proposedBiopsyTime < enforcedBiopsyTime){
          
          if(proposedBiopsyTime < curVisitTime){
            proposedBiopsyTime = curVisitTime
          }
          
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