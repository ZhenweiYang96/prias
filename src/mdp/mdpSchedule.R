#riskMethodName is either F1score
#These values should match the column names that Dimitris returns in the cutpoints object
runMDPSchedule_1biopsyonly = function(jointModelData){
  
  #Remember survival estimates are valid only upto the last progression time in the training dataset
  #doesnt matter if it is censoring time or progression time. This because baseline hazard is only valid upto that point.
  #Check when this gives a vector of all zeroes to confirm
  #splineDesign(jointModelData$mvJoint_dre_psa_simDs$control$knots, <time>, ord = 4, outer.ok = T)
  
  max_surv_time = max(jointModelData$trainingData$trainingDs$progression_time) - 10e-2
  
  getSurvThreshold = function(lastBiopsyTime, curVisitTime){
    #Dt = curVisitTime - lastBiopsyTime
    Dt = 0.5
    
    availableDt = as.numeric(names(thresholdsList))
    Dt = availableDt[which.min(abs(Dt-availableDt))]
    
    cutpoints = thresholdsList[[as.character(Dt)]]$cut_points
    
    lastBiopsyTime = min(max(lastBiopsyTime, cutpoints[1,1]), cutpoints[nrow(cutpoints),1])
    
    #1st column in cutpoints is the column of times
    upperTimeIndex = which(lastBiopsyTime <= cutpoints[,1])[1]
    lowerTimeIndex = tail(which(lastBiopsyTime >= cutpoints[,1]),1)
    
    probDiffPerUnitTimeDiff = (cutpoints[upperTimeIndex, riskMethodName] - cutpoints[lowerTimeIndex, riskMethodName]) / (cutpoints[upperTimeIndex, 1] - cutpoints[lowerTimeIndex,1])
    
    diffProb = probDiffPerUnitTimeDiff * (lastBiopsyTime - cutpoints[lowerTimeIndex,1])
    
    if(is.nan(diffProb)){
      diffProb = 0
    }
    
    #2nd column in cutpoints is the column for F1 score thresholds
    return(cutpoints[lowerTimeIndex, riskMethodName] + diffProb)
  }
  
  expectedCondFailureTime = function(newdata, lastBiopsyTime){
    
    dynamicPredProb = function(futureTimes, newdata, lastBiopsyTime){
      return(survfitJM(jointModelData$mvJoint_dre_psa_simDs, newdata, last.time = lastBiopsyTime,
                       idVar="P_ID", survTimes = futureTimes)$summaries[[1]][, "Mean"])
    }
    
    upper = max_surv_time
    #Change rel.tol = 0.05 to make it faster
    repeat{
      integration_result = try(integrate(dynamicPredProb, lower=lastBiopsyTime, 
                                         upper=upper, newdata=newdata, lastBiopsyTime = lastBiopsyTime, rel.tol = 0.05)$value, silent = T)
      
      if(!inherits(integration_result, "try-error")){
        break
      }else{
        upper = upper - 0.1
        if(upper < max(9, lastBiopsyTime)){
          #Just a large number
          integration_result = 20
          break
        }
      }
    }
    
    lastBiopsyTime + integration_result
  }
  
  pDynSurvTimeOptimal = function(newdata, lastBiopsyTime, survProb = NA){
    lastVisitTime = max(newdata$visitTimeYears)
    
    round1SurvTimes = c(seq(lastBiopsyTime + 0.0001, max_surv_time, 1), max_surv_time)
    
    round1 = survfitJM(jointModelData$mvJoint_dre_psa_simDs, newdata, idVar = "P_ID", last.time=lastBiopsyTime,
                       survTimes = round1SurvTimes)
    
    times = round1$summaries[[1]][, "times"]
    probs = round1$summaries[[1]][, "Mean"]
    
    nearest_time = times[which.min(abs(probs-survProb))[1]]
    upper = min(max_surv_time, nearest_time + 1)
    lower = max(lastBiopsyTime + 1e-4, nearest_time - 1)
    
    round2SurvTimes = seq(lower, upper, 0.05)
    
    round2 = survfitJM(jointModelData$mvJoint_dre_psa_simDs, newdata, idVar = "P_ID", last.time=lastBiopsyTime,
                       survTimes = round2SurvTimes)
    
    times = round2$summaries[[1]][, "times"]
    probs = round2$summaries[[1]][, "Mean"]
    
    times[which.min(abs(probs-survProb))[1]]
  }
  
  testDs.id = jointModelData$testData$testDs.id
  testDs = jointModelData$testData$testDs[jointModelData$testData$testDs$visitTimeYears <= max_surv_time,]
  
  getWeightedOffset = function(visitTime, lastBiomarkerTime, lastBiopsyTime,
                               medianFailureTime){
    survFit_j = survfitJM(object=jointModelData$mvJoint_dre_psa_simDs, 
                          newdata = testDs_j[testDs_j$visitTimeYears<=lastBiomarkerTime,], idVar = "P_ID",
                          survTimes = visitTime, last.time = lastBiopsyTime)
    
    meanSurvProb = round(survFit_j$summaries[[1]][, "Mean"], 2)
    #survThreshold = getSurvThreshold(lastBiopsyTime, curVisitTime)
    meanFailProb = 1 - meanSurvProb
    
    expectedOffsetIfFail = (visitTime - lastBiopsyTime)/2
    #expectedOffsetIfSurv = curVisitTime - expectedCondFailureTime(newdata = testDs_j[1:j,], lastBiopsyTime = curVisitTime)
    
    ifsurvOffset = meanSurvProb * (max_surv_time - medianFailureTime)
    
    if(ifsurvOffset == 0){
      return(Inf)
    }else{
      weightedOffset = ifsurvOffset + meanFailProb * expectedOffsetIfFail
    }
    return(weightedOffset)
  }
  
  
  cl = makeCluster(max_cores)
  registerDoParallel(cl)
  
  results = foreach(pid = testDs.id$P_ID[1:100], .combine = "rbind", .packages = c("splines", "JMbayes"),
                    .export=c("MIN_FIXED_ROWS")) %dopar%{
                      progressionTime = testDs.id$progression_time[testDs.id$P_ID == pid]
                      
                      lastBiopsyTime = 0
                      nb = 0
                      offset = lastBiopsyTime - progressionTime
                      
                      testDs_j = testDs[testDs$P_ID == pid,]
                      
                      #the first row is time 0
                      for(j in MIN_FIXED_ROWS:nrow(testDs_j)){
                        curVisitTime = testDs_j$visitTimeYears[j]
                        
                        if((curVisitTime - lastBiopsyTime) >= 1){
                          medianFailureTime = pDynSurvTimeOptimal(newdata = testDs_j[1:j,], 
                                                                  lastBiopsyTime = lastBiopsyTime,
                                                                  survProb = 0.15)
                          
                          weightedOffsets = sapply(j:(min(j+6, 24)), function(k){
                            weightedOffset = getWeightedOffset(visitTime = testDs_j$visitTimeYears[k], 
                                            lastBiomarkerTime = curVisitTime, 
                                            lastBiopsyTime =  lastBiopsyTime, 
                                            medianFailureTime = medianFailureTime)
                          })
                          
                          #print(weightedOffsets)
                          #print(which.min(weightedOffsets))
                          
                          minIndex = which.min(weightedOffsets)
                          minoffset = weightedOffsets[minIndex]
                                                    
                          if(minIndex == 1 & minoffset!=Inf){
                            lastBiopsyTime = curVisitTime
                            nb = nb + 1
                            offset = curVisitTime - progressionTime
                            
                            if(!is.na(offset) & offset > 0){
                              break
                            }
                          }
                        }
                      }
                      
                      if(offset<0){
                        nb = nb+1
                        offset = 10 - progressionTime
                      }
                      
                      return(c("nb"=nb, "offset"=offset))
                    }
  stopCluster(cl)
  
  rownames(results) = NULL
  
  return(results)
}

#riskMethodName is either F1score
#These values should match the column names that Dimitris returns in the cutpoints object
runMDPSchedule_dynsurv = function(jointModelData){
  
  #Remember survival estimates are valid only upto the last progression time in the training dataset
  #doesnt matter if it is censoring time or progression time. This because baseline hazard is only valid upto that point.
  #Check when this gives a vector of all zeroes to confirm
  #splineDesign(jointModelData$mvJoint_dre_psa_simDs$control$knots, <time>, ord = 4, outer.ok = T)
  
  max_surv_time = max(jointModelData$trainingData$trainingDs$progression_time) - 10e-2
  
  getSurvThreshold = function(lastBiopsyTime, curVisitTime){
    #Dt = curVisitTime - lastBiopsyTime
    Dt = 0.5
    
    availableDt = as.numeric(names(thresholdsList))
    Dt = availableDt[which.min(abs(Dt-availableDt))]
    
    cutpoints = thresholdsList[[as.character(Dt)]]$cut_points
    
    lastBiopsyTime = min(max(lastBiopsyTime, cutpoints[1,1]), cutpoints[nrow(cutpoints),1])
    
    #1st column in cutpoints is the column of times
    upperTimeIndex = which(lastBiopsyTime <= cutpoints[,1])[1]
    lowerTimeIndex = tail(which(lastBiopsyTime >= cutpoints[,1]),1)
    
    probDiffPerUnitTimeDiff = (cutpoints[upperTimeIndex, riskMethodName] - cutpoints[lowerTimeIndex, riskMethodName]) / (cutpoints[upperTimeIndex, 1] - cutpoints[lowerTimeIndex,1])
    
    diffProb = probDiffPerUnitTimeDiff * (lastBiopsyTime - cutpoints[lowerTimeIndex,1])
    
    if(is.nan(diffProb)){
      diffProb = 0
    }
    
    #2nd column in cutpoints is the column for F1 score thresholds
    return(cutpoints[lowerTimeIndex, riskMethodName] + diffProb)
  }
  
  expectedCondFailureTime = function(newdata, lastBiopsyTime){
    
    dynamicPredProb = function(futureTimes, newdata, lastBiopsyTime){
      return(survfitJM(jointModelData$mvJoint_dre_psa_simDs, newdata, last.time = lastBiopsyTime,
                       idVar="P_ID", survTimes = futureTimes)$summaries[[1]][, "Mean"])
    }
    
    upper = max_surv_time
    #Change rel.tol = 0.05 to make it faster
    repeat{
      integration_result = try(integrate(dynamicPredProb, lower=lastBiopsyTime, 
                                         upper=upper, newdata=newdata, lastBiopsyTime = lastBiopsyTime, rel.tol = 0.05)$value, silent = T)
      
      if(!inherits(integration_result, "try-error")){
        break
      }else{
        upper = upper - 0.1
        if(upper < max(9, lastBiopsyTime)){
          #Just a large number
          integration_result = 20
          break
        }
      }
    }
    
    lastBiopsyTime + integration_result
  }
  
  pDynSurvTimeOptimal = function(newdata, lastBiopsyTime, survProb = NA){
    lastVisitTime = max(newdata$visitTimeYears)
    
    round1SurvTimes = c(seq(lastBiopsyTime + 0.0001, max_surv_time, 1), max_surv_time)
    
    round1 = survfitJM(jointModelData$mvJoint_dre_psa_simDs, newdata, idVar = "P_ID", last.time=lastBiopsyTime,
                       survTimes = round1SurvTimes)
    
    times = round1$summaries[[1]][, "times"]
    probs = round1$summaries[[1]][, "Mean"]
    
    nearest_time = times[which.min(abs(probs-survProb))[1]]
    upper = min(max_surv_time, nearest_time + 1)
    lower = max(lastBiopsyTime + 1e-4, nearest_time - 1)
    
    round2SurvTimes = seq(lower, upper, 0.05)
    
    round2 = survfitJM(jointModelData$mvJoint_dre_psa_simDs, newdata, idVar = "P_ID", last.time=lastBiopsyTime,
                       survTimes = round2SurvTimes)
    
    times = round2$summaries[[1]][, "times"]
    probs = round2$summaries[[1]][, "Mean"]
    
    times[which.min(abs(probs-survProb))[1]]
  }
  
  testDs.id = jointModelData$testData$testDs.id
  testDs = jointModelData$testData$testDs[jointModelData$testData$testDs$visitTimeYears <= max_surv_time,]
  
  getWeightedOffset = function(visitTime, lastBiomarkerTime, lastBiopsyTime,
                               medianFailureTime){
    survFit_j = survfitJM(object=jointModelData$mvJoint_dre_psa_simDs, 
                          newdata = testDs_j[testDs_j$visitTimeYears<=lastBiomarkerTime,], idVar = "P_ID",
                          survTimes = visitTime, last.time = lastBiopsyTime)
    
    meanSurvProb = round(survFit_j$summaries[[1]][, "Mean"], 2)
    #survThreshold = getSurvThreshold(lastBiopsyTime, curVisitTime)
    meanFailProb = 1 - meanSurvProb
    
    expectedOffsetIfFail = (visitTime - lastBiopsyTime)/2
    #expectedOffsetIfSurv = curVisitTime - expectedCondFailureTime(newdata = testDs_j[1:j,], lastBiopsyTime = curVisitTime)
    
    ifsurvOffset = meanSurvProb * (max_surv_time - medianFailureTime)
    
    if(ifsurvOffset == 0){
      return(Inf)
    }else{
      weightedOffset = ifsurvOffset + meanFailProb * expectedOffsetIfFail
    }
    return(weightedOffset)
  }
  
  
  cl = makeCluster(max_cores)
  registerDoParallel(cl)
  
  results = foreach(pid = testDs.id$P_ID[1:100], .combine = "rbind", .packages = c("splines", "JMbayes"),
                    .export=c("MIN_FIXED_ROWS")) %dopar%{
                      progressionTime = testDs.id$progression_time[testDs.id$P_ID == pid]
                      
                      lastBiopsyTime = 0
                      nb = 0
                      offset = lastBiopsyTime - progressionTime
                      
                      testDs_j = testDs[testDs$P_ID == pid,]
                      
                      #the first row is time 0
                      for(j in MIN_FIXED_ROWS:nrow(testDs_j)){
                        curVisitTime = testDs_j$visitTimeYears[j]
                        
                        if((curVisitTime - lastBiopsyTime) >= 1){
                          medianFailureTime = pDynSurvTimeOptimal(newdata = testDs_j[1:j,], 
                                                                  lastBiopsyTime = lastBiopsyTime,
                                                                  survProb = 0.5)
                          
                          weightedOffsets = sapply(j:(min(j+6, 24)), function(k){
                            weightedOffset = getWeightedOffset(visitTime = testDs_j$visitTimeYears[k], 
                                                               lastBiomarkerTime = curVisitTime, 
                                                               lastBiopsyTime =  lastBiopsyTime, 
                                                               medianFailureTime = medianFailureTime)
                          })
                          
                          #print(weightedOffsets)
                          #print(which.min(weightedOffsets))
                          
                          minIndex = which.min(weightedOffsets)
                          minoffset = weightedOffsets[minIndex]
                          
                          if(minIndex == 1 & minoffset!=Inf){
                            lastBiopsyTime = curVisitTime
                            nb = nb + 1
                            offset = curVisitTime - progressionTime
                            
                            if(!is.na(offset) & offset > 0){
                              break
                            }
                          }
                        }
                      }
                      
                      if(offset<0){
                        nb = nb+1
                        offset = 10 - progressionTime
                      }
                      
                      return(c("nb"=nb, "offset"=offset))
                    }
  stopCluster(cl)
  
  rownames(results) = NULL
  
  return(results)
}

temp = runMDPSchedule_1biopsyonly(jointModelData)
