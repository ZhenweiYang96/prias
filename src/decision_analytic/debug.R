runHybridSchedule_mod = function(jointModelData){
  
  max_surv_time = max(jointModelData$trainingData$trainingDs$progression_time)
  
  ## Expected failure time:
  #Since we dont have the entire PPD of the patient, we only know the "at least" estimate
  # of the mean failure time.
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
  
  pMedianTimeOptimal = function(newdata, lastBiopsyTime){
    lastVisitTime = max(newdata$visitTimeYears)
    
    round1SurvTimes = c(seq(lastBiopsyTime + 0.0001, max_surv_time, 1), max_surv_time)
    
    round1 = survfitJM(jointModelData$mvJoint_dre_psa_simDs, newdata, idVar = "P_ID", last.time=lastBiopsyTime,
                       survTimes = round1SurvTimes)
    
    times = round1$summaries[[1]][, "times"]
    probs = round1$summaries[[1]][, "Mean"]
    
    #adjusted median prob for max_surv_time year follow up
    #min(probs) corresponds to prob at year max_surv_time
    survProb = 1 - ((1-min(probs))*0.5)  
    
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
  
  #Default is survTime = 0.5, median
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
  
  getSurvThreshold = function(lastBiopsyTime, Dt){
    availableDt = as.numeric(names(thresholdsList))
    Dt = availableDt[which.min(abs(Dt-availableDt))]
    
    cutpoints = thresholdsList[[as.character(Dt)]]$cut_points
    
    lastBiopsyTime = min(max(lastBiopsyTime, cutpoints[1,1]), cutpoints[nrow(cutpoints),1])
    
    #1st column in cutpoints is the column of times
    upperTimeIndex = which(lastBiopsyTime <= cutpoints[,1])[1]
    lowerTimeIndex = tail(which(lastBiopsyTime >= cutpoints[,1]),1)
    
    probDiffPerUnitTimeDiff = (cutpoints[upperTimeIndex, 2] - cutpoints[lowerTimeIndex, 2]) / (cutpoints[upperTimeIndex, 1] - cutpoints[lowerTimeIndex,1])
    
    diffProb = probDiffPerUnitTimeDiff * (lastBiopsyTime - cutpoints[lowerTimeIndex,1])
    
    if(is.nan(diffProb)){
      diffProb = 0
    }
    
    #2nd column in cutpoints is the column for F1 score thresholds
    return(cutpoints[lowerTimeIndex, 2] + diffProb)
  }
  
  testDs.id = jointModelData$testData$testDs.id
  testDs = jointModelData$testData$testDs[jointModelData$testData$testDs$visitTimeYears <= max_surv_time,]
  
  cl = makeCluster(max_cores)
  registerDoParallel(cl)
  
  #testDs.id$P_ID
  results = foreach(pid = testDs.id$P_ID, .combine = "rbind", .packages = c("splines", "JMbayes"),
                    .export=c("MIN_FIXED_ROWS", "MAX_FAIL_TIME", "thresholdsList")) %dopar%{
                      progressionTime = testDs.id$progression_time[testDs.id$P_ID == pid]
                      
                      # print("#######################")
                      print(paste(pid, ":", progressionTime))
                      # print("#######################")
                      
                      lastBiopsyTime = 0
                      nb = 0
                      offset = lastBiopsyTime - progressionTime
                      
                      testDs_j = testDs[testDs$P_ID == pid,]
                      
                      #the first row is time 0
                      for(j in MIN_FIXED_ROWS:nrow(testDs_j)){
                        curVisitTime = testDs_j$visitTimeYears[j]
                        
                        if((curVisitTime - lastBiopsyTime) >= 1){
                          
                          if(curVisitTime<=3.5){
                            sfit2 = survfitJM(object=jointModelData$mvJoint_dre_psa_simDs, 
                                              newdata = testDs_j[1:j,], idVar = "P_ID",
                                              survTimes = curVisitTime, last.time = lastBiopsyTime)
                            
                            meanSurv2 = round(sfit2$summaries[[1]][, "Mean"], 2)
                            proposedBiopsyTime = ifelse(meanSurv2>0.85, Inf, curVisitTime)
                          }else{
                            if(lastBiopsyTime < 9){
                              if(lastBiopsyTime < 3.5){
                                survFit_j = survfitJM(object=jointModelData$mvJoint_dre_psa_simDs, 
                                                      newdata = testDs_j[1:j,], idVar = "P_ID",
                                                      survTimes = 3.5, last.time = lastBiopsyTime)
                                
                                meanSurvProb = round(survFit_j$summaries[[1]][, "Mean"], 2)
                                
                                Dt = 3.5 - lastBiopsyTime
                                if(meanSurvProb <= getSurvThreshold(lastBiopsyTime, Dt)){
                                  sfit2 = survfitJM(object=jointModelData$mvJoint_dre_psa_simDs, 
                                                    newdata = testDs_j[1:j,], idVar = "P_ID",
                                                    survTimes = curVisitTime, last.time = lastBiopsyTime)
                                  
                                  meanSurv2 = round(sfit2$summaries[[1]][, "Mean"], 2)
                                  proposedBiopsyTime = ifelse(meanSurv2>0.85, Inf, curVisitTime)
                                }else{
                                  survFit_j = survfitJM(object=jointModelData$mvJoint_dre_psa_simDs, 
                                                        newdata = testDs_j[1:j,], idVar = "P_ID",
                                                        survTimes = 9, last.time = lastBiopsyTime)
                                  
                                  meanSurvProb = round(survFit_j$summaries[[1]][, "Mean"], 2)
                                  
                                  Dt = 9 - lastBiopsyTime
                                  if(meanSurvProb <= getSurvThreshold(lastBiopsyTime, Dt)){
                                    proposedBiopsyTime = pMedianTimeOptimal(testDs_j[1:j,], lastBiopsyTime)
                                  }else{
                                    proposedBiopsyTime = expectedCondFailureTime(newdata = testDs_j[1:j,], lastBiopsyTime = lastBiopsyTime)
                                  }
                                }
                              }else{
                                survFit_j = survfitJM(object=jointModelData$mvJoint_dre_psa_simDs, 
                                                      newdata = testDs_j[1:j,], idVar = "P_ID",
                                                      survTimes = , last.time = lastBiopsyTime)
                                
                                meanSurvProb = round(survFit_j$summaries[[1]][, "Mean"], 2)
                                
                                Dt = 9 - lastBiopsyTime
                                if(meanSurvProb <= getSurvThreshold(lastBiopsyTime, Dt)){
                                  proposedBiopsyTime = pMedianTimeOptimal(testDs_j[1:j,], lastBiopsyTime)
                                }else{
                                  proposedBiopsyTime = expectedCondFailureTime(newdata = testDs_j[1:j,], lastBiopsyTime = lastBiopsyTime)
                                }
                              }
                            }else{
                              proposedBiopsyTime = expectedCondFailureTime(newdata = testDs_j[1:j,], lastBiopsyTime = lastBiopsyTime)
                            }
                          }
                          
                          #nearest 0.25 years is 3 months
                          if(curVisitTime >= proposedBiopsyTime){
                            lastBiopsyTime = curVisitTime
                            nb = nb + 1
                            offset = curVisitTime - progressionTime
                            
                            if(!is.na(offset) & offset > 0){
                              break
                            }
                          }
                        }
                      }
                      
                      return(c("nb"=nb, "offset"=offset, "lastBiopsyTime"=lastBiopsyTime))
                    }
  stopCluster(cl)
  
  rownames(results) = NULL
  
  return(results)
}
