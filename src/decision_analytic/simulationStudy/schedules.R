runFixedSchedule = function(testDs.id, biopsyTimes){
  progressionTimes = pmin(max(biopsyTimes), testDs.id$progression_time)
  
  nb = sapply(progressionTimes, function(x){which(biopsyTimes>=x)[1]})
  detectionTimes = biopsyTimes[nb]
  offset = detectionTimes - progressionTimes
  
  return(cbind(nb, offset))
}

runPRIASSchedule = function(testDs.id, testDs){
  fixedSchedule = c(1, 4, 7, 10)
  
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

runFixedRiskGRSchedule = function(jointModelData, survThreshold){
  
  max_surv_time = max(jointModelData$trainingData$trainingDs$progression_time) - 10e-2
  
  testDs.id = jointModelData$testData$testDs.id
  testDs = jointModelData$testData$testDs[jointModelData$testData$testDs$visitTimeYears <= max_surv_time,]
  
  cl = makeCluster(max_cores)
  registerDoParallel(cl)
  
  results = foreach(pid = testDs.id$P_ID, .combine = "rbind", .packages = c("splines", "JMbayes"),
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
                          survFit_j = survfitJM(object=jointModelData$mvJoint_dre_psa_simDs, 
                                                newdata = testDs_j[1:j,], idVar = "P_ID",
                                                survTimes = curVisitTime, last.time = lastBiopsyTime)
                          
                          meanSurvProb = round(survFit_j$summaries[[1]][, "Mean"], 2)
                          if(meanSurvProb <= survThreshold){
                            lastBiopsyTime = curVisitTime
                            nb = nb + 1
                            offset = curVisitTime - progressionTime
                            
                            if(!is.na(offset) & offset > 0){
                              break
                            }
                          }
                        }
                      }
                      
                      return(c("nb"=nb, "offset"=offset))
                    }
  stopCluster(cl)
  
  rownames(results) = NULL
  
  return(results)
}

#riskMethodName is either F1score, or Youden. 
#These values should match the column names that Dimitris returns in the cutpoints object
runDynRiskGRSchedule = function(jointModelData, riskMethodName="F1score"){
  
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
  
  testDs.id = jointModelData$testData$testDs.id
  testDs = jointModelData$testData$testDs[jointModelData$testData$testDs$visitTimeYears <= max_surv_time,]
  
  cl = makeCluster(max_cores)
  registerDoParallel(cl)
  
  results = foreach(pid = testDs.id$P_ID, .combine = "rbind", .packages = c("splines", "JMbayes"),
                    .export=c("MIN_FIXED_ROWS", "thresholdsList")) %dopar%{
                      progressionTime = testDs.id$progression_time[testDs.id$P_ID == pid]
                      
                      lastBiopsyTime = 0
                      nb = 0
                      offset = lastBiopsyTime - progressionTime
                      
                      testDs_j = testDs[testDs$P_ID == pid,]
                      
                      #the first row is time 0
                      for(j in MIN_FIXED_ROWS:nrow(testDs_j)){
                        curVisitTime = testDs_j$visitTimeYears[j]
                        
                        if((curVisitTime - lastBiopsyTime) >= 1){
                          survFit_j = survfitJM(object=jointModelData$mvJoint_dre_psa_simDs, 
                                                newdata = testDs_j[1:j,], idVar = "P_ID",
                                                survTimes = curVisitTime, last.time = lastBiopsyTime)
                          
                          meanSurvProb = round(survFit_j$summaries[[1]][, "Mean"], 2)
                          survThreshold = getSurvThreshold(lastBiopsyTime, curVisitTime)
                          
                          if(meanSurvProb <= survThreshold){
                            lastBiopsyTime = curVisitTime
                            nb = nb + 1
                            offset = curVisitTime - progressionTime
                            
                            if(!is.na(offset) & offset > 0){
                              break
                            }
                          }
                        }
                      }
                      
                      return(c("nb"=nb, "offset"=offset))
                    }
  stopCluster(cl)
  
  rownames(results) = NULL
  
  return(results)
}

runExpectedFailureTimeSchedule = function(jointModelData){
  
  max_surv_time = max(jointModelData$trainingData$trainingDs$progression_time) - 10e-2
  
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
  
  testDs.id = jointModelData$testData$testDs.id
  testDs = jointModelData$testData$testDs[jointModelData$testData$testDs$visitTimeYears <= max_surv_time,]
  
  cl = makeCluster(max_cores)
  registerDoParallel(cl)
  
  #testDs.id$P_ID
  results = foreach(pid = testDs.id$P_ID, .combine = "rbind", .packages = c("splines", "JMbayes"),
                    .export=c("MIN_FIXED_ROWS", "MAX_FAIL_TIME", "thresholdsList")) %dopar%{
                      progressionTime = testDs.id$progression_time[testDs.id$P_ID == pid]
                      
                      lastBiopsyTime = 0
                      nb = 0
                      offset = lastBiopsyTime - progressionTime
                      
                      testDs_j = testDs[testDs$P_ID == pid,]
                      
                      #the first row is time 0
                      for(j in MIN_FIXED_ROWS:nrow(testDs_j)){
                        curVisitTime = testDs_j$visitTimeYears[j]
                        
                        if((curVisitTime - lastBiopsyTime) >= 1){
                          proposedBiopsyTime = expectedCondFailureTime(newdata = testDs_j[1:j,], lastBiopsyTime = lastBiopsyTime)
                          
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
                      
                      return(c("nb"=nb, "offset"=offset))
                    }
  stopCluster(cl)
  
  rownames(results) = NULL
  
  return(results)
}

runMedianFailureTimeSchedule = function(jointModelData){
  
  max_surv_time = max(jointModelData$trainingData$trainingDs$progression_time) - 10e-2
  
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
  
  testDs.id = jointModelData$testData$testDs.id
  testDs = jointModelData$testData$testDs[jointModelData$testData$testDs$visitTimeYears <= max_surv_time,]
  
  cl = makeCluster(max_cores)
  registerDoParallel(cl)
  
  #testDs.id$P_ID
  results = foreach(pid = testDs.id$P_ID, .combine = "rbind", .packages = c("splines", "JMbayes"),
                    .export=c("MIN_FIXED_ROWS", "MAX_FAIL_TIME", "thresholdsList")) %dopar%{
                      progressionTime = testDs.id$progression_time[testDs.id$P_ID == pid]
                      
                      lastBiopsyTime = 0
                      nb = 0
                      offset = lastBiopsyTime - progressionTime
                      
                      testDs_j = testDs[testDs$P_ID == pid,]
                      
                      #the first row is time 0
                      for(j in MIN_FIXED_ROWS:nrow(testDs_j)){
                        curVisitTime = testDs_j$visitTimeYears[j]
                        
                        if((curVisitTime - lastBiopsyTime) >= 1){
                          proposedBiopsyTime = pMedianTimeOptimal(newdata = testDs_j[1:j,], lastBiopsyTime = lastBiopsyTime)
                          
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
                      
                      return(c("nb"=nb, "offset"=offset))
                    }
  stopCluster(cl)
  
  rownames(results) = NULL
  
  return(results)
}

#####################################################
# Hybrid schedule number 1
# we use 15% rule if the last biopsy is before 3.5 years and if the ROC tells us he can fail in 3.5 years
# we use Median failure time if the last biopsy is before 8 years and the if ROC tells us he can fail in 8 years
# otherwise we use mean failure time
#####################################################
runHybridSchedule = function(jointModelData){
  
  max_surv_time = max(jointModelData$trainingData$trainingDs$progression_time) - 10e-2
  
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
                      
                      lastBiopsyTime = 0
                      nb = 0
                      offset = lastBiopsyTime - progressionTime
                      
                      testDs_j = testDs[testDs$P_ID == pid,]
                      
                      #the first row is time 0
                      for(j in MIN_FIXED_ROWS:nrow(testDs_j)){
                        curVisitTime = testDs_j$visitTimeYears[j]
                        
                        if((curVisitTime - lastBiopsyTime) >= 1){
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
                      
                      return(c("nb"=nb, "offset"=offset))
                    }
  stopCluster(cl)
  
  rownames(results) = NULL
  
  return(results)
}

#####################################################
# Hybrid schedule number 2
# If you can fail in next 1 year use 15% rule, else check if the last biopsy is more than 8 years. if it is after 8 years use mean, if it is before 8 years use median
# we use Median failure time if the last biopsy is before 8 years and the if ROC tells us he can fail in 8 years
# otherwise we use mean failure time
#####################################################
runHybridSchedule2 = function(jointModelData){
  
  max_surv_time = max(jointModelData$trainingData$trainingDs$progression_time) - 10e-2
  
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
  
  testDs.id = jointModelData$testData$testDs.id
  testDs = jointModelData$testData$testDs[jointModelData$testData$testDs$visitTimeYears <= max_surv_time,]
  
  cl = makeCluster(max_cores)
  registerDoParallel(cl)
  
  #testDs.id$P_ID
  results = foreach(pid = testDs.id$P_ID, .combine = "rbind", .packages = c("splines", "JMbayes"),
                    .export=c("MIN_FIXED_ROWS", "MAX_FAIL_TIME", "thresholdsList")) %dopar%{
                      progressionTime = testDs.id$progression_time[testDs.id$P_ID == pid]
                      
                      lastBiopsyTime = 0
                      nb = 0
                      offset = lastBiopsyTime - progressionTime
                      
                      testDs_j = testDs[testDs$P_ID == pid,]
                      
                      #the first row is time 0
                      for(j in MIN_FIXED_ROWS:nrow(testDs_j)){
                        curVisitTime = testDs_j$visitTimeYears[j]
                        if((curVisitTime - lastBiopsyTime) >= 1){
                          if(lastBiopsyTime < 9){
                            Dt = 1
                            cutpoints = thresholdsList[[as.character(Dt)]]$cut_points
                            
                            nearestTimeIndex = which.min(abs(lastBiopsyTime - cutpoints[,1]))
                            
                            survFit_j = survfitJM(object=jointModelData$mvJoint_dre_psa_simDs, 
                                                  newdata = testDs_j[1:j,], idVar = "P_ID",
                                                  survTimes = lastBiopsyTime + Dt, last.time = lastBiopsyTime)
                            
                            meanSurvProb = round(survFit_j$summaries[[1]][, "Mean"], 2)
                            
                            if(meanSurvProb <= cutpoints[nearestTimeIndex, 2]){
                              sfit2 = survfitJM(object=jointModelData$mvJoint_dre_psa_simDs, 
                                                newdata = testDs_j[1:j,], idVar = "P_ID",
                                                survTimes = curVisitTime, last.time = lastBiopsyTime)
                              
                              meanSurv2 = round(sfit2$summaries[[1]][, "Mean"], 2)
                              proposedBiopsyTime = ifelse(meanSurv2>0.85, Inf, curVisitTime)
                            }else{
                              proposedBiopsyTime = pMedianTimeOptimal(testDs_j[1:j,], lastBiopsyTime)
                            }
                          }else{
                            proposedBiopsyTime = expectedCondFailureTime(newdata = testDs_j[1:j,], lastBiopsyTime = lastBiopsyTime)
                          }

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
                      
                      return(c("nb"=nb, "offset"=offset))
                    }
  stopCluster(cl)
  
  rownames(results) = NULL
  
  return(results)
}
