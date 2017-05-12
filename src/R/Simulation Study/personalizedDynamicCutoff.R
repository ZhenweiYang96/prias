#Calculate the prevalence of the disease
computeRoc=function(dsId, Dt = 1){
  ct = makeCluster(cores)
  registerDoParallel(ct)
  
  rocList = foreach(tstart = simulatedDsList[[dsId]]$dynamicCutOffTimes, .packages = c("splines", "JMbayes"),
                    .export=c("simulatedDsList", "rocJM_mod")) %dopar%{
                      res <- tryCatch({
                        rocJM_mod(simulatedDsList[[dsId]]$models$simJointModel_replaced, simulatedDsList[[dsId]]$trainingDs, 
                              Tstart=tstart, Dt=Dt, idVar = "P_ID")
                      }, error=function(e) NULL)
                    }
  stopCluster(ct)
  return(rocList)
}

computeRocDataDrivenDt=function(dsId){
  ct = makeCluster(cores)
  registerDoParallel(ct)
  
  rocList = foreach(tstart = simulatedDsList[[dsId]]$dynamicCutOffTimes, .packages = c("splines", "JMbayes"),
                    .export=c("simulatedDsList", "rocJM_mod")) %dopar%{
                      res <- tryCatch({
                        aucDt1 = aucJM(simulatedDsList[[dsId]]$models$simJointModel_replaced, simulatedDsList[[dsId]]$trainingDs, 
                                  Tstart=tstart, Dt=1, idVar = "P_ID")$auc
                        
                        aucDtpt5 = aucJM(simulatedDsList[[dsId]]$models$simJointModel_replaced, simulatedDsList[[dsId]]$trainingDs, 
                                       Tstart=tstart, Dt=0.5, idVar = "P_ID")$auc
                        
                        Dt = 1
                        if(aucDt1 < aucDtpt5){
                          Dt = 0.5
                        }
                        
                        rocJM_mod(simulatedDsList[[dsId]]$models$simJointModel_replaced, simulatedDsList[[dsId]]$trainingDs, 
                                  Tstart=tstart, Dt=Dt, idVar = "P_ID")
                      }, error=function(e) NULL)
                    }
  stopCluster(ct)
  return(rocList)
}

computeCutOffValues = function(dsId){
  
  rocList = simulatedDsList[[dsId]]$rocList
  
  cutoffValues = lapply(1:length(rocList), function(i){
    
    #For youden and sens
    youden = maxTPR = minFPR = minFNR = maxTNR = NA
    k=0
    while((i-k)>0){
      rocRes = rocList[[i-k]]
      
      if(is.null(rocRes) | any(is.nan(rocRes$TP)) | any(is.nan(rocRes$FP))){
        k = k + 1
      }else{
        youden = rocRes$thrs[which.max(rocRes$TP - rocRes$FP)]
        maxTPR = rocRes$thrs[which.max(rocRes$TP)]
        minFPR = rocRes$thrs[which.min(rocRes$FP)]
        minFNR = rocRes$thrs[which.min(1-rocRes$TP)]
        maxTNR = rocRes$thrs[which.max(1-rocRes$FP)]
        break
      }
    }
    
    #For youden and sens
    accuracy = f1score = mcc = NA
    k=0
    while((i-k)>0){
      rocRes = rocList[[i-k]]
      
      if(is.null(rocRes) | any(is.nan(rocRes$nTP)) | any(is.nan(rocRes$nFP)) | any(is.nan(rocRes$nFN)) | any(is.nan(rocRes$nTN))){
        k = k + 1
      }else{
        accuracy = rocRes$thrs[which.max((rocRes$nTP + rocRes$nTN)/(rocRes$nTP + rocRes$nTN + rocRes$nFN + rocRes$nFP))]
        f1score = rocRes$thrs[which.max(2*rocRes$nTP/(2*rocRes$nTP + rocRes$nFN + rocRes$nFP))]
        mcc = rocRes$thrs[which.max((rocRes$nTP * rocRes$nTN - rocRes$nFP * rocRes$nFN)/((rocRes$nTP + rocRes$nFP)*(rocRes$nTP + rocRes$nFN)*(rocRes$nTN + rocRes$nFP)*(rocRes$nTN + rocRes$nFN)))]
        break
      }
    }
    
    kmfit = survfit(Surv(progression_time, progressed)~1, conf.type="log-log", data=prias.id)
    kmprob <- stepfun(kmfit$time, c(1, kmfit$surv))
    prevalence = kmprob(simulatedDsList[[dsId]]$dynamicCutOffTimes) - kmprob(simulatedDsList[[dsId]]$dynamicCutOffTimes + 1)

    #For markedness and npv
    markedness = maxNPV = NA
    k=0
    while((i-k)>0){
      rocRes = rocList[[i-k]]
      
      if(!is.null(rocRes)){
        ppv = (rocRes$TP * prevalence[i])/(rocRes$TP * prevalence[i] + (1-prevalence[i]) * rocRes$FP)
        npv = ((1-prevalence[i]) * (1-rocRes$FP))/((1-prevalence[i]) * (1-rocRes$FP) + prevalence[i] * (1-rocRes$TP))
        
        if(all(is.nan(ppv + npv)) | all(is.nan(npv))){
           k = k + 1
        }else{
          markedness = rocRes$thrs[which.max(ppv + npv - 1)]
          maxNPV = rocRes$thrs[which.max(npv)]
          break
        }
      }else{
        k = k + 1
      }
    }
    
    c(youden = youden, maxTPR = maxTPR, minFPR = minFPR, minFNR = minFNR, maxTNR = maxTNR, accuracy=accuracy, f1score = f1score, mcc=mcc, markedness = markedness, maxNPV = maxNPV)
  })
}

computeBiopsyTimes = function(dsId, visitCount){
  
  testDs = simulatedDsList[[dsId]]$testDs
  testDs = testDs[testDs$visitNumber<=visitCount,]
  patientDsList = split(testDs, testDs$P_ID)
  
  dynamicCutOffTimes = simulatedDsList[[dsId]]$dynamicCutOffTimes
  
  #All patients share the last visit time and cutoff values here. So I randomly choose one of the patients
  lastVisitTime = max(patientDsList[[1]]$visitTimeYears)
  nearest_time_index = which(abs(dynamicCutOffTimes-lastVisitTime)==min(abs(dynamicCutOffTimes-lastVisitTime)))[1]
  cutoffValues = c(0.85, simulatedDsList[[dsId]]$cutoffValues[[nearest_time_index]][c("youden", "maxTPR", "accuracy", "f1score")])
    
  round1Granularity = 1
  if(15-lastVisitTime < round1Granularity){
    round1Granularity = 15-lastVisitTime
  }
  
  biopsyTimes = foreach(i=1:length(patientDsList), .packages = c("splines", "JMbayes", "coda"),
          .export=c("timesPerSubject", "dynamicCutOffTimes",
                    "expectedCondFailureTime", "dynamicPredProb",
                    "simulatedDsList", "pDynSurvTime", "invDynSurvival")) %dopar%{
                      patientDs_i = patientDsList[[i]]
                      
                      round1 = survfitJM(simulatedDsList[[dsId]]$models$simJointModel_replaced, 
                                         newdata = patientDs_i, idVar = "P_ID", 
                                         survTimes = seq(lastVisitTime, 15, round1Granularity))
                      
                      times = round1$summaries[[1]][, "times"]
                      probs = round1$summaries[[1]][, "Median"]
                      
                      newTimes = do.call(c, lapply(cutoffValues[!is.na(cutoffValues)], function(p){
                        nearest_time = times[which(abs(probs-p)==min(abs(probs-p)))[1]]
                        upper = nearest_time + 1
                        lower = nearest_time - 1
                        if(upper > 15){
                          upper = 15
                        }
                        if(lower < 0){
                          lower = 0
                        }
                        seq(lower, upper, 0.05)
                      }))
                      
                      round2 = survfitJM(simulatedDsList[[dsId]]$models$simJointModel_replaced, 
                                         newdata = patientDs_i, idVar = "P_ID", 
                                         survTimes = unique(newTimes))
                      
                      times = round2$summaries[[1]][, "times"]
                      probs = round2$summaries[[1]][, "Median"]
                      
                      c(expectedCondFailureTime(dsId, patientDs_i), unlist(lapply(cutoffValues, function(prob){
                        if(is.na(prob)){
                          NA
                        }else{
                          times[which(abs(probs-prob)==min(abs(probs-prob)))[1]]
                        }
                      })))
                    }
  return(biopsyTimes)
  
}


# computeBiopsyTimes = function(minVisits = 5, dsId, patientRowNum){
#   
#   cutoffValues = simulatedDsList[[dsId]]$cutoffValues
#   
#   testDs = simulatedDsList[[dsId]]$testDs
#   patientDs_i = split(testDs, testDs$P_ID)[[patientRowNum]]
#   
#   patientDs_i$expectedFailureTime = NA
#   patientDs_i$survTimeYouden = NA
#   patientDs_i$survTimeAccuracy = NA
#   patientDs_i$survTimeMaxTPR = NA
#   patientDs_i$survTimeMinFPR = NA
#   patientDs_i$survTimeF1Score = NA
#   patientDs_i$survTime85 = NA
#   
#   #instead of times per subject I choose 28 because 28th time is 11 years.
#   #expected time of failure is not available for last time point if faiure time for all subejccts is less than tht last time point
#   #training max progression time is 11.216 and test is 10. something
#   visitsOfInterest = minVisits:30
#   res = foreach(j=visitsOfInterest, .packages = c("splines", "JMbayes", "coda"),
#           .export=c("timesPerSubject", "dynamicCutOffTimes",
#                     "expectedCondFailureTime", "dynamicPredProb",
#                     "simulatedDsList", "pDynSurvTime", "invDynSurvival")) %dopar%{
#                       persTestDs = patientDs_i[patientDs_i$visitNumber <= j,]
#                       temp_lasttime = max(persTestDs$visitTimeYears)
#                       nearest_time_index = which(abs(dynamicCutOffTimes-temp_lasttime)==min(abs(dynamicCutOffTimes-temp_lasttime)))
#                       nearest_time_index = nearest_time_index[1]
#                       
#                       expectedFailureTime = expectedCondFailureTime(dsId, persTestDs)
#                       survTimeYouden = pDynSurvTime(survProb = cutoffValues[[nearest_time_index]]["youden"], dsId, persTestDs)
#                       survTimeAccuracy = pDynSurvTime(survProb =  cutoffValues[[nearest_time_index]]["accuracy"], dsId, persTestDs)
#                       survTimeMaxTPR = pDynSurvTime(survProb =  cutoffValues[[nearest_time_index]]["maxTPR"], dsId, persTestDs)
#                       survTimeF1Score = pDynSurvTime(survProb =  cutoffValues[[nearest_time_index]]["f1score"], dsId, persTestDs)
#                       survTime85 = pDynSurvTime(survProb = 0.85, dsId, persTestDs)
#                       
#                       list(expectedFailureTime, survTimeYouden, survTimeAccuracy, 
#                            survTimeMaxTPR, survTimeF1Score, survTime85)
#                     }
#   
#   patientDs_i$expectedFailureTime[visitsOfInterest] = sapply(res, function(x){x[[1]]}, simplify = T)
#   patientDs_i$survTimeYouden[visitsOfInterest] = sapply(res, function(x){x[[2]]}, simplify = T)
#   patientDs_i$survTimeAccuracy[visitsOfInterest] = sapply(res, function(x){x[[3]]}, simplify = T)
#   patientDs_i$survTimeMaxTPR[visitsOfInterest] = sapply(res, function(x){x[[4]]}, simplify = T)
#   patientDs_i$survTimeF1Score[visitsOfInterest] = sapply(res, function(x){x[[5]]}, simplify = T)
#   patientDs_i$survTime85[visitsOfInterest] = sapply(res, function(x){x[[6]]}, simplify = T)
#   
#   return(patientDs_i)
# }

computeBiopsyProbs = function(minVisits = 5, dsId, patientRowNum){
  
  cutoffValues = simulatedDsList[[dsId]]$cutoffValues
  
  testDs = simulatedDsList[[dsId]]$testDs
  patientDs_i = split(testDs, testDs$P_ID)[[patientRowNum]]
  
  patientDs_i$survProbYouden = NA
  patientDs_i$survProbAccuracy = NA
  patientDs_i$survProbF1Score = NA
  patientDs_i$survProbmaxTPR = NA
  patientDs_i$survProb85 = NA
  
  #instead of times per subject I choose 28 because 28th time is 11 years.
  #expected time of failure is not available for last time point if faiure time for all subejccts is less than tht last time point
  #training max progression time is 11.216 and test is 10. something
  visitsOfInterest = minVisits:30
  res = foreach(j=visitsOfInterest, .packages = c("splines", "JMbayes", "coda"),
                .export=c("timesPerSubject", "dynamicCutOffTimes", "simulatedDsList")) %do%{
                  persTestDs = patientDs_i[patientDs_i$visitNumber <= j,]
                  temp_lasttime = max(persTestDs$visitTimeYears)
                  nearest_time_index = which(abs(dynamicCutOffTimes-temp_lasttime)==min(abs(dynamicCutOffTimes-temp_lasttime)))
                  nearest_time_index = nearest_time_index[1]
                  
                  survProbYouden = cutoffValues[[nearest_time_index]]["youden"]
                  survProbAccuracy = cutoffValues[[nearest_time_index]]["accuracy"]
                  survProbF1Score = cutoffValues[[nearest_time_index]]["f1score"]
                  survProbmaxTPR = cutoffValues[[nearest_time_index]]["maxTPR"]
                  survProb85 = 0.85
                  
                  list(survProbYouden, survProbAccuracy, 
                       survProbF1Score, survProbmaxTPR, survProb85)
                }
  
  patientDs_i$survProbYouden[visitsOfInterest] = sapply(res, function(x){x[[1]]}, simplify = T)
  patientDs_i$survProbAccuracy[visitsOfInterest] = sapply(res, function(x){x[[2]]}, simplify = T)
  patientDs_i$survProbF1Score[visitsOfInterest] = sapply(res, function(x){x[[3]]}, simplify = T)
  patientDs_i$survProbmaxTPR[visitsOfInterest] = sapply(res, function(x){x[[4]]}, simplify = T)
  patientDs_i$survProb85[visitsOfInterest] = sapply(res, function(x){x[[5]]}, simplify = T)
  
  return(patientDs_i)
}