#Calculate the prevalence of the disease
computeRoc=function(dsId, Dt = 1){
  ct = makeCluster(4)
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

computeCutOffValues = function(dsId){
  
  rocList = simulatedDsList[[dsId]]$rocList_2
  
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

computeBiopsyTimes = function(minVisits = 5, dsId, patientRowNum){
  
  cutoffValues = simulatedDsList[[dsId]]$cutoffValues
  
  testDs = simulatedDsList[[dsId]]$testDs
  patientDs_i = split(testDs, testDs$P_ID)[[patientRowNum]]
  
  patientDs_i$expectedFailureTime = NA
  patientDs_i$survTimeYouden = NA
  patientDs_i$survTimeAccuracy = NA
  patientDs_i$survTimeMaxTPR = NA
  patientDs_i$survTimeMinFPR = NA
  patientDs_i$survTimeF1Score = NA
  patientDs_i$survTime85 = NA
  
  #instead of times per subject I choose 28 because 28th time is 11 years.
  #expected time of failure is not available for last time point if faiure time for all subejccts is less than tht last time point
  #training max progression time is 11.216 and test is 10. something
  visitsOfInterest = minVisits:30
  res = foreach(j=visitsOfInterest, .packages = c("splines", "JMbayes", "coda"),
          .export=c("timesPerSubject", "dynamicCutOffTimes",
                    "expectedCondFailureTime", "dynamicPredProb",
                    "simulatedDsList", "pDynSurvTime", "invDynSurvival")) %dopar%{
                      persTestDs = patientDs_i[patientDs_i$visitNumber <= j,]
                      temp_lasttime = max(persTestDs$visitTimeYears)
                      nearest_time_index = which(abs(dynamicCutOffTimes-temp_lasttime)==min(abs(dynamicCutOffTimes-temp_lasttime)))
                      nearest_time_index = nearest_time_index[1]
                      
                      expectedFailureTime = expectedCondFailureTime(dsId, persTestDs)
                      survTimeYouden = pDynSurvTime(survProb = cutoffValues[[nearest_time_index]]["youden"], dsId, persTestDs)
                      survTimeAccuracy = pDynSurvTime(survProb =  cutoffValues[[nearest_time_index]]["accuracy"], dsId, persTestDs)
                      survTimeMaxTPR = pDynSurvTime(survProb =  cutoffValues[[nearest_time_index]]["maxTPR"], dsId, persTestDs)
                      survTimeF1Score = pDynSurvTime(survProb =  cutoffValues[[nearest_time_index]]["f1score"], dsId, persTestDs)
                      survTime85 = pDynSurvTime(survProb = 0.85, dsId, persTestDs)
                      
                      list(expectedFailureTime, survTimeYouden, survTimeAccuracy, 
                           survTimeMaxTPR, survTimeF1Score, survTime85)
                    }
  
  patientDs_i$expectedFailureTime[visitsOfInterest] = sapply(res, function(x){x[[1]]}, simplify = T)
  patientDs_i$survTimeYouden[visitsOfInterest] = sapply(res, function(x){x[[2]]}, simplify = T)
  patientDs_i$survTimeAccuracy[visitsOfInterest] = sapply(res, function(x){x[[3]]}, simplify = T)
  patientDs_i$survTimeMaxTPR[visitsOfInterest] = sapply(res, function(x){x[[4]]}, simplify = T)
  patientDs_i$survTimeF1Score[visitsOfInterest] = sapply(res, function(x){x[[5]]}, simplify = T)
  patientDs_i$survTime85[visitsOfInterest] = sapply(res, function(x){x[[6]]}, simplify = T)
  
  return(patientDs_i)
}

computeBiopsyProbs = function(minVisits = 5, dsId, patientRowNum){
  
  cutoffValues = simulatedDsList[[dsId]]$cutoffValues_2 
  
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

summarizeBiopsyResults =  function(dsId, patientRowNum, minVisits, biopsyIfLessThanTime=1, 
                                   methodName = "expectedFailureTime", biopsyEveryKYears = NA){
  trueProgressionTime = simulatedDsList[[dsId]]$testDs.id[patientRowNum,]$progression_time
     
  visitTimeYears = simulatedDsList[[dsId]]$biopsyTimes[[patientRowNum]]$visitTimeYears
  biopsyTimes = simulatedDsList[[dsId]]$biopsyTimes[[patientRowNum]][, methodName]
  
  nb = 0
  biopsyTimeOffset = NA
  
  lastBiopsyTime = -Inf

  for(j in minVisits:timesPerSubject){
    curVisitTime = visitTimeYears[j]
    predBiopsyTime = biopsyTimes[j]
    
    biopsy_gap_needed = (curVisitTime - lastBiopsyTime) < 1
    
    condition1 = !is.na(predBiopsyTime) & ((predBiopsyTime - curVisitTime) <= biopsyIfLessThanTime)
    condition2 = !is.na(biopsyEveryKYears) & !is.na(predBiopsyTime) & ((predBiopsyTime - lastBiopsyTime) >= biopsyEveryKYears)
    if(biopsy_gap_needed==FALSE & (condition1 | condition2)){
      biopsyTimesOfInterest = biopsyTimes[visitTimeYears >= curVisitTime & visitTimeYears<=predBiopsyTime]
      biopsyIndexOfInterest = which.min(biopsyTimesOfInterest) - 1 + j
      
      nb = nb + 1
      biopsyTimeOffset = biopsyTimes[biopsyIndexOfInterest] - trueProgressionTime
      lastBiopsyTime = biopsyTimes[biopsyIndexOfInterest]
      
      if(biopsyTimeOffset >= 0){
        break
      }
    }
  }
  
  return(c(patientRowNum = patientRowNum, nb = nb, biopsyTimeOffset = biopsyTimeOffset))
}

getBiopsyResults = function(dsId, biopsyIfLessThanTime = 1, biopsyEveryKYears = NA , minVisits, methodNames = c("expectedFailureTime", "survTime85",
                                                  "survTimeYouden", "survTimeAccuracy", "survTimeMaxTPR", "survTimeF1Score")){
  subjectCount = length(simulatedDsList[[dsId]]$biopsyTimes)
  
  biopsyResults = data.frame(patientRowNum = numeric(), nb = numeric(), biopsyTimeOffset = numeric(), methodName=factor(levels = methodNames))
  for(methodName in methodNames){
    summary = sapply(1:subjectCount, summarizeBiopsyResults, dsId = dsId, methodName = methodName, minVisits = minVisits, biopsyIfLessThanTime=biopsyIfLessThanTime)  
    biopsyResults = rbind(biopsyResults, data.frame(t(summary), methodName = rep(methodName,ncol(summary))))
  }
  
  if(!any(is.na(biopsyEveryKYears))){
    for(K in biopsyEveryKYears){
      for(methodName in methodNames){
        extMethodName = paste(methodName, K, sep="")
        summary = sapply(1:subjectCount, summarizeBiopsyResults, dsId = dsId, methodName = methodName, minVisits = minVisits, biopsyIfLessThanTime=biopsyIfLessThanTime, biopsyEveryKYears = K)
        biopsyResults = rbind(biopsyResults, data.frame(t(summary), methodName = rep(extMethodName,ncol(summary))))
      }
    }
  }
  
  #Fixed schedule
  fixedSummary = sapply(1:subjectCount, function(patientRowNum){
    progressionTime = simulatedDsList[[dsId]]$testDs.id$progression_time[patientRowNum]
    
    nb = 0
    if(longTimes[minVisits]>=1){
      nb = -1
    }
    
    if(progressionTime<=1){
      return(c(patientRowNum = patientRowNum, nb = nb+1, biopsyTimeOffset = 1-progressionTime))
    }else if(progressionTime<=4){
      return(c(patientRowNum = patientRowNum, nb = nb+2, biopsyTimeOffset = 4-progressionTime))
    }else if(progressionTime<=7){
      return(c(patientRowNum = patientRowNum, nb = nb+3, biopsyTimeOffset = 7-progressionTime))
    }else if(progressionTime<=10){
      return(c(patientRowNum = patientRowNum, nb = nb+4, biopsyTimeOffset = 10-progressionTime))
    }else if(progressionTime<=15){
      return(c(patientRowNum = patientRowNum, nb = nb+5, biopsyTimeOffset = 15-progressionTime))
    }else if(progressionTime<=20){
      return(c(patientRowNum = patientRowNum, nb = nb+6, biopsyTimeOffset = 20-progressionTime))
    }
  }, simplify = T)
  biopsyResults = rbind(biopsyResults, data.frame(t(fixedSummary), methodName = rep("fixed",ncol(fixedSummary))))
  
  #Johns hopkins schedule
  johnsHopkinsSummary =  sapply(1:subjectCount, function(patientRowNum){
    progressionTime = simulatedDsList[[dsId]]$testDs.id$progression_time[patientRowNum]
    detectionTime = ceiling(progressionTime)
    return(c(patientRowNum = patientRowNum, nb = detectionTime-1, biopsyTimeOffset = detectionTime-progressionTime))
  })
  biopsyResults = rbind(biopsyResults, data.frame(t(johnsHopkinsSummary), methodName = rep("johnsSummary",ncol(johnsHopkinsSummary))))
  
  biopsyResults
}



# approach2res = approach2res[!is.na(approach2res[,5]) & !is.na(approach2res[,6]) & !is.na(approach2res[,7]) & approach2res[,5]>0 & approach2res[,6]>0 & approach2res[,7]>0,]
# apply(approach2res, MARGIN = 2, mean, na.rm=T)
