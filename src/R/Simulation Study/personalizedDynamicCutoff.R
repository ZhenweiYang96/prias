dynamicCutOffTimes = seq(generateLongtiudinalTimeBySchedule()[minVisits], 15, 0.1)
#Calculate the prevalence of the disease
kmfit = survfit(Surv(progression_time, progressed)~1, conf.type="log-log", data=prias.id)
kmprob <- stepfun(kmfit$time, c(1, kmfit$surv))
prevalence = kmprob(dynamicCutOffTimes) - kmprob(dynamicCutOffTimes + 1)

computeCutOffValues = function(dsId){
  
  rocList = foreach(tstart = dynamicCutOffTimes, .packages = c("splines", "JMbayes"),
                    .export=c("simulatedDsList")) %dopar%{
    res <- tryCatch({
      rocJM(simulatedDsList[[dsId]]$models$simJointModel_replaced, simulatedDsList[[dsId]]$trainingDs, 
            Tstart=tstart, Dt=1, idVar = "P_ID")
    }, error=function(e) NULL)
  }
  
  cutoffValues = lapply(1:length(rocList), function(i){
    
    #For youden and sens
    youden = maxSens = NA
    k=0
    while((i-k)>0){
      rocRes = rocList[[i-k]]
      
      if(is.null(rocRes) | any(is.nan(rocRes$TP)) | any(is.nan(rocRes$FP))){
        k = k + 1
      }else{
        youden = rocRes$thrs[which.max(rocRes$TP - rocRes$FP)]
        maxSens = rocRes$thrs[which.max(rocRes$TP)]
        break
      }
    }
    
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
    
    c(youden = youden, maxSens = maxSens, markedness = markedness, maxNPV = maxNPV)
  })
}

computeBiopsyTimes = function(minVisits = 5, dsId, patientRowNum){
  
  cutoffValues = simulatedDsList[[dsId]]$cutoffValues
  
  testDs = simulatedDsList[[dsId]]$testDs
  patientDs_i = split(testDs, testDs$P_ID)[[patientRowNum]]
  
  patientDs_i$expectedFailureTime = NA
  patientDs_i$survTimeYouden = NA
  patientDs_i$survTimeMarkedness = NA
  patientDs_i$survTimeMaxSens = NA
  patientDs_i$survTimeMaxNPV = NA
  patientDs_i$survTime90 = NA
  patientDs_i$survTime80 = NA
  
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
                      survTimeMarkedness = pDynSurvTime(survProb =  cutoffValues[[nearest_time_index]]["markedness"], dsId, persTestDs)
                      survTimeMaxSens = pDynSurvTime(survProb =  cutoffValues[[nearest_time_index]]["maxSens"], dsId, persTestDs)
                      survTimeMaxNPV = pDynSurvTime(survProb =  cutoffValues[[nearest_time_index]]["maxNPV"], dsId, persTestDs)
                      survTime90 = pDynSurvTime(survProb = 0.9, dsId, persTestDs)
                      survTime80 = pDynSurvTime(survProb = 0.8, dsId, persTestDs)
                      
                      list(expectedFailureTime, survTimeYouden, survTimeMarkedness, survTimeMaxSens, survTimeMaxNPV, survTime90, survTime80)
                    }
  
  patientDs_i$expectedFailureTime[visitsOfInterest] = sapply(res, function(x){x[[1]]}, simplify = T)
  patientDs_i$survTimeYouden[visitsOfInterest] = sapply(res, function(x){x[[2]]}, simplify = T)
  patientDs_i$survTimeMarkedness[visitsOfInterest] = sapply(res, function(x){x[[3]]}, simplify = T)
  patientDs_i$survTimeMaxSens[visitsOfInterest] = sapply(res, function(x){x[[4]]}, simplify = T)
  patientDs_i$survTimeMaxNPV[visitsOfInterest] = sapply(res, function(x){x[[5]]}, simplify = T)
  patientDs_i$survTime90[visitsOfInterest] = sapply(res, function(x){x[[6]]}, simplify = T)
  patientDs_i$survTime80[visitsOfInterest] = sapply(res, function(x){x[[7]]}, simplify = T)
  
  return(patientDs_i)
}

# biopsy_times = list(expectedFailureTime = c(), survTimeYouden = c(), survTimeMarkedness = c(),
#                     survTimeMaxSens = c(), survTimeMaxNPV=c(), survTime90 = c(), survTime80=c())
# #Expected failure time
# if(patientDs_i$expectedFailureTime[j] < patientDs_i$visitTimeYears[j + 1]){
#   biopsy_times$expectedFailureTime = c(biopsy_times$expectedFailureTime, patientDs_i$expectedFailureTime[j])
# }
# 
# #Youden
# if(patientDs_i$survTimeYouden[j] < patientDs_i$visitTimeYears[j + 1]){
#   biopsy_times$survTimeYouden = c(biopsy_times$survTimeYouden, patientDs_i$survTimeYouden[j])
# }
# 
# #Markedness
# if(patientDs_i$survTimeMarkedness[j] < patientDs_i$visitTimeYears[j + 1]){
#   biopsy_times$survTimeMarkedness = c(biopsy_times$survTimeMarkedness, patientDs_i$survTimeMarkedness[j])
# }
# 
# #Max sensitivity
# if(patientDs_i$survTimeMaxSens[j] < patientDs_i$visitTimeYears[j + 1]){
#   biopsy_times$survTimeMaxSens = c(biopsy_times$survTimeMaxSens, patientDs_i$survTimeMaxSens[j])
# }
# 
# #Max NPV
# if(patientDs_i$survTimeMaxNPV[j] < patientDs_i$visitTimeYears[j + 1]){
#   biopsy_times$survTimeMaxNPV = c(biopsy_times$survTimeMaxNPV, patientDs_i$survTimeMaxNPV[j])
# }
# 
# #90% survival
# if(patientDs_i$survTime90[j] < patientDs_i$visitTimeYears[j + 1]){
#   biopsy_times$survTime90 = c(biopsy_times$survTime90, patientDs_i$survTime90[j])
# }
# 
# #80% survival
# if(patientDs_i$survTime80[j] < patientDs_i$visitTimeYears[j + 1]){
#   biopsy_times$survTime80 = c(biopsy_times$survTime80, patientDs_i$survTime80[j])
# }

# 
# #Store original results as they are
# orig_res = res
# 
# ###########################################################
# #Approach #1 do biopsy as given by the estimates from the above model
# ############################################################
# approach1res = foreach(i=1:length(res), .combine = "rbind")%do%{
#   trueProgressionTime = testDs.id[i,]$progression_time
#   
#   biopsy_times = res[[i]][[2]]
#   #Expected value
#   nb1 = which(biopsy_times$expectedFailureTime >= trueProgressionTime)[1]
#   nb2 = which(biopsy_times$survTimeYouden >= trueProgressionTime)[1]
#   nb3 = which(biopsy_times$survTimeMarkedness >= trueProgressionTime)[1]
#   
#   biopsytimeOffset1 = biopsytimeOffset2 = biopsytimeOffset3 = NA
#   
#   if(!is.na(nb1)){
#     biopsytimeOffset1 = biopsy_times$expectedFailureTime[nb1] - trueProgressionTime
#   }
#   
#   if(!is.na(nb2)){
#     biopsytimeOffset2 = biopsy_times$survTimeYouden[nb2] - trueProgressionTime
#   }
#   
#   if(!is.na(nb3)){
#     biopsytimeOffset3 = biopsy_times$survTimeMarkedness[nb2] - trueProgressionTime
#   }
#   
#   c(nb1, nb2, nb3, biopsytimeOffset1, biopsytimeOffset2, biopsytimeOffset3)
# }
# 
# #############################################################
# #Approach #2 do biopsy if the failure time is less than an year then do biopsy at that proposed time
# ############################################################

summarizeBiopsyResults =  function(dsId, patientRowNum, methodName = "expectedFailureTime", biopsyEvery3Years = F){
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
    
    condition1 = !is.na(predBiopsyTime) & ((predBiopsyTime - curVisitTime) <= 1)
    condition2 = biopsyEvery3Years & !is.na(predBiopsyTime) & ((predBiopsyTime - lastBiopsyTime) >= 3)
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

getBiopsyResults = function(dsId, methodNames = c("expectedFailureTime", "survTime90", "survTime80", 
                                                  "survTimeYouden", "survTimeMarkedness", "survTimeMaxSens", "survTimeMaxNPV")){
  subjectCount = length(simulatedDsList[[dsId]]$biopsyTimes)
  
  biopsyResults = data.frame(patientRowNum = numeric(), nb = numeric(), biopsyTimeOffset = numeric(), methodName=factor(levels = methodNames))
  for(methodName in methodNames){
    summary = sapply(1:subjectCount, summarizeBiopsyResults, dsId = dsId, methodName = methodName)  
    biopsyResults = rbind(biopsyResults, data.frame(t(summary), methodName = rep(methodName,ncol(summary))))
  }
  
  for(methodName in methodNames){
    extMethodName = paste(methodName, 3, sep="")
    summary = sapply(1:subjectCount, summarizeBiopsyResults, dsId = dsId, methodName = methodName, biopsyEvery3Years = T)
    biopsyResults = rbind(biopsyResults, data.frame(t(summary), methodName = rep(extMethodName,ncol(summary))))
  }
  
  #Fixed schedule
  fixedSummary = sapply(1:subjectCount, function(patientRowNum){
    progressionTime = simulatedDsList[[dsId]]$testDs.id$progression_time[patientRowNum]
    if(progressionTime<=4){
      return(c(patientRowNum = patientRowNum, nb = 1, biopsyTimeOffset = 4-progressionTime))
    }else if(progressionTime<=7){
      return(c(patientRowNum = patientRowNum, nb = 2, biopsyTimeOffset = 7-progressionTime))
    }else if(progressionTime<=10){
      return(c(patientRowNum = patientRowNum, nb = 3, biopsyTimeOffset = 10-progressionTime))
    }else if(progressionTime<=15){
      return(c(patientRowNum = patientRowNum, nb = 4, biopsyTimeOffset = 15-progressionTime))
    }
  })
  biopsyResults = rbind(biopsyResults, data.frame(t(fixedSummary), methodName = rep("fixed",ncol(summary))))
  
  biopsyResults
}



# approach2res = approach2res[!is.na(approach2res[,5]) & !is.na(approach2res[,6]) & !is.na(approach2res[,7]) & approach2res[,5]>0 & approach2res[,6]>0 & approach2res[,7]>0,]
# apply(approach2res, MARGIN = 2, mean, na.rm=T)
