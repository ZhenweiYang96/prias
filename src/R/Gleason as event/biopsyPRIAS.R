#So instead of writing the code from scratch, I will try to use the existing stuff
dynamicCutoffTimes_PRIAS = seq(1/12, max(prias.id$progression_time)-0.0001, 0.1)
ct = makeCluster(detectCores())
registerDoParallel(ct)

rocList_PRIAS = foreach(tstart = dynamicCutoffTimes_PRIAS, .packages = c("splines", "JMbayes"),
                  .export=c("rocJM_mod")) %dopar%{
                    res <- tryCatch({
                      rocJM_mod(joint_psa_replaced, psa_data_set,
                                Tstart=tstart, Dt=1, idVar = "P_ID")
                    }, error=function(e) NULL)
                  }
stopCluster(ct)

#################################################################################
cutoffValues_PRIAS = lapply(1:length(rocList_PRIAS), function(i){

  #For youden and sens
  youden = maxTPR = minFPR = minFNR = maxTNR = NA
  k=0
  while((i-k)>0){
    rocRes = rocList_PRIAS[[i-k]]

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
    rocRes = rocList_PRIAS[[i-k]]

    if(is.null(rocRes) | any(is.nan(rocRes$nTP)) | any(is.nan(rocRes$nFP)) | any(is.nan(rocRes$nFN)) | any(is.nan(rocRes$nTN))){
      k = k + 1
    }else{
      accuracy = rocRes$thrs[which.max((rocRes$nTP + rocRes$nTN)/(rocRes$nTP + rocRes$nTN + rocRes$nFN + rocRes$nFP))]
      f1score = rocRes$thrs[which.max(2*rocRes$nTP/(2*rocRes$nTP + rocRes$nFN + rocRes$nFP))]
      mcc = rocRes$thrs[which.max((rocRes$nTP * rocRes$nTN - rocRes$nFP * rocRes$nFN)/((rocRes$nTP + rocRes$nFP)*(rocRes$nTP + rocRes$nFN)*(rocRes$nTN + rocRes$nFP)*(rocRes$nTN + rocRes$nFN)))]
      break
    }
  }

  c(youden = youden, maxTPR = maxTPR, minFPR = minFPR, minFNR = minFNR, maxTNR = maxTNR, accuracy=accuracy, f1score = f1score, mcc=mcc)
})

#How to search for the right patient profiles
demo_pids = prias.id$P_ID[order(by(psa_data_set$logpsa1, psa_data_set$P_ID, max), decreasing = T)][50:1]
for(demoPid in demo_pids){
  dataset = psa_data_set[psa_data_set$P_ID == demoPid,]
  print(dataset$progressed[1])
  p = ggplot(data=dataset, aes(x=visitTimeYears, y=logpsa1)) + geom_line() + geom_point() + ggtitle(demoPid)
  print(p)
}

demoPatientPID = c(911, 1573, 2340, 3225, 3174)

maxPossibleFailureTime = 20
for(patientId in demoPatientPID){
  print(paste("Patient:", patientId))
  
  demoPatient = psa_data_set[psa_data_set$P_ID == patientId,]

  minVisits = 1
  plotList = vector("list", nrow(demoPatient) - minVisits)
  for(j in minVisits:nrow(demoPatient)){
    subDataSet = demoPatient[1:j, ]
    
    subDataSet$expectedFailureTime = NA
    subDataSet$survTime85 = NA
    subDataSet$survProbYouden = NA
    subDataSet$survProbAccuracy = NA
    subDataSet$survProbF1Score = NA
    subDataSet$survProbMaxTPR = NA
    
    lastVisitTime = tail(subDataSet$visitTimeYears, 1)
    lastBiopsyTime = getLastBiopsyTime(patientId, upperLimitTime = lastVisitTime)
    
    subDataSet$lastBiopsyTime = lastBiopsyTime
    
    print(paste("Last biopsy time:", lastBiopsyTime))
    survTimes = seq(lastBiopsyTime, maxPossibleFailureTime, 0.1)
    survProbs = c(1,survfitJM(joint_psa_replaced, subDataSet,
                             idVar="P_ID", last.time = lastBiopsyTime,
                             survTimes = survTimes)$summaries[[1]][, "Median"])

    nearest_time_index = which(abs(dynamicCutoffTimes_PRIAS-lastVisitTime)==min(abs(dynamicCutoffTimes_PRIAS-lastVisitTime)))[1]
     
    survProbYouden = cutoffValues_PRIAS[[nearest_time_index]]["youden"]
    survProbAccuracy = cutoffValues_PRIAS[[nearest_time_index]]["accuracy"]
    survProbF1Score = cutoffValues_PRIAS[[nearest_time_index]]["f1score"]
    survProbMaxTPR = cutoffValues_PRIAS[[nearest_time_index]]["maxTPR"]
    survProb85 = 0.85
       
    subDataSet$survTime85[j] = survTimes[which(abs(survProbs-survProb85)==min(abs(survProbs-survProb85)))[1]]
    subDataSet$survTimeYouden[j] = if(is.na(survProbYouden)){NA}else{survTimes[which(abs(survProbs-survProbYouden)==min(abs(survProbs-survProbYouden)))[1]]}
    subDataSet$survTimeAccuracy[j] = if(is.na(survProbAccuracy)){NA}else{survTimes[which(abs(survProbs-survProbAccuracy)==min(abs(survProbs-survProbAccuracy)))[1]]}
    subDataSet$survTimeF1Score[j] = if(is.na(survProbF1Score)){NA}else{survTimes[which(abs(survProbs-survProbF1Score)==min(abs(survProbs-survProbF1Score)))[1]]}
    subDataSet$survTimeMaxTPR[j] = if(is.na(survProbMaxTPR)){NA}else{survTimes[which(abs(survProbs-survProbMaxTPR)==min(abs(survProbs-survProbMaxTPR)))[1]]}
    
    subDataSet$expectedFailureTime[j]= expectedCondFailureTime(joint_psa_replaced, subDataSet, idVar = "P_ID", 
                                                                 lastBiopsyTime, maxPossibleFailureTime = maxPossibleFailureTime)
    
    print(subDataSet$expectedFailureTime[j])
    #print(sqrt(varCondFailureTime(joint_psa_replaced, subDataSet, 
    #                                 idVar = "P_ID", lastBiopsyTime, maxPossibleFailureTime = maxPossibleFailureTime)))
     
    plotList[[j-minVisits + 1]] = ggplot(data=subDataSet) + geom_point(aes(x = visitTimeYears, y=logpsa1)) +
       geom_line(aes(x = visitTimeYears, y=logpsa1)) +
       geom_vline(aes(xintercept = max(expectedFailureTime, na.rm = T), color="Expected FailureTime")) +
       geom_vline(aes(xintercept = max(lastBiopsyTime, na.rm = T), color="Last Biopsy")) +
       geom_vline(aes(xintercept = max(survTime85, na.rm = T), color="SurvTime85")) +
       geom_vline(aes(xintercept = max(survTimeYouden, na.rm = T), color="SurvTimeYouden")) +
       geom_vline(aes(xintercept = max(survTimeAccuracy, na.rm = T), color="SurvTimeAccuracy")) +
       geom_vline(aes(xintercept = max(survTimeF1Score, na.rm = T), color="SurvTimeF1Score")) +
       geom_vline(aes(xintercept = max(survTimeMaxTPR, na.rm = T), color="SurvTimeMaxTPR")) +
       ticksX(0, 20, 1) +
       #theme(legend.position="none") +
       xlab("Time(years)") + ylab("log(PSA + 1)")
     
    # print(plotList[[j-minVisits + 1]])
  }
   
  png(width=1280, height=960, filename = paste("images/prias_demo/prias_demo_pid_new", patientId, ".png", sep=""))
  multiplot(plotlist =  plotList, cols = 2)
  dev.off()
}

plotVarianceOverTime = function(modelObject, prias_P_ID){
  source("../JMBayes/Anirudh/dev/varCondFailureTime.R")
  environment(varCondFailureTime) <- asNamespace('JMbayes')
  
  ct= makeCluster(detectCores())
  registerDoParallel(ct)
  
  prias_long_i = prias_long[prias_long$P_ID==prias_P_ID,]
  
  prias_long_i$variances=foreach(k=1:nrow(prias_long_i),.combine='c', .export=c("varCondFailureTime"),
                                 .packages=c("JMbayes", "splines")) %dopar%{
  #variances=foreach(k=1:5,.combine='rbind', .packages=c("JMbayes", "splines")) %dopar%{
    subset = prias_long_i[1:k,]
    last.time = tail(subset$visitTimeYears[!is.na(subset$gleason)],1)
    last.time = 0
    return(varCondFailureTime(modelObject, subset[!is.na(subset$psa),], "P_ID", last.time, maxPossibleFailureTime = 20))
  }
  
  stopCluster(ct)
  
  biopsyTimes = prias_long_i$visitTimeYears[!is.na(prias_long_i$gleason)]
  prias_long_i$biopsyTimes = c(biopsyTimes, rep(NA, nrow(prias_long_i)-length(biopsyTimes)))
  
  pp = ggplot(data=prias_long_i) + geom_point(aes(x=visitTimeYears, y=variances)) + 
    geom_line(aes(x=visitTimeYears, y=variances)) + 
    geom_vline(aes(xintercept = biopsyTimes, na.rm=T), color="blue") + xlab("Time (years)") + 
    ylab("Variance") + ggtitle(paste("Patient ID:", prias_P_ID)) + ticksX(0, max = 20, 1) +
    ticksY(0, 100, 5)
  
  print(pp)  
}

