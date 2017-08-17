#So instead of writing the code from scratch, I will try to use the existing stuff
dynamicCutoffTimes_PRIAS = seq(1/12, max(training.prias.id.rightCens$progression_time)-0.0001, 0.1)
ct = makeCluster(detectCores())
registerDoParallel(ct)

rocList_PRIAS = foreach(tstart = dynamicCutoffTimes_PRIAS, .packages = c("splines", "JMbayes"),
                  .export=c("rocJM_mod")) %dopar%{
                    res <- tryCatch({
                      training_psa_data_set$progressed = ifelse(training_psa_data_set$progressed>0,yes = 1, no = 0) 
                      rocJM_mod(joint_psa_replaced_prias, training_psa_data_set,
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


modifyScheduledBiopsyTime = function(proposedTime, curVisitTime, lastBiopsyTime){
  if(proposedTime < curVisitTime){
    if(curVisitTime - lastBiopsyTime <= 1){
      return(lastBiopsyTime + 1)
    }else{
      return(curVisitTime)
    }
  }else{
    if(proposedTime - lastBiopsyTime <= 1){
      return(lastBiopsyTime + 1)
    }else{
      return(proposedTime)
    }
  }
}

source("../JMBayes/Anirudh/dev/expectedCondFailureTime.R")
source("../JMBayes/Anirudh/dev/varCondFailureTime.R")

#demoPatientPID = c(911, 1573, 2340, 3225, 3174)
demoPatientPID = c(2340)

maxPossibleFailureTime = 20
for(patientId in demoPatientPID){
  print(paste("Patient:", patientId))
  
  demoPatient = prias_long[prias_long$P_ID == patientId,]

  minVisits = 1
  plotList = vector("list", nrow(demoPatient) - minVisits)
  #for(j in minVisits:nrow(demoPatient)){
  for(j in c(13,18)){
    subDataSet = demoPatient[1:j, ]
    
    subDataSet$expectedFailureTime = NA
    subDataSet$survProbF1Score = NA
    
    curVisitTime = tail(subDataSet$visitTimeYears, 1)
    #lastBiopsyTime = getLastBiopsyTime(patientId, upperLimitTime = curVisitTime)
    lastBiopsyTime = tail(subDataSet$visitTimeYears[!is.na(subDataSet$gleason)],1)
    lastBiopsyTime  = 0
    
    subDataSet$lastBiopsyTime = lastBiopsyTime

    print(paste("Last biopsy time:", lastBiopsyTime))
    survTimes = seq(lastBiopsyTime, maxPossibleFailureTime, 0.1)
    survProbs = c(1,survfitJM(joint_psa_replaced_prias, subDataSet[!is.na(subDataSet$psa),],
                             idVar="P_ID", last.time = lastBiopsyTime,
                             survTimes = survTimes)$summaries[[1]][, "Median"])

    nearest_time_index = which(abs(dynamicCutoffTimes_PRIAS-curVisitTime)==min(abs(dynamicCutoffTimes_PRIAS-curVisitTime)))[1]
    survProbF1Score = cutoffValues_PRIAS[[nearest_time_index]]["f1score"]
    if(is.na(survProbF1Score)){
      subDataSet$survTimeF1Score[1] = NA
    }else{
      subDataSet$survTimeF1Score[1] = survTimes[which(abs(survProbs-survProbF1Score)==min(abs(survProbs-survProbF1Score)))[1]]
      subDataSet$survTimeF1Score[1] = modifyScheduledBiopsyTime(subDataSet$survTimeF1Score[1], curVisitTime, lastBiopsyTime)
    }
    subDataSet$expectedFailureTime[1]= expectedCondFailureTime(joint_psa_replaced_prias, subDataSet[!is.na(subDataSet$psa),], idVar = "P_ID", 
                                                                 lastBiopsyTime, maxPossibleFailureTime = maxPossibleFailureTime)
    subDataSet$expectedFailureTime[1] = modifyScheduledBiopsyTime(subDataSet$expectedFailureTime[1], curVisitTime, lastBiopsyTime)
    
    print(survProbF1Score)
    print(subDataSet$expectedFailureTime[1])
    print(subDataSet$survTimeF1Score[1])
    #print(sqrt(varCondFailureTime(joint_psa_replaced, subDataSet, 
    #                                 idVar = "P_ID", lastBiopsyTime, maxPossibleFailureTime = maxPossibleFailureTime)))
     
    plotList[[j-minVisits + 1]] = ggplot(data=subDataSet[!is.na(subDataSet$psa),]) + geom_point(aes(x = visitTimeYears, y=psa)) +
       geom_line(aes(x = visitTimeYears, y=psa)) +
      #geom_vline(aes(xintercept = max(lastBiopsyTime, na.rm = T), color="Last Biopsy", linetype="solid")) + 
      #geom_vline(aes(xintercept = max(expectedFailureTime, na.rm = T), color="Expected GR Time", linetype="dashed")) +
      # geom_vline(aes(xintercept = max(survTimeF1Score, na.rm = T), color="Dynamic risk of GR (F1)", linetype="dotted")) +
      geom_vline(aes(xintercept = max(expectedFailureTime, na.rm = T),  linetype="longdash")) +
      geom_vline(aes(xintercept = max(survTimeF1Score, na.rm = T), linetype="dotdash")) +
      geom_vline(aes(xintercept = max(lastBiopsyTime, na.rm = T),linetype="solid")) + 
       ticksX(0, 20, 2) + ylim(6, 22.5) + scale_linetype_identity(guide="legend") + 
       theme(text = element_text(size=14), axis.text=element_text(size=15), legend.position="none")+
       xlab("Time(years)") + ylab("PSA (ng/mL)")
     
    # print(plotList[[j-minVisits + 1]])
  }
  
  #png(width=1280, height=960, filename = paste("images/prias_demo/prias_demo_pid_new", patientId, ".png", sep=""))
  #plotList = lapply(plotList, function(x){x + theme(legend.position="none")})
  multiplot(plotList[[13]], plotList[[18]], cols = 2)
  
  #dev.off()
}

plotVarianceOverTime = function(modelObject, prias_P_ID, sd=T){
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
    return(varCondFailureTime(modelObject, subset[!is.na(subset$psa),], "P_ID", last.time, maxPossibleFailureTime = 20))
  }
  
  stopCluster(ct)
  
  biopsyTimes = prias_long_i$visitTimeYears[!is.na(prias_long_i$gleason)]
  prias_long_i$biopsyTimes = c(biopsyTimes, rep(NA, nrow(prias_long_i)-length(biopsyTimes)))
  
  pp = ggplot(data=prias_long_i) + geom_point(aes(x=visitTimeYears, y=sqrt(variances))) + 
    geom_line(aes(x=visitTimeYears, y=sqrt(variances))) + 
    geom_vline(aes(xintercept = biopsyTimes, na.rm=T), color="blue", linetype="dashed") + xlab("Time (years)") + 
    ylab(TeX('$SD\\left[T^*_j\\right]$')) + ggtitle(paste("Patient ID:", 3174)) + ticksX(0, max = 20, 1) +
    ticksY(0, 10, 1) + theme(text = element_text(size=14), axis.text=element_text(size=15), plot.title = element_text(hjust = 0.5))
  
  print(pp)  
}

