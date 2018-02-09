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
demoPatientPID = c(3174)

maxPossibleFailureTime = 20
for(patientId in demoPatientPID){
  print(paste("Patient:", patientId))
  
  demoPatient = prias_long[prias_long$P_ID == patientId,]

  minVisits = 1
  plotList = vector("list", nrow(demoPatient) - minVisits)
  #for(j in minVisits:nrow(demoPatient)){
  #for(j in c(12,18)){ # for 2340
  #for(j in c(12,15)){ # for 911
  for(j in c(3,8,9)){ # for 3174
    subDataSet = demoPatient[1:j, ]
    
    subDataSet$expectedFailureTime = NA
    subDataSet$survProbF1Score = NA
    
    curVisitTime = tail(subDataSet$visitTimeYears, 1)
    #lastBiopsyTime = getLastBiopsyTime(patientId, upperLimitTime = curVisitTime)
    lastBiopsyTime = tail(subDataSet$visitTimeYears[!is.na(subDataSet$gleason)],1)
    #lastBiopsyTime  = 1.03561643835616
    
    subDataSet$lastBiopsyTime = lastBiopsyTime

    print(paste("Last biopsy time:", lastBiopsyTime))
    survTimes = seq(lastBiopsyTime, maxPossibleFailureTime, 0.1)
    survProbs = c(1,survfitJM(joint_psa_replaced_prias_norm, subDataSet[!is.na(subDataSet$psa),],
                             idVar="P_ID", last.time = lastBiopsyTime,
                             survTimes = survTimes)$summaries[[1]][, "Median"])

    nearest_time_index = which(abs(dynamicCutoffTimes_PRIAS-lastBiopsyTime)==min(abs(dynamicCutoffTimes_PRIAS-lastBiopsyTime)))[1]
    survProbF1Score = cutoffValues_PRIAS[[nearest_time_index]]["f1score"]
    
    if(is.na(survProbF1Score)){
      subDataSet$survTimeF1Score[1] = NA
    }else{
      subDataSet$survTimeF1Score[1] = survTimes[which(abs(survProbs-survProbF1Score)==min(abs(survProbs-survProbF1Score)))[1]]

      subDataSet$survTimeF1Score[1] = expectedCondFailureTime(joint_psa_replaced_prias_norm, subDataSet[!is.na(subDataSet$psa),], idVar = "P_ID", 
                                                              lastBiopsyTime, maxPossibleFailureTime = maxPossibleFailureTime)
      subDataSet$survTimeF1Score[1] = modifyScheduledBiopsyTime(subDataSet$survTimeF1Score[1], curVisitTime, lastBiopsyTime)
    }
    subDataSet$expectedFailureTime[1]= expectedCondFailureTime(joint_psa_replaced_prias_t3, subDataSet[!is.na(subDataSet$psa),], idVar = "P_ID", 
                                                                 lastBiopsyTime, maxPossibleFailureTime = maxPossibleFailureTime)
    subDataSet$expectedFailureTime[1] = modifyScheduledBiopsyTime(subDataSet$expectedFailureTime[1], curVisitTime, lastBiopsyTime)
    
    print(survProbF1Score)
    print(subDataSet$survTimeF1Score[1])
    print(subDataSet$expectedFailureTime[1])
    #print(sqrt(varCondFailureTime(joint_psa_replaced_prias, subDataSet[!is.na(subDataSet$log2psa),], idVar = "P_ID", lastBiopsyTime, maxPossibleFailureTime = maxPossibleFailureTime)))
     
    breaks=1
    if(subDataSet$expectedFailureTime[1]>8 | subDataSet$survTimeF1Score[1]>8){
      breaks = 2
    }
    
    plotList[[j-minVisits + 1]] = ggplot(data=subDataSet[!is.na(subDataSet$psa),]) + geom_point(aes(x = visitTimeYears, y=psa)) +
       geom_line(aes(x = visitTimeYears, y=psa)) +
      # geom_vline(aes(xintercept = max(expectedFailureTime, na.rm = T),  linetype="Exp. GR Time")) +
      # geom_vline(aes(xintercept = max(survTimeF1Score, na.rm = T), linetype="Dyn. risk GR")) +
      geom_vline(aes(xintercept = max(expectedFailureTime, na.rm = T),  linetype="Exp. GR Time (t-distributed, df=3 errors)")) +
      geom_vline(aes(xintercept = max(survTimeF1Score, na.rm = T), linetype="Exp. GR Time (Normally distributed errors)")) +
      geom_vline(aes(xintercept = max(lastBiopsyTime, na.rm = T),linetype="Latest biopsy")) + 
       ticksX(0, 20, breaks) + ylim(0, 100) + 
      scale_linetype_manual(values=c("dotted", "dashed", "solid")) +
       theme(text = element_text(size=11), axis.text=element_text(size=11), 
             legend.text = element_blank(size=11),
             plot.title = element_text(hjust = 0.5, size=13))+
       xlab("Time(years)") + ylab("PSA (ng/mL)")
     
    # print(plotList[[j-minVisits + 1]])
  }
  
  #png(width=1280, height=960, filename = paste("images/prias_demo/prias_demo_pid_new", patientId, ".png", sep=""))
  #plotList = lapply(plotList, function(x){x + theme(legend.position="none")})
  #multiplot(plotList[[13]], plotList[[18]], cols = 2)
  
  #dev.off()
}

ggsave(file="report/pers_schedule/biometrics_submission/images/model_fit/demo_expfail_norm_t3.eps", 
       grid_arrange_shared_legend(plotlist = temp, ncol=3, nrow=3),
       width=8.27, height=8.27)

ggsave(file="report/pers_schedule/biometrics_submission/images/prias_demo/case_2340.eps", 
       grid_arrange_shared_legend(plotList[[3]], plotList[[18]], nrow = 1, ncol = 2),
       width=8.27, height=8.27/2)

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
    geom_vline(aes(xintercept = 0, linetype="Exp. GR Time")) +
    geom_vline(aes(xintercept = 0, linetype="Dyn. risk GR")) +
    geom_vline(aes(xintercept = biopsyTimes, na.rm=T, linetype="Latest biopsy")) + xlab("Time (years)") + 
    scale_linetype_manual(values=c("dotted", "twodash", "solid"))  +
    ylab(TeX('$SD\\left[T^*_j\\right]$')) + ticksX(0, max = 20, 1) +
    ticksY(0, 10, 1) + 
    theme(text = element_text(size=11), axis.text=element_text(size=11), 
          legend.title = element_blank(), legend.text = element_text(size=11),
          plot.title = element_text(hjust = 0.5)) 
  
  print(pp)
  return(pp) 
}

ggsave(file="report/pers_schedule/biometrics_submission/images/prias_demo/case_2340.eps", 
       grid_arrange_shared_legend(p1,p2, plotList[[18]], nrow = 2, ncol = 2),
       width=8.27, height=9.69/1.25)


compareVarianceOverTime = function(jmfitNorm, jmfitT3, prias_P_ID, sd=T){
  source("../JMBayes/Anirudh/dev/varCondFailureTime.R")
  environment(varCondFailureTime) <- asNamespace('JMbayes')
  
  ct= makeCluster(detectCores())
  registerDoParallel(ct)
  
  prias_long_i = prias_long[prias_long$P_ID==prias_P_ID,]
  biopsyTimes = prias_long_i$visitTimeYears[!is.na(prias_long_i$gleason)]
  prias_long_i$biopsyTimes = c(biopsyTimes, rep(NA, nrow(prias_long_i)-length(biopsyTimes)))
  
  nrowOrig = nrow(prias_long_i)
  
  variances_norm=foreach(k=1:nrow(prias_long_i),.combine='c', .export=c("varCondFailureTime"),
                         .packages=c("JMbayes", "splines")) %dopar%{
                           #variances=foreach(k=1:5,.combine='rbind', .packages=c("JMbayes", "splines")) %dopar%{
                           subset = prias_long_i[1:k,]
                           last.time = tail(subset$visitTimeYears[!is.na(subset$gleason)],1)
                           return(varCondFailureTime(jmfitNorm, subset[!is.na(subset$psa),], "P_ID", last.time, maxPossibleFailureTime = 20))
                         }
  
  variances_t3=foreach(k=1:nrow(prias_long_i),.combine='c', .export=c("varCondFailureTime"),
                       .packages=c("JMbayes", "splines")) %dopar%{
                         #variances=foreach(k=1:5,.combine='rbind', .packages=c("JMbayes", "splines")) %dopar%{
                         subset = prias_long_i[1:k,]
                         last.time = tail(subset$visitTimeYears[!is.na(subset$gleason)],1)
                         return(varCondFailureTime(jmfitT3, subset[!is.na(subset$psa),], "P_ID", last.time, maxPossibleFailureTime = 20))
                       }
  
  stopCluster(ct)
  
  plotDf = data.frame(sd = sqrt(c(variances_norm, variances_t3)), 
                      type=rep(c("Normally distributed errors", "t-distributed (df=3) errors"),each=nrowOrig), 
                      visitTimeYears=rep(prias_long_i$visitTimeYears,2))
  plotDf$biopsyTimes = c(biopsyTimes, rep(NA, nrow(plotDf)-length(biopsyTimes)))
  
  pp = ggplot(data=plotDf) + geom_point(aes(x=visitTimeYears, y=sd)) + 
    geom_line(aes(x=visitTimeYears, y=sd, linetype=type)) + 
    scale_linetype_manual(values=c("solid", "dotted", "dashed"))  +
    ylab(TeX('$SD\\left[T^*_j\\right]$')) + ticksX(0, max = 20, 1) +
    geom_vline(aes(xintercept = biopsyTimes, na.rm=T, linetype="Latest biopsy")) + xlab("Time (years)") + 
    ticksY(0, 10, 1) + 
    theme(text = element_text(size=11), axis.text=element_text(size=11), 
          legend.title = element_blank(), legend.text = element_text(size=11),
          plot.title = element_text(hjust = 0.5, size=11)) 
  
  print(pp)
  return(pp) 
}

# variance_911 = compareVarianceOverTime(joint_psa_replaced_prias_norm, joint_psa_replaced_prias_t3, 911) + ggtitle("Demo patient 1")
ggsave(file="report/pers_schedule/biometrics_submission/images/model_fit/variance_demo_patients.eps", 
       grid_arrange_shared_legend(variance_911, variance_3174, variance_2340, ncol=2, nrow=2),
       width=8.27, height=8.27 * 2.45/3)

#Plotting to show model fit
plotLog2PSAJMFitT3VsNorm = function(jmfitNorm, jmfitT3, patientId, maxRowNum=NA){
  
  ds = prias_long[prias_long$P_ID == patientId, ]
  
  if(!is.na(maxRowNum)){
    ds=ds[1:maxRowNum,]
  }
  
  lastBiopsyTime = tail(ds$visitTimeYears[!is.na(ds$gleason)],1)
  
  predNorm = predict(jmfitNorm, ds[!is.na(ds$psa),], type="Subject", idVar="P_ID", 
                     last.time = lastBiopsyTime,
                     FtTimes = seq(min(ds$visitTimeYears), max(ds$visitTimeYears), length.out = 50))
  
  predT3 = predict(jmfitT3, ds[!is.na(ds$psa),], type="Subject", idVar="P_ID", 
                   FtTimes = seq(min(ds$visitTimeYears), max(ds$visitTimeYears), length.out = 50))
  
  ds = ds[,c("visitTimeYears", "log2psa")]
  ds$type="Observed"
  
  ds = rbind(ds, data.frame("visitTimeYears"=attributes(predNorm)$time.to.pred[[1]], 
                            "log2psa"=c(predNorm), type="Fitted (Normally distributed errors)"))
  
  ds = rbind(ds, data.frame("visitTimeYears"=attributes(predT3)$time.to.pred[[1]], 
                            "log2psa"=c(predT3), type="Fitted (t-distributed df=3, errors)"))
  
  
  breaks = if(round(max(ds$visitTimeYears))>6){
    seq(0, round(max(ds$visitTimeYears)), by=2)
  }else if(round(max(ds$visitTimeYears))>2){
    seq(0, round(max(ds$visitTimeYears)), by=1)
  }else if(round(max(ds$visitTimeYears))==2){
    seq(0, round(max(ds$visitTimeYears)), length.out=5)
  }else{
    seq(0, round(max(ds$visitTimeYears)), length.out=3)
  }
  
  p=ggplot(data=ds[!is.na(ds$log2psa),]) + 
    geom_line(aes(x=visitTimeYears, y=log2psa, linetype=type)) + 
    scale_x_continuous(breaks=breaks) + 
    geom_vline(xintercept = lastBiopsyTime, linetype="twodash") + 
    scale_linetype_manual(values=c("dotted", "dashed","solid")) + 
    theme(text = element_text(size=11), axis.text=element_text(size=11),
          legend.text = element_text(size=11), legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size=13)) +
    xlab("Time(years)") + ylab(expression('log'[2]*'(PSA)'))
  
  print(p)
  
  return(p)
}


set.seed(2017)
pidSample = training.prias.id$P_ID[sample(1:nrow(training.prias.id), size = 30)]

plots = lapply(pidSample, plotLog2PSAJMFitT3VsNorm, 
               jmfitNorm=joint_psa_replaced_prias, jmfitT3=joint_psa_replaced_prias_t3)

plots[sapply(plots, function(x){
  dat = x$data
  if(nrow(dat[dat$type=="Observed",]) <= 3){
    T
  }else{
    F
  }
})] <- NULL


ggsave(file="report/pers_schedule/biometrics_submission/images/model_fit/subject_fittedVsObserved_psa_norm_t3.eps", 
       grid_arrange_shared_legend(plotlist = plots[1:9], ncol=3, nrow=3),
       width=8.27, height=8.27)

#Demo patients, dynamic fit of the profile
demoPlots = vector("list", 6)
demoPlots[[1]] = plotLog2PSAJMFitT3VsNorm(jmfitNorm = joint_psa_replaced_prias_norm,
                                          jmfitT3 = joint_psa_replaced_prias_t3,
                                          patientId = 911, maxRowNum = 12) + 
  ggtitle("Demo patient 1")
demoPlots[[2]] = plotLog2PSAJMFitT3VsNorm(jmfitNorm = joint_psa_replaced_prias_norm,
                                          jmfitT3 = joint_psa_replaced_prias_t3,
                                          patientId = 3174, maxRowNum = 3) + 
  ggtitle("Demo patient 2")
demoPlots[[3]] = plotLog2PSAJMFitT3VsNorm(jmfitNorm = joint_psa_replaced_prias_norm,
                                          jmfitT3 = joint_psa_replaced_prias_t3,
                                          patientId = 2340, maxRowNum = 12) + 
  ggtitle("Demo patient 3")

demoPlots[[4]] = plotLog2PSAJMFitT3VsNorm(jmfitNorm = joint_psa_replaced_prias_norm,
                                          jmfitT3 = joint_psa_replaced_prias_t3,
                                          patientId = 911, maxRowNum = 15) + 
  ggtitle("Demo patient 1")
demoPlots[[5]] = plotLog2PSAJMFitT3VsNorm(jmfitNorm = joint_psa_replaced_prias_norm,
                                          jmfitT3 = joint_psa_replaced_prias_t3,
                                          patientId = 3174, maxRowNum = 8) + 
  ggtitle("Demo patient 2")
demoPlots[[6]] = plotLog2PSAJMFitT3VsNorm(jmfitNorm = joint_psa_replaced_prias_norm,
                                          jmfitT3 = joint_psa_replaced_prias_t3,
                                          patientId = 2340, maxRowNum = 18) + 
  ggtitle("Demo patient 3")

ggsave(file="report/pers_schedule/biometrics_submission/images/model_fit/fitted_demo_patients.eps", 
       grid_arrange_shared_legend(plotlist = demoPlots, ncol=3, nrow=2),
       width=8.27, height=8.27 * 2.25/3)


