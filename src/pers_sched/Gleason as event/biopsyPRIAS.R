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
print(paste("Patient:", demoPatientPID))

demoPatient = prias_long[prias_long$P_ID == demoPatientPID,]

minVisits = 1
plotList = vector("list", nrow(demoPatient) - minVisits)
#for(j in minVisits:nrow(demoPatient)){
for(j in c(12,18)){ # for 2340
#for(j in c(12,15)){ # for 911
#for(j in c(3,8)){ # for 3174
  subDataSet = demoPatient[1:j, ]
  
  subDataSet$expectedFailureTime = NA
  subDataSet$survProbF1Score = NA
  
  curVisitTime = tail(subDataSet$visitTimeYears, 1)
  #lastBiopsyTime = getLastBiopsyTime(demoPatientPID, upperLimitTime = curVisitTime)
  lastBiopsyTime = tail(subDataSet$visitTimeYears[!is.na(subDataSet$gleason)],1)
  #lastBiopsyTime  = 1.03561643835616
  
  subDataSet$lastBiopsyTime = lastBiopsyTime
  
  print(paste("Last biopsy time:", lastBiopsyTime))
  survTimes = seq(lastBiopsyTime, maxPossibleFailureTime, 0.1)
  survFitObj = survfitJM(mvjoint_log2psa_plus1_spline_pt1pt54_pt1, subDataSet[!is.na(subDataSet$psa),],
                         idVar="P_ID", last.time = lastBiopsyTime,
                         survTimes = survTimes)
  survProbs = c(1,survFitObj$summaries[[1]][, "Median"])
  
  nearest_time_index = which(abs(dynamicCutoffTimes_PRIAS-lastBiopsyTime)==min(abs(dynamicCutoffTimes_PRIAS-lastBiopsyTime)))[1]
  survProbF1Score = cutoffValues_PRIAS[[nearest_time_index]]["f1score"]
  
  if(is.na(survProbF1Score)){
    subDataSet$survTimeF1Score[1] = NA
  }else{
    subDataSet$survTimeF1Score[1] = survTimes[which(abs(survProbs-survProbF1Score)==min(abs(survProbs-survProbF1Score)))[1]]
    
    #subDataSet$survTimeF1Score[1] = expectedCondFailureTime(joint_psa_replaced_prias_norm, subDataSet[!is.na(subDataSet$psa),], idVar = "P_ID", 
    #                                                        lastBiopsyTime, maxPossibleFailureTime = maxPossibleFailureTime)
    subDataSet$survTimeF1Score[1] = modifyScheduledBiopsyTime(subDataSet$survTimeF1Score[1], curVisitTime, lastBiopsyTime)
  }
  subDataSet$expectedFailureTime[1]= expectedCondFailureTime(mvjoint_log2psa_plus1_spline_pt1pt54_pt1, subDataSet[!is.na(subDataSet$psa),], idVar = "P_ID", 
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
  
  fitted.psa = survFitObj$fitted.y[[1]][[1]]
  
  subDataSet$rowType=rep("Observed", nrow(subDataSet))
  newRows = subDataSet[rep(1, length(fitted.psa)),]
  newRows$rowType="Fitted"
  newRows$visitTimeYears = survFitObj$fitted.times[[1]]
  newRows$log2psa_plus1 = fitted.psa
  
  subDataSet = rbind(subDataSet, newRows)
  
  plotList[[j-minVisits + 1]] = ggplot() + 
    #geom_point(aes(x = visitTimeYears, y=psa)) +
    geom_line(data=subDataSet[!is.na(subDataSet$psa),], mapping=aes(x = visitTimeYears, y=log2psa_plus1, linetype=rowType)) +
    geom_point(size=2,data=subDataSet[!is.na(subDataSet$psa) & subDataSet$rowType=="Observed",], mapping=aes(x = visitTimeYears, y=log2psa_plus1, linetype=rowType)) +
    geom_vline(data=subDataSet[!is.na(subDataSet$psa),], mapping=aes(xintercept = max(expectedFailureTime, na.rm = T),  linetype="Exp. GR Time"), show.legend=F) +
    geom_vline(data=subDataSet[!is.na(subDataSet$psa),], mapping=aes(xintercept = max(survTimeF1Score, na.rm = T), linetype="Dyn. risk GR"), show.legend=F) +
    geom_vline(data=subDataSet[!is.na(subDataSet$psa),], mapping=aes(xintercept = max(lastBiopsyTime, na.rm = T),linetype="Latest biopsy"), show.legend=F) +
    #geom_vline(aes(xintercept = max(expectedFailureTime, na.rm = T),  linetype="Exp. GR Time (t-distributed, df=3 errors)")) +
    #geom_vline(aes(xintercept = max(survTimeF1Score, na.rm = T), linetype="Exp. GR Time (Normally distributed errors)")) +
    ticksX(0, 20, breaks) + ylim(0, 6) + 
    scale_linetype_manual(values=c("dotted", "dashed", "longdash", "solid", "solid")) +
    theme(text = element_text(size=11), axis.text=element_text(size=11), 
          legend.text = element_text(size=11),
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size=13))+
    xlab("Time (Years)") + ylab(expression('log'[2]*' (PSA + 1)'))
  #ylab("PSA (ng/mL)")
  
  # print(plotList[[j-minVisits + 1]])
}

#png(width=1280, height=960, filename = paste("images/prias_demo/prias_demo_pid_new", patientId, ".png", sep=""))
#plotList = lapply(plotList, function(x){x + theme(legend.position="none")})
#multiplot(plotList[[13]], plotList[[18]], cols = 2)

#dev.off()


ggsave(file="report/pers_schedule/biometrics_submission/images/model_fit/demo_expfail_norm_t3.eps", 
       grid_arrange_shared_legend(plotlist = temp, ncol=3, nrow=3),
       width=8.27, height=8.27)

ggsave(file="report/pers_schedule/biometrics_submission/images/prias_demo/case_2340.eps", 
       grid_arrange_shared_legend(plotList[[3]], plotList[[18]], nrow = 1, ncol = 2),
       width=8.27, height=8.27/2)

plotVarianceOverTime = function(modelObject, prias_P_ID, sd=T){
  library(latex2exp)
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
  
  pp = ggplot(data=prias_long_i) + 
    geom_line(aes(x=visitTimeYears, y=sqrt(variances), linetype="Fitted")) + 
    geom_vline(aes(xintercept = 0, linetype="Exp. GR Time"), show.legend = F) +
    geom_vline(aes(xintercept = 0, linetype="Dyn. risk GR"), show.legend = F) +
    geom_vline(aes(xintercept = biopsyTimes, na.rm=T, linetype="Latest biopsy"), show.legend = F) + 
    xlab("Time (Years)") + 
    scale_linetype_manual(values=c("dotted", "twodash", "dashed", "solid"))  +
    ylab(TeX('$SD_g\\left[T^*_j\\right]$')) + ticksX(0, max = 20, 1) +
    ticksY(0, 10, 1) + ylim(0,7) +
    theme(text = element_text(size=11), axis.text=element_text(size=11), 
          legend.title = element_blank(), legend.text = element_text(size=11),
          plot.title = element_text(hjust = 0.5)) 
  
  print(pp)
  return(pp) 
}
ggsave(file="report/pers_schedule/biometrics_submission/images/prias_demo/case_911_t3.eps", 
       grid_arrange_shared_legend(plotList[[12]],plotList[[15]],pp, nrow = 2, ncol = 2),
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
  
  pp = ggplot(data=plotDf) +
    geom_line(aes(x=visitTimeYears, y=sd, linetype=type)) + 
    geom_vline(aes(xintercept = biopsyTimes, na.rm=T, linetype="Latest biopsy"), show.legend = F) + 
    xlab("Time (years)") + ylab(TeX('$SD_g\\left[T^*_j\\right]$')) + ticksX(0, max = 20, 1) +
    ticksY(0, 10, 1) + 
    scale_linetype_manual(values=c("solid", "dotted", "dashed"))  +
    theme(text = element_text(size=11), axis.text=element_text(size=11), 
          legend.title = element_blank(), legend.text = element_text(size=11),
          plot.title = element_text(hjust = 0.5, size=11)) 
  
  print(pp)
  return(pp) 
}

variance_3174 = compareVarianceOverTime(joint_psa_replaced_prias_norm, joint_psa_replaced_prias_t3, 3174) + ggtitle("Demo patient 2")
ggsave(file="report/pers_schedule/biometrics_submission/images/model_fit/variance_demo_patients_norm_t3.eps", 
       grid_arrange_shared_legend(variance_911, variance_3174, variance_2340, ncol=2, nrow=2),
       width=8.27, height=8.27 * 2.45/3)

set.seed(2017)
pidSample = as.numeric(as.character(training.prias.id$P_ID[sample(1:nrow(training.prias.id), size = 30)]))
pidSample = sapply(pidSample, function(pid){
  print(pid)
  if(nrow(training_psa_data_set[training_psa_data_set$P_ID==pid,])<=3){
    NA
  }else{
    pid
  }
})
pidSample = pidSample[!is.na(pidSample)]

plots = lapply(pidSample, plotLog2PSAJMFitT3VsNorm, 
               jmfitNorm=joint_psa_replaced_prias_norm, jmfitT3=joint_psa_replaced_prias_t3,
               modelNames=c("Fitted", "Fitted"))

ggsave(file="report/pers_schedule/biometrics_submission/images/model_fit/subject_fittedVsObserved_psa_norm_t3.eps", 
       grid_arrange_shared_legend(plotlist = plots[1:9], ncol=3, nrow=3),
       width=8.27, height=8.27)

#Demo patients, dynamic fit of the profile
demoPlots = vector("list", 9)
modelNames = c("Fitted", "Fitted")
modelNames=c("Fitted (Normally distributed errors)","Fitted (t-distributed df=3, errors)")
demoPlots[[1]] = plotLog2PSAJMFitT3VsNorm(jmfitNorm = joint_psa_replaced_prias_norm,
                                          jmfitT3 = joint_psa_replaced_prias_t3,
                                          patientId = 911, maxRowNum = 12,
                                          modelNames = modelNames) + 
  ggtitle("Demo patient 1")
demoPlots[[2]] = plotLog2PSAJMFitT3VsNorm(jmfitNorm = joint_psa_replaced_prias_norm,
                                          jmfitT3 = joint_psa_replaced_prias_t3,
                                          patientId = 911, maxRowNum = 15,
                                          modelNames = modelNames) + 
  ggtitle("Demo patient 1")
demoPlots[[3]] = plotLog2PSAJMFitT3VsNorm(jmfitNorm = joint_psa_replaced_prias_norm,
                                          jmfitT3 = joint_psa_replaced_prias_t3,
                                          patientId = 911, maxRowNum = NA,
                                          modelNames = modelNames) + 
  ggtitle("Demo patient 1")


demoPlots[[4]] = plotLog2PSAJMFitT3VsNorm(jmfitNorm = joint_psa_replaced_prias_norm,
                                          jmfitT3 = joint_psa_replaced_prias_t3,
                                          patientId = 3174, maxRowNum = 3, 
                                          modelNames = modelNames) + 
  ggtitle("Demo patient 2")
demoPlots[[5]] = plotLog2PSAJMFitT3VsNorm(jmfitNorm = joint_psa_replaced_prias_norm,
                                          jmfitT3 = joint_psa_replaced_prias_t3,
                                          patientId = 3174, maxRowNum = 8,
                                          modelNames = modelNames) + 
  ggtitle("Demo patient 2")
demoPlots[[6]] = plotLog2PSAJMFitT3VsNorm(jmfitNorm = joint_psa_replaced_prias_norm,
                                          jmfitT3 = joint_psa_replaced_prias_t3,
                                          patientId = 3174, maxRowNum = NA,
                                          modelNames = modelNames) + 
  ggtitle("Demo patient 2")

demoPlots[[7]] = plotLog2PSAJMFitT3VsNorm(jmfitNorm = joint_psa_replaced_prias_norm,
                                          jmfitT3 = joint_psa_replaced_prias_t3,
                                          patientId = 2340, maxRowNum = 12,
                                          modelNames = modelNames) + 
  ggtitle("Demo patient 3")
demoPlots[[8]] = plotLog2PSAJMFitT3VsNorm(jmfitNorm = joint_psa_replaced_prias_norm,
                                          jmfitT3 = joint_psa_replaced_prias_t3,
                                          patientId = 2340, maxRowNum = 18,
                                          modelNames = modelNames) + 
  ggtitle("Demo patient 3")

demoPlots[[9]] = plotLog2PSAJMFitT3VsNorm(jmfitNorm = joint_psa_replaced_prias_norm,
                                          jmfitT3 = joint_psa_replaced_prias_t3,
                                          patientId = 2340, maxRowNum = NA,
                                          modelNames = modelNames) + 
  ggtitle("Demo patient 3")

ggsave(file="report/pers_schedule/biometrics_submission/images/model_fit/fitted_demo_patients_norm_t3.eps", 
       grid_arrange_shared_legend(plotlist = demoPlots, ncol=3, nrow=3),
       width=8.27, height=8.27)

ggsave(file="report/pers_schedule/biometrics_submission/images/model_fit/fitted_demo_patients_t3.eps", 
       grid_arrange_shared_legend(plotlist =  demoPlots[c(3,6,9)], ncol=3, nrow=1),
       width=8.27, height=8.27/2.5)

