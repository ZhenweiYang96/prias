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
  for(j in c(3,8)){ # for 3174
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
    #subDataSet$expectedFailureTime[1]= expectedCondFailureTime(joint_psa_replaced_prias_t3, subDataSet[!is.na(subDataSet$psa),], idVar = "P_ID", 
    #                                                             lastBiopsyTime, maxPossibleFailureTime = maxPossibleFailureTime)
    #subDataSet$expectedFailureTime[1] = modifyScheduledBiopsyTime(subDataSet$expectedFailureTime[1], curVisitTime, lastBiopsyTime)
    
    survProbs = c(1,survfitJM(joint_psa_replaced_prias_t3_2, subDataSet[!is.na(subDataSet$psa),],
                              idVar="P_ID", last.time = lastBiopsyTime,
                              survTimes = survTimes)$summaries[[1]][, "Median"])
    
    nearest_time_index = which(abs(dynamicCutoffTimes_PRIAS-lastBiopsyTime)==min(abs(dynamicCutoffTimes_PRIAS-lastBiopsyTime)))[1]
    survProbF1Score = cutoffValues_PRIAS[[nearest_time_index]]["f1score"]
    
    if(is.na(survProbF1Score)){
      subDataSet$expectedFailureTime[1] = NA
    }else{
      subDataSet$expectedFailureTime[1] = survTimes[which(abs(survProbs-survProbF1Score)==min(abs(survProbs-survProbF1Score)))[1]]
      subDataSet$expectedFailureTime[1] = expectedCondFailureTime(joint_psa_replaced_prias_t3_2, subDataSet[!is.na(subDataSet$psa),], idVar = "P_ID", 
                                                              lastBiopsyTime, maxPossibleFailureTime = maxPossibleFailureTime)
      subDataSet$expectedFailureTime[1] = modifyScheduledBiopsyTime(subDataSet$expectedFailureTime[1], curVisitTime, lastBiopsyTime)
    }
    
    
    print(survProbF1Score)
    print(subDataSet$survTimeF1Score[1])
    print(subDataSet$expectedFailureTime[1])
    #print(sqrt(varCondFailureTime(joint_psa_replaced_prias, subDataSet[!is.na(subDataSet$log2psa),], idVar = "P_ID", lastBiopsyTime, maxPossibleFailureTime = maxPossibleFailureTime)))
    
    plotList[[j-minVisits + 1]] = ggplot(data=subDataSet[!is.na(subDataSet$psa),]) + geom_point(aes(x = visitTimeYears, y=psa)) +
      geom_line(aes(x = visitTimeYears, y=psa)) +
      # geom_vline(aes(xintercept = max(expectedFailureTime, na.rm = T),  linetype="Exp. GR Time")) +
      # geom_vline(aes(xintercept = max(survTimeF1Score, na.rm = T), linetype="Dyn. risk GR")) +
      geom_vline(aes(xintercept = max(expectedFailureTime, na.rm = T),  linetype="Exp. GR Time (t-distributed, df=3 errors)")) +
      geom_vline(aes(xintercept = max(survTimeF1Score, na.rm = T), linetype="Exp. GR Time (Normally disttributed Errors)")) +
      geom_vline(aes(xintercept = max(lastBiopsyTime, na.rm = T),linetype="Latest biopsy")) + 
      ticksX(0, 20, 1) + ylim(0, 100) + 
      scale_linetype_manual(values=c("dotted", "dashed", "solid")) +
      theme(text = element_text(size=11), axis.text=element_text(size=11), 
            legend.title = element_blank(),
            plot.title = element_text(hjust = 0.5))+
      xlab("Time(years)") + ylab("PSA (ng/mL)")
    
    # print(plotList[[j-minVisits + 1]])
  }
  
  #png(width=1280, height=960, filename = paste("images/prias_demo/prias_demo_pid_new", patientId, ".png", sep=""))
  #plotList = lapply(plotList, function(x){x + theme(legend.position="none")})
  #multiplot(plotList[[13]], plotList[[18]], cols = 2)
  
  #dev.off()
}
