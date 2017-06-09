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