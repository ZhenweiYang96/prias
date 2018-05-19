plotFittedDRE = function(modelObject, pid=NA){
  if(is.na(pid)){
    pid = sample(prias.id$P_ID, size = 1)
    print(paste("Choosing patient", pid, "randomly because pid not provided"))
  }
  
  data = modelObject$model_info$mvglmer_components$data
  data = data[!is.na(data$high_dre),] 

  rowNums = modelObject$model_info$mvglmer_components$id1
  dreFit = fitted(modelObject, process = "Longitudinal", type="Subject")[[1]]
  dreFit = dreFit[rowNums==which(prias.id$P_ID==pid)]
  dreFit = 1/(1 + exp(-(dreFit)))
    
  dreObserved = data$high_dre[data$P_ID==pid]
  
  plotdf = data.frame(time = data$visitTimeYears[data$P_ID == pid], dreObserved=dreObserved, dreFit = dreFit)
  ggplot(data=plotdf) + geom_point(aes(x=time,y=dreObserved)) + geom_line(aes(x=time, y=dreFit)) +
    ylim(c(0,1)) + xlab("Time (years)") + ylab("DRE > T1c") + ggtitle(paste("Patient", pid)) 
  
}

plotFittedPSA = function(modelObject, pid=NA){
  if(is.na(pid)){
    pid = sample(prias.id$P_ID, size = 1)
    print(paste("Choosing patient", pid, "randomly because pid not provided"))
  }
  
  data = modelObject$model_info$mvglmer_components$data
  data = data[!is.na(data$log2psa),] 
  
  rowNums = modelObject$model_info$mvglmer_components$id2
  psaFit = fitted(modelObject, process = "Longitudinal", type="Subject")[[2]]
  psaFit = psaFit[rowNums==which(prias.id$P_ID==pid)]
  
  log2PSAObserved = data$log2psa[data$P_ID==pid]
  
  plotdf = data.frame(time = data$visitTimeYears[data$P_ID == pid], log2PSAObserved=log2PSAObserved, psaFit = psaFit)
  ggplot(data=plotdf) + geom_point(aes(x=time,y=log2PSAObserved)) + geom_line(aes(x=time, y=psaFit)) +
    xlab("Time (years)") + ylab(expression('log'[2]*' (PSA)')) + ggtitle(paste("Patient", pid)) 
}

plotPSAResidualVsFitted = function(modelObject, pid=NA){
  data = modelObject$model_info$mvglmer_components$data
  data = data[!is.na(data$log2psa),] 
  
  rowNums = modelObject$model_info$mvglmer_components$id2
  psaFit = fitted(modelObject, process = "Longitudinal", type="Subject")[[2]]
  log2PSAObserved = data$log2psa
  residualPSA = log2PSAObserved - psaFit
  
  qplot(y=residualPSA, x=log2PSAObserved, xlab = "Fitted", ylab = "Residuals")
}